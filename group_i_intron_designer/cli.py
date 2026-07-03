"""Command-line interface for the Group I Intron P1/P1ex Primer Designer."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .designer import IntronPrimerDesigner
from .intron_analyzer import analyse_intron
from .models import LibraryMode
from .report import format_json_report, format_text_report
from .sequence_utils import parse_fasta
from .validators import validate_inputs


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="group_i_intron_designer",
        description=(
            "Design site-directed mutagenesis primers to modify the P1/P1ex "
            "region of a group I intron for splicing in a new exonic context."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  # Rational design (single construct)\n"
            "  %(prog)s --intron intron.fa --five-prime-exon GTCCGT "
            "--three-prime-exon GACAGGCATG\n\n"
            "  # Library mode — both strands randomised\n"
            "  %(prog)s --intron intron.fa --five-prime-exon GTCCGT "
            "--three-prime-exon GACAGGCATG --library both\n\n"
            "  # JSON output with custom parameters\n"
            "  %(prog)s --intron intron.fa --five-prime-exon GTCCGT "
            "--three-prime-exon GACAGGCATG --library 5prime "
            "--p1-length 5 --target-tm 62 --format json\n"
        ),
    )

    parser.add_argument(
        "--intron",
        required=True,
        help="Path to FASTA file containing the intron sequence.",
    )
    parser.add_argument(
        "--five-prime-exon",
        required=True,
        help=(
            "5' exon sequence (DNA). Provide at least the last 6 nt "
            "adjacent to the splice site."
        ),
    )
    parser.add_argument(
        "--three-prime-exon",
        required=True,
        help=(
            "3' exon sequence (DNA). Provide at least the first "
            "p10_offset + p10_length nt from the splice site."
        ),
    )
    parser.add_argument(
        "--library",
        choices=["none", "5prime", "3prime", "both"],
        default="none",
        help=(
            "Library mode: 'none' for rational design (default), "
            "'5prime'/'3prime'/'both' to randomise P1ex-stem strands."
        ),
    )
    parser.add_argument(
        "--p1-length",
        type=int,
        default=None,
        help="P1 helix length in bp (4–6). Default: use native length.",
    )
    parser.add_argument(
        "--p10-length",
        type=int,
        default=None,
        help="P10 helix length in bp (2–5). Default: use native length.",
    )
    parser.add_argument(
        "--p10-offset",
        type=int,
        default=4,
        help=(
            "Nucleotides to skip in the 3' exon before the P10 pairing "
            "region (default: 4)."
        ),
    )
    parser.add_argument(
        "--target-tm",
        type=float,
        default=60.0,
        help="Target annealing Tm for primers in °C (default: 60.0).",
    )
    parser.add_argument(
        "--no-thermo",
        dest="thermo",
        action="store_false",
        help=(
            "Disable ViennaRNA thermodynamic scoring of the P1/P10 helices "
            "(on by default; skipped if ViennaRNA is absent)."
        ),
    )
    parser.add_argument(
        "--validate-construct",
        action="store_true",
        help=(
            "Rebuild the retargeted intron and re-run cmscan to confirm the "
            "subtype (requires Infernal; adds one CM scan)."
        ),
    )
    parser.add_argument(
        "--format",
        choices=["text", "json"],
        default="text",
        help="Output format: 'text' (default) or 'json'.",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output file path. Default: stdout.",
    )

    return parser


def main(argv: list[str] | None = None) -> int:
    """Entry point for the CLI.

    Returns 0 on success, 1 on validation error, 2 on runtime error.
    """
    parser = _build_parser()
    args = parser.parse_args(argv)

    # Read intron FASTA
    intron_path = Path(args.intron)
    if not intron_path.exists():
        print(f"Error: intron file not found: {intron_path}", file=sys.stderr)
        return 1

    try:
        _, intron_seq = parse_fasta(intron_path.read_text())
    except ValueError as exc:
        print(f"Error reading intron FASTA: {exc}", file=sys.stderr)
        return 1

    five_prime_exon = args.five_prime_exon.upper()
    three_prime_exon = args.three_prime_exon.upper()

    library_mode = LibraryMode(args.library)

    # Validate inputs
    errors = validate_inputs(
        intron_seq=intron_seq,
        five_prime_exon=five_prime_exon,
        three_prime_exon=three_prime_exon,
        library_mode=library_mode,
        p1_length=args.p1_length,
        p10_length=args.p10_length,
        p10_offset=args.p10_offset,
        target_tm=args.target_tm,
    )
    if errors:
        print("Validation errors:", file=sys.stderr)
        for err in errors:
            print(f"  - {err}", file=sys.stderr)
        return 1

    # Analyse intron
    try:
        print("Analysing intron structure with Infernal...", file=sys.stderr)
        analysis = analyse_intron(intron_seq)
    except (RuntimeError, FileNotFoundError, ValueError) as exc:
        print(f"Error during intron analysis: {exc}", file=sys.stderr)
        return 2

    if not analysis.is_group_i:
        print(
            "Warning: sequence may not be a group I intron.", file=sys.stderr
        )

    # Design
    designer = IntronPrimerDesigner(
        intron_sequence=intron_seq,
        five_prime_exon=five_prime_exon,
        three_prime_exon=three_prime_exon,
        library_mode=library_mode,
        p1_length=args.p1_length,
        p10_length=args.p10_length,
        p10_offset=args.p10_offset,
        target_tm=args.target_tm,
        analysis=analysis,
    )

    try:
        if args.validate_construct:
            print("Re-validating retargeted construct with cmscan...", file=sys.stderr)
        report = designer.run(
            assess_thermodynamics=args.thermo,
            validate_construct=args.validate_construct,
        )
    except Exception as exc:
        print(f"Error during design: {exc}", file=sys.stderr)
        return 2

    # Format output
    if args.format == "json":
        output_text = format_json_report(report)
    else:
        output_text = format_text_report(report)

    # Write output
    if args.output:
        Path(args.output).write_text(output_text + "\n")
        print(f"Report written to {args.output}", file=sys.stderr)
    else:
        print(output_text)

    return 0
