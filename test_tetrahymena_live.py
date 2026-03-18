#!/usr/bin/env python3
"""Live end-to-end test with the Tetrahymena thermophila rDNA group I intron.

Runs the full pipeline (Infernal analysis → design → primers) with
several exonic contexts and all four library modes.

Usage:
    python3 test_tetrahymena_live.py
"""

from __future__ import annotations

import sys
from pathlib import Path

# Add project to path
sys.path.insert(0, str(Path(__file__).parent))

from group_i_intron_designer.intron_analyzer import analyse_intron
from group_i_intron_designer.designer import IntronPrimerDesigner
from group_i_intron_designer.models import LibraryMode
from group_i_intron_designer.report import format_text_report, format_json_report
from group_i_intron_designer.sequence_utils import parse_fasta


def main():
    # --- Load the intron ---
    fasta_path = Path(__file__).parent / "Tetrahymena_intron_IC1.fasta"
    _, intron_seq = parse_fasta(fasta_path.read_text())
    print(f"Intron loaded: {len(intron_seq)} nt")
    print(f"  5' end: {intron_seq[:20]}...")
    print(f"  3' end: ...{intron_seq[-20:]}")
    print(f"  Omega G: {intron_seq[-1]}")
    print()

    # --- Analyse with Infernal ---
    print("=" * 67)
    print("  STEP 1: Infernal structural analysis")
    print("=" * 67)
    print("Running cmscan + cmalign (this may take a minute)...")
    analysis = analyse_intron(intron_seq)

    print(f"  Classification: Group I, subtype {analysis.subtype}")
    print(f"  Score: {analysis.score:.1f} bits, E-value: {analysis.e_value:.1e}")
    print(f"  IGS positions: {analysis.igs_positions}")
    print(f"  P1ex-P10 positions: {analysis.p1ex_p10_positions}")
    print(f"  P1ex-stem 5' positions: {analysis.p1ex_stem.five_prime_positions}")
    print(f"  P1ex-stem 3' positions: {analysis.p1ex_stem.three_prime_positions}")
    print(f"  Wobble G position: {analysis.wobble_g_position}")
    print(f"  Omega G position: {analysis.omega_g_position}")
    print(f"  Native IGS: {analysis.native_igs}")
    print(f"  Native P1ex-P10: {analysis.native_p1ex_p10}")
    print(f"  Native stem 5': {analysis.p1ex_stem.five_prime_native}")
    print(f"  Native stem 3': {analysis.p1ex_stem.three_prime_native}")
    print()

    # --- Test several exonic contexts ---
    test_cases = [
        {
            "name": "Native rDNA context (should produce no/minimal mutations)",
            "five_prime_exon": "CCCTTAAAAA",  # native 26S rDNA context
            "three_prime_exon": "TCTAAGTATA",
            "library": LibraryMode.NONE,
            "p10_offset": 4,
            "p10_length": 3,
        },
        {
            "name": "GFP insertion site (rational design)",
            "five_prime_exon": "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCAT",
            "three_prime_exon": "CCTGGTCGAGCTGGACGGCGACGTAAACGGC",
            "library": LibraryMode.NONE,
            "p10_offset": 4,
            "p10_length": 3,
        },
        {
            "name": "GFP insertion site (library: both strands)",
            "five_prime_exon": "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCAT",
            "three_prime_exon": "CCTGGTCGAGCTGGACGGCGACGTAAACGGC",
            "library": LibraryMode.BOTH,
            "p10_offset": 4,
            "p10_length": 3,
        },
        {
            "name": "GFP insertion site (library: 5' strand only)",
            "five_prime_exon": "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCAT",
            "three_prime_exon": "CCTGGTCGAGCTGGACGGCGACGTAAACGGC",
            "library": LibraryMode.FIVE_PRIME,
            "p10_offset": 4,
            "p10_length": 3,
        },
        {
            "name": "LacZ insertion site (library: 3' strand only)",
            "five_prime_exon": "GTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAAT",
            "three_prime_exon": "CGCTCACTGACTCGCTGCGCTCGGTCGTTCGGCTGCG",
            "library": LibraryMode.THREE_PRIME,
            "p10_offset": 4,
            "p10_length": 3,
        },
    ]

    for i, tc in enumerate(test_cases, 1):
        print()
        print("=" * 67)
        print(f"  TEST CASE {i}: {tc['name']}")
        print("=" * 67)
        print(f"  5' exon: ...{tc['five_prime_exon'][-15:]}")
        print(f"  3' exon: {tc['three_prime_exon'][:15]}...")
        print(f"  Library mode: {tc['library'].value}")
        print()

        designer = IntronPrimerDesigner(
            intron_sequence=intron_seq,
            five_prime_exon=tc["five_prime_exon"],
            three_prime_exon=tc["three_prime_exon"],
            library_mode=tc["library"],
            p10_offset=tc["p10_offset"],
            p10_length=tc["p10_length"],
            target_tm=60.0,
            analysis=analysis,
        )

        report = designer.run()
        text = format_text_report(report)
        print(text)
        print()

    # --- Also write one JSON example ---
    print("\n--- JSON output example (GFP, both strands) ---\n")
    designer = IntronPrimerDesigner(
        intron_sequence=intron_seq,
        five_prime_exon="ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCAT",
        three_prime_exon="CCTGGTCGAGCTGGACGGCGACGTAAACGGC",
        library_mode=LibraryMode.BOTH,
        p10_offset=4,
        p10_length=3,
        target_tm=60.0,
        analysis=analysis,
    )
    report = designer.run()
    print(format_json_report(report))


if __name__ == "__main__":
    main()
