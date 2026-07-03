"""Report generation for the Group I Intron P1/P1ex Primer Designer.

Produces both human-readable text reports and machine-readable JSON output.
"""

from __future__ import annotations

import json
from dataclasses import asdict
from typing import Any

from .models import DesignReport, LibraryMode


_HEADER_WIDTH = 67


def _header(title: str) -> str:
    """Produce a section header."""
    return f"── {title} " + "─" * max(0, _HEADER_WIDTH - len(title) - 4)


def _format_library_mode_label(mode: LibraryMode) -> str:
    """Human-readable label for the library mode."""
    labels = {
        LibraryMode.NONE: "Rational design (single construct)",
        LibraryMode.FIVE_PRIME: "Library — randomise 5' strand of P1ex-stem",
        LibraryMode.THREE_PRIME: "Library — randomise 3' strand of P1ex-stem",
        LibraryMode.BOTH: "Library — randomise BOTH strands of P1ex-stem",
    }
    return labels.get(mode, str(mode.value))


def _format_positions(positions: list[int]) -> str:
    """Format a position list as a compact range string."""
    if not positions:
        return "(none)"
    if len(positions) == 1:
        return f"pos {positions[0]}"
    return f"pos {positions[0]}–{positions[-1]}"


def format_text_report(report: DesignReport) -> str:
    """Format the design report as a human-readable text string."""
    lines: list[str] = []
    a = report.intron_analysis
    p1 = report.p1
    p1ex = report.p1ex

    # Banner
    lines.append("=" * _HEADER_WIDTH)
    lines.append(
        "  GROUP I INTRON P1/P1ex PRIMER DESIGNER — DESIGN REPORT"
    )
    lines.append("=" * _HEADER_WIDTH)
    lines.append(
        f"  Mode: {_format_library_mode_label(report.p1ex.library_mode)}"
    )
    lines.append("")

    # ── Intron Analysis ──
    lines.append(_header("Intron Analysis"))
    lines.append("")
    lines.append(
        f"  Classification:     Group I, subtype {a.subtype} "
        f"({a.score:.1f} bits, E={a.e_value:.1e})"
    )
    lines.append(f"  Intron length:      {a.intron_length} nt")
    lines.append("")
    lines.append("  Native 5' structure:")

    # P1ex-P10
    p10_pos_str = _format_positions(a.p1ex_p10_positions)
    lines.append(
        f"    P1ex-P10 ({p10_pos_str}):        "
        f"{a.native_p1ex_p10:8s}  pairs with 3' exon (P10)"
    )

    # P1ex-stem 5' strand
    stem = a.p1ex_stem
    s5_pos = _format_positions(stem.five_prime_positions)
    lines.append(
        f"    P1ex-stem 5' strand ({s5_pos}): "
        f"{stem.five_prime_native:4s}  ┐"
    )
    lines.append(
        f"    {'':38s}├── internal stem ({stem.stem_length} bp)"
    )
    s3_pos = _format_positions(stem.three_prime_positions)
    lines.append(
        f"    P1ex-stem 3' strand ({s3_pos}): "
        f"{stem.three_prime_native:4s}  ┘"
    )

    # Wobble G
    lines.append(
        f"    Wobble G (pos {a.wobble_g_position}):           G"
    )

    # IGS
    igs_pos = _format_positions(a.igs_positions)
    lines.append(
        f"    IGS / P1 ({igs_pos}):       "
        f"{a.native_igs:8s}  pairs with 5' exon (P1)"
    )
    lines.append("")

    # ── New P1 Helix ──
    lines.append(_header("New P1 Helix"))
    lines.append("")
    lines.append(
        f"  5' exon segment:  {p1.exon_segment}  "
        f"(last {p1.p1_length} nt of 5' exon)"
    )
    lines.append(f"  New IGS:          {p1.new_igs}  (reverse complement)")
    wob_status = "valid (G·U)" if p1.wobble_valid else "DISRUPTED"
    lines.append(
        f"  Wobble pair:      G · {p1.wobble_pair[1]}  ({wob_status})"
    )
    lines.append("")
    lines.append(p1.helix_diagram)
    lines.append("")

    # ── New P10 Helix / P1ex Design ──
    lines.append(_header("P1ex Design"))
    lines.append("")
    lines.append(
        f"  3' exon segment:  {p1ex.exon_segment}  "
        f"(offset {p1ex.p10_offset}, length {p1ex.p10_length})"
    )
    lines.append(
        f"  New P1ex-P10:     {p1ex.new_p1ex_p10}  (reverse complement)"
    )

    # Stem display
    s5_label = "randomised" if "N" in p1ex.stem_5prime_sequence else "native"
    s3_label = "randomised" if "N" in p1ex.stem_3prime_sequence else "native"
    lines.append(
        f"  Stem 5' strand:   {p1ex.stem_5prime_sequence}  ({s5_label})"
    )
    lines.append(
        f"  Stem 3' strand:   {p1ex.stem_3prime_sequence}  ({s3_label})"
    )
    lines.append(f"  Wobble G:         G  (fixed)")
    lines.append("")
    lines.append(p1ex.helix_diagram)
    lines.append("")

    # ── Library ──
    if report.library is not None:
        lib = report.library
        lines.append(_header("Library"))
        lines.append("")
        mode_labels = {
            LibraryMode.FIVE_PRIME: "5' strand only",
            LibraryMode.THREE_PRIME: "3' strand only",
            LibraryMode.BOTH: "Both strands",
        }
        lines.append(
            f"  Randomisation mode:     {mode_labels.get(lib.mode, str(lib.mode.value))}"
        )
        if lib.stem_5prime_randomised:
            lines.append(
                f"  5' strand positions:    {lib.n_positions_5prime}"
            )
        if lib.stem_3prime_randomised:
            lines.append(
                f"  3' strand positions:    {lib.n_positions_3prime}"
            )
        lines.append(
            f"  Total degenerate (N):   {lib.total_degenerate}"
        )
        lines.append(
            f"  Library complexity:     4^{lib.total_degenerate} = "
            f"{lib.complexity:,} variants"
        )
        lines.append("")

    # ── Mutations ──
    lines.append(_header("Mutations"))
    lines.append("")
    if report.mutations:
        for m in report.mutations:
            lines.append(
                f"  pos {m.intron_position:4d}:  "
                f"{m.original_nt} → {m.new_nt}  ({m.region})"
            )
    else:
        lines.append("  No mutations required — native matches target context.")
    lines.append("")

    # ── Primers ──
    lines.append(_header("Primers"))
    lines.append("")

    fwd = report.forward_primer
    lines.append(f"  Forward: 5'- {fwd.sequence} -3'")
    if fwd.mutation_region and fwd.annealing_region:
        mut_len = len(fwd.mutation_region)
        ann_len = len(fwd.annealing_region)
        lines.append(
            f"               {'─' * mut_len} {'─' * ann_len}"
        )
        lines.append(
            f"               ├─ modified ─┤├── annealing ─┤"
        )
    lines.append(
        f"    Length: {fwd.length} nt | "
        f"Annealing Tm: {fwd.annealing_tm:.1f} °C | "
        f"GC: {fwd.annealing_gc_percent:.0f}%"
    )
    if fwd.contains_degenerate:
        lines.append(
            f"    Degenerate positions: {fwd.degenerate_count}"
        )
    lines.append("")

    rev = report.reverse_primer
    lines.append(f"  Reverse: 5'- {rev.sequence} -3'")
    lines.append(
        f"    Length: {rev.length} nt | "
        f"Annealing Tm: {rev.annealing_tm:.1f} °C | "
        f"GC: {rev.annealing_gc_percent:.0f}%"
    )
    lines.append("")

    if fwd.contains_degenerate:
        lines.append(
            "  Note: Order forward primer with mixed bases at N positions."
        )
        lines.append("")

    # ── Thermodynamic Assessment ──
    if report.thermodynamics is not None:
        t = report.thermodynamics
        lines.append(_header("Thermodynamic Assessment (P1 / P10)"))
        lines.append("")
        if not t.viennarna_available:
            lines.append(
                "  ViennaRNA not installed; scores unavailable "
                "(pip install ViennaRNA to enable)."
            )
        else:
            def _dg(v: float | None) -> str:
                return f"{v:+.1f} kcal/mol" if v is not None else "n/a"

            def _margin(m: float | None, frac: float | None) -> str:
                if m is None:
                    return "n/a"
                return f"{m:+.1f} kcal/mol vs random (P_random≤real = {frac:.2f})"

            lines.append(
                f"  P1  (IGS·5' exon):   {_dg(t.p1_dg):>16s}   [{t.p1_stability}]"
            )
            lines.append(
                f"    exon-specificity:  {_margin(t.p1_specificity_margin, t.p1_frac_random_stronger)}"
            )
            lines.append(
                f"  P10 (P1ex·3' exon):  {_dg(t.p10_dg):>16s}   [{t.p10_stability}]"
            )
            lines.append(
                f"    exon-specificity:  {_margin(t.p10_specificity_margin, t.p10_frac_random_stronger)}"
            )
            lines.append("")
            lines.append(
                "  Note: scores rank designs; they do not predict splicing. "
                "Use the P1ex"
            )
            lines.append("  library selection as the functional test.")
        lines.append("")

    # ── Construct Re-validation ──
    if report.construct_validation is not None and report.construct_validation.performed:
        cv = report.construct_validation
        lines.append(_header("Retargeted Construct Re-validation (cmscan)"))
        lines.append("")
        status = "PASS" if cv.still_group_i else "FAIL"
        lines.append(f"  Still a group I intron:  {status}")
        if cv.subtype is not None:
            same = "same as native" if cv.same_subtype else "native was different"
            lines.append(
                f"  Best CM match:           {cv.subtype} "
                f"({cv.score:.1f} bits, E={cv.e_value:.1e}); {same}"
            )
        if cv.delta_bits is not None:
            lines.append(
                f"  Bit-score vs native:     {cv.delta_bits:+.1f} bits"
            )
        lines.append(f"  {cv.note}")
        lines.append("")

    # ── Warnings ──
    if report.warnings:
        lines.append(_header("Warnings"))
        lines.append("")
        for w in report.warnings:
            lines.append(f"  ⚠  {w}")
        lines.append("")

    lines.append("=" * _HEADER_WIDTH)
    return "\n".join(lines)


def _serialise_for_json(obj: Any) -> Any:
    """Recursively convert dataclasses and enums for JSON serialisation."""
    if hasattr(obj, "__dataclass_fields__"):
        result = {}
        for field_name in obj.__dataclass_fields__:
            result[field_name] = _serialise_for_json(getattr(obj, field_name))
        return result
    if isinstance(obj, (list, tuple)):
        return [_serialise_for_json(item) for item in obj]
    if hasattr(obj, "value"):  # Enum
        return obj.value
    return obj


def format_json_report(report: DesignReport) -> str:
    """Format the design report as a JSON string."""
    data = _serialise_for_json(report)
    return json.dumps(data, indent=2)
