#!/usr/bin/env python3
"""Test the designer with a custom intron and exonic contexts."""

from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from group_i_intron_designer.intron_analyzer import analyse_intron
from group_i_intron_designer.designer import IntronPrimerDesigner
from group_i_intron_designer.models import LibraryMode
from group_i_intron_designer.report import format_text_report

INTRON_SEQ = (
    "aaatagcaatatttacctttggagggaaaagttatcaggcatgcacctggtagctagtctttaaaccaatagattgcatcggtttaaaaggcaagaccgtcaaattgcgggaaaggggtcaacagccgttcagtaccaagtctcaggggaaactttgagatggccttgcaaagggtatggtaataagctgacggacatggtcctaaccacgcagccaagtcctaagtcaacagatcttctgttgatatggatgcagttcacaGactaaatgtcggtcggggaagatgtattcttctcataagatatagtcggacctctccttaatgggagctagcggatgaagtgatgcaacactggagccgctgggaactaatttgtatgcgaaagtatattgattagttttggag"
).upper()

# 5' exon ends in T — maintains G·U wobble pair
FIVE_PRIME_EXON = "GAAGTCAAGT"
THREE_PRIME_EXON = "TTGAAGGTG"


def main():
    print(f"Intron: {len(INTRON_SEQ)} nt")
    print(f"  5' end: {INTRON_SEQ[:30]}...")
    print(f"  3' end: ...{INTRON_SEQ[-20:]}")
    print(f"  Last nt (omega G): {INTRON_SEQ[-1]}")
    print(f"  5' exon: {FIVE_PRIME_EXON}  (last nt = {FIVE_PRIME_EXON[-1]} → G·U wobble)")
    print(f"  3' exon: {THREE_PRIME_EXON}")
    print()

    # --- Analyse with Infernal ---
    print("Running Infernal analysis...")
    analysis = analyse_intron(INTRON_SEQ)

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

    # --- Rational design and library (both strands) ---
    for mode in [LibraryMode.NONE, LibraryMode.BOTH]:
        print("=" * 67)
        print(f"  Mode: {mode.value}")
        print("=" * 67)

        designer = IntronPrimerDesigner(
            intron_sequence=INTRON_SEQ,
            five_prime_exon=FIVE_PRIME_EXON,
            three_prime_exon=THREE_PRIME_EXON,
            library_mode=mode,
            p10_offset=4,
            p10_length=3,
            target_tm=60.0,
            analysis=analysis,
        )

        report = designer.run()
        print(format_text_report(report))
        print()


if __name__ == "__main__":
    main()
