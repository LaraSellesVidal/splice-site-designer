"""ASCII helix diagrams for P1 and P10 helices.

Draws text-based representations of the RNA helices formed between
the intron and exon sequences, showing Watson–Crick and wobble pairs.
"""

from __future__ import annotations

from .models import LibraryMode

_WC_PAIRS = {
    ("A", "T"): True, ("T", "A"): True,
    ("G", "C"): True, ("C", "G"): True,
    ("A", "U"): True, ("U", "A"): True,
    ("G", "U"): True, ("U", "G"): True,
    ("G", "T"): True, ("T", "G"): True,  # DNA wobble equivalent
}


def _pair_symbol(nt1: str, nt2: str) -> str:
    """Return the pairing symbol between two nucleotides.

    '|' for Watson–Crick, '·' for G·U wobble, ' ' for mismatch,
    '?' for degenerate positions.
    """
    if "N" in (nt1, nt2):
        return "?"
    nt1, nt2 = nt1.upper(), nt2.upper()
    if (nt1, nt2) in (("G", "T"), ("T", "G"), ("G", "U"), ("U", "G")):
        return "·"
    if (nt1, nt2) in _WC_PAIRS:
        return "|"
    return " "


def draw_p1_helix(
    exon_segment: str,
    igs_sequence: str,
    wobble_valid: bool,
) -> str:
    """Draw the P1 helix between the 5' exon and the IGS.

    The P1 helix is antiparallel:
    - 5' exon reads 5'→3' (left to right)
    - IGS reads 3'→5' (right to left, shown bottom strand)

    The wobble G·U pair is shown at the splice site boundary.

    ::

        P1 helix:
        5' exon  5'─ ...G T C C G T ─3'  ← splice site
                          | | | | |
        IGS      3'─    C A G G C A ─5'
                                     ↑ wobble G on intron
    """
    exon = exon_segment.upper()
    igs = igs_sequence.upper()

    # The exon is 5'→3'; the IGS pairs antiparallel
    # So exon[-1] pairs with wobble G, exon[-2] pairs with IGS[0], etc.
    # Actually: IGS is the reverse complement of the exon segment,
    # so IGS[0] pairs with exon[-1], IGS[1] pairs with exon[-2], etc.

    lines: list[str] = []
    lines.append("  P1 helix:")

    # Top strand: 5' exon (5'→3')
    exon_display = " ".join(exon)
    lines.append(f"  5' exon  5'─ ...{exon_display} ─3'  splice site")

    # Pairing symbols
    pairs: list[str] = []
    for i in range(len(exon)):
        # exon[i] pairs with igs[len-1-i] (antiparallel)
        if i < len(igs):
            igs_partner = igs[len(igs) - 1 - i]
            pairs.append(_pair_symbol(exon[i], igs_partner))
        else:
            pairs.append(" ")
    pair_display = " ".join(pairs)
    lines.append(f"                  {pair_display}")

    # Bottom strand: IGS (3'→5', displayed left to right as 3'→5')
    igs_rev = igs[::-1]  # reverse to show 3'→5'
    igs_display = " ".join(igs_rev)
    wobble_note = "  G·U wobble" if wobble_valid else "  wobble disrupted"
    lines.append(f"  IGS      3'─   {igs_display} ─5'{wobble_note}")

    return "\n".join(lines)


def draw_p10_helix(
    exon_segment: str,
    p1ex_p10: str,
    stem_5prime: str,
    stem_3prime: str,
    library_mode: LibraryMode,
) -> str:
    """Draw the P10 helix and P1ex-stem overview.

    ::

        P10 helix:
        3' exon   5'─ G A C A ─3'
                      | | | |
        P1ex-P10  3'─ C T G T ─5'  (intron 5' end)

        P1ex-stem:
          5' strand: A C G T  ┐
                               ├── stem (4 bp)
          3' strand: T G C A  ┘
    """
    exon = exon_segment.upper()
    p10 = p1ex_p10.upper()

    lines: list[str] = []
    lines.append("  P10 helix:")

    # Top strand: 3' exon (5'→3')
    exon_display = " ".join(exon)
    lines.append(f"  3' exon   5'─ {exon_display} ─3'")

    # Pairing symbols (antiparallel: exon[0] pairs with p10[-1])
    pairs: list[str] = []
    for i in range(len(exon)):
        if i < len(p10):
            p10_partner = p10[len(p10) - 1 - i]
            pairs.append(_pair_symbol(exon[i], p10_partner))
        else:
            pairs.append(" ")
    pair_display = " ".join(pairs)
    lines.append(f"                {pair_display}")

    # Bottom strand: P1ex-P10 (3'→5')
    p10_rev = p10[::-1]
    p10_display = " ".join(p10_rev)
    lines.append(f"  P1ex-P10  3'─ {p10_display} ─5'  (intron 5' end)")

    # P1ex-stem
    lines.append("")
    lines.append("  P1ex-stem:")

    stem5_display = " ".join(stem_5prime.upper()) if stem_5prime else "(none)"
    stem3_display = " ".join(stem_3prime.upper()) if stem_3prime else "(none)"
    bp_count = min(len(stem_5prime), len(stem_3prime)) if stem_5prime and stem_3prime else 0

    label_5p = _stem_label(stem_5prime, library_mode, "5prime")
    label_3p = _stem_label(stem_3prime, library_mode, "3prime")

    lines.append(f"    5' strand: {stem5_display}  ┐  {label_5p}")
    lines.append(f"               {' ' * len(stem5_display)}  ├── stem ({bp_count} bp)")
    lines.append(f"    3' strand: {stem3_display}  ┘  {label_3p}")

    return "\n".join(lines)


def _stem_label(seq: str, mode: LibraryMode, strand: str) -> str:
    """Label for a stem strand depending on randomisation status."""
    if not seq:
        return ""
    if "N" in seq:
        return "← randomised"
    if mode == LibraryMode.NONE:
        return "← native"
    if strand == "5prime" and mode in (LibraryMode.FIVE_PRIME, LibraryMode.BOTH):
        return "← randomised"
    if strand == "3prime" and mode in (LibraryMode.THREE_PRIME, LibraryMode.BOTH):
        return "← randomised"
    return "← native (unchanged)"
