"""Nucleotide sequence utilities.

Provides complement, reverse complement, GC content, and nearest-neighbour
Tm calculation.  All functions operate on uppercase DNA strings.
"""

from __future__ import annotations

import math

# ---------------------------------------------------------------------------
# Complement tables
# ---------------------------------------------------------------------------

_COMPLEMENT: dict[str, str] = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
    "N": "N",
}


def complement(seq: str) -> str:
    """Return the complement of a DNA sequence (5'->3' preserved)."""
    try:
        return "".join(_COMPLEMENT[nt] for nt in seq.upper())
    except KeyError as exc:
        bad = exc.args[0]
        raise ValueError(f"Non-DNA character '{bad}' in sequence") from None


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return complement(seq)[::-1]


# ---------------------------------------------------------------------------
# GC content
# ---------------------------------------------------------------------------

def gc_content(seq: str) -> float:
    """Return GC fraction (0.0–1.0) of *seq*, ignoring ambiguous bases."""
    seq = seq.upper()
    gc = sum(1 for nt in seq if nt in ("G", "C"))
    total = sum(1 for nt in seq if nt in ("A", "T", "G", "C"))
    if total == 0:
        return 0.0
    return gc / total


# ---------------------------------------------------------------------------
# Nearest-neighbour Tm (SantaLucia 1998, unified parameters)
# ---------------------------------------------------------------------------

# ΔH in kcal/mol, ΔS in cal/(mol·K)
_NN_PARAMS: dict[str, tuple[float, float]] = {
    "AA": (-7.9, -22.2),
    "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4),
    "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7),
    "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4),
    "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0),
    "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2),
    "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9),
    "CC": (-8.0, -19.9),
}

# Initiation parameters
_INIT_GC: tuple[float, float] = (0.1, -2.8)   # terminal G or C
_INIT_AT: tuple[float, float] = (2.3, 4.1)     # terminal A or T

# Gas constant in cal/(mol·K)
_R = 1.987


def tm_nearest_neighbour(
    seq: str,
    oligo_conc_nm: float = 250.0,
    na_conc_m: float = 0.05,
) -> float:
    """Nearest-neighbour Tm for a DNA oligonucleotide (SantaLucia 1998).

    Parameters
    ----------
    seq : str
        DNA sequence (only A/T/G/C; no degenerate bases).
    oligo_conc_nm : float
        Total oligonucleotide concentration in nM (default 250 nM).
    na_conc_m : float
        Monovalent cation concentration in M (default 50 mM).

    Returns
    -------
    float
        Melting temperature in °C.
    """
    seq = seq.upper()
    if len(seq) < 2:
        raise ValueError("Sequence must be at least 2 nt for Tm calculation")

    # Sum ΔH and ΔS over dinucleotide steps
    dh = 0.0
    ds = 0.0
    for i in range(len(seq) - 1):
        dinuc = seq[i : i + 2]
        if dinuc not in _NN_PARAMS:
            raise ValueError(
                f"Cannot compute Tm: non-standard dinucleotide '{dinuc}'"
            )
        h, s = _NN_PARAMS[dinuc]
        dh += h
        ds += s

    # Initiation
    for terminal in (seq[0], seq[-1]):
        if terminal in ("G", "C"):
            dh += _INIT_GC[0]
            ds += _INIT_GC[1]
        else:
            dh += _INIT_AT[0]
            ds += _INIT_AT[1]

    # Convert oligo concentration to M; self-complementary correction
    ct = oligo_conc_nm * 1e-9
    # For non-self-complementary (general case): Ct/4
    ct_term = ct / 4.0

    # Tm in Kelvin
    dh_cal = dh * 1000.0  # kcal → cal
    tm_k = dh_cal / (ds + _R * math.log(ct_term))

    # Salt correction (Owczarzy et al. simplified)
    tm_c = tm_k - 273.15
    if na_conc_m > 0:
        # SantaLucia 1998 salt correction
        tm_c += 16.6 * math.log10(na_conc_m)

    return round(tm_c, 1)


# ---------------------------------------------------------------------------
# FASTA parsing
# ---------------------------------------------------------------------------

def parse_fasta(text: str) -> tuple[str, str]:
    """Parse a single-record FASTA string.

    Returns
    -------
    tuple of (header, sequence)
        *header* includes the '>' prefix stripped; *sequence* is uppercased
        with whitespace removed.

    Raises
    ------
    ValueError
        If the text does not contain a valid FASTA record.
    """
    lines = text.strip().splitlines()
    if not lines or not lines[0].startswith(">"):
        raise ValueError("Input does not look like FASTA (no '>' header)")

    header = lines[0][1:].strip()
    seq_lines: list[str] = []
    for line in lines[1:]:
        line = line.strip()
        if line.startswith(">"):
            break  # only first record
        seq_lines.append(line)

    sequence = "".join(seq_lines).upper().replace(" ", "")
    if not sequence:
        raise ValueError("FASTA record contains no sequence data")

    return header, sequence
