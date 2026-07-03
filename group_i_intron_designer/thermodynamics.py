"""Thermodynamic scoring of designed P1/P10 helices via ViennaRNA.

The IGS (P1) and P1ex-P10 are the reverse complement of the target exon, so they
pair it by construction. This module scores two further properties:

- Stability: minimum free energy of the intended helix (GC-rich exons give strong
  P1 helices, AU-rich exons weak ones).
- Specificity: whether the guide pairs its target better than composition-matched
  shuffles of the same exon (a short wobble-permissive P1 can otherwise be
  satisfied by base composition alone).

The scores rank designs; they do not predict splicing (P1 exon-pairing does not
covary above phylogeny, so it cannot be validated from sequence). The P1ex library
selection remains the functional test. If ViennaRNA is not installed every
function returns None and the pipeline is unaffected.
"""

from __future__ import annotations

import random

try:  # ViennaRNA is an optional dependency (pip install ViennaRNA)
    import RNA  # type: ignore

    _HAVE_VIENNA = True
except Exception:  # pragma: no cover - exercised only without ViennaRNA
    _HAVE_VIENNA = False


def viennarna_available() -> bool:
    """True if ViennaRNA (import RNA) is importable."""
    return _HAVE_VIENNA


def _rna(seq: str) -> str:
    return seq.upper().replace("T", "U")


def duplex_dg(seq_a: str, seq_b: str) -> float | None:
    """Minimum free energy of the intermolecular *seq_a*:*seq_b* duplex (kcal/mol).

    Returns None if ViennaRNA is unavailable or either sequence is too short.
    """
    if not _HAVE_VIENNA or len(seq_a) < 2 or len(seq_b) < 2:
        return None
    return round(float(RNA.duplexfold(_rna(seq_a), _rna(seq_b)).energy), 2)


def classify_stability(dg: float | None) -> str:
    """Qualitative label for a helix free energy (kcal/mol)."""
    if dg is None:
        return "unknown"
    if dg <= -8.0:
        return "strong"
    if dg <= -4.0:
        return "moderate"
    return "weak"


def specificity_margin(
    exon_segment: str,
    guide_seq: str,
    n_random: int = 60,
    seed: int = 0,
) -> tuple[float | None, float | None]:
    """Composition-controlled exon-specificity of a designed guide.

    Compares the intended *exon_segment*:*guide_seq* duplex energy against
    *n_random* composition-matched shuffles of the exon paired to the same guide.
    Returns (margin_kcal, frac_random_stronger), where margin = median(random dG)
    minus real dG (positive means the real exon pairs the guide better than a
    composition-matched random sequence) and frac_random_stronger is the fraction
    of shuffles at least as stable as the real pair (small means specific).
    Returns (None, None) without ViennaRNA.
    """
    real = duplex_dg(exon_segment, guide_seq)
    if real is None:
        return None, None
    bases = list(_rna(exon_segment))
    rng = random.Random(seed)
    rand_dgs: list[float] = []
    for _ in range(n_random):
        shuffled = bases[:]
        rng.shuffle(shuffled)
        d = duplex_dg("".join(shuffled), guide_seq)
        if d is not None:
            rand_dgs.append(d)
    if not rand_dgs:
        return None, None
    rand_dgs.sort()
    median = rand_dgs[len(rand_dgs) // 2]
    frac_stronger = sum(1 for d in rand_dgs if d <= real) / len(rand_dgs)
    return round(median - real, 2), round(frac_stronger, 3)


def assess(
    p1_exon_segment: str,
    new_igs: str,
    p10_exon_segment: str,
    new_p1ex_p10: str,
) -> dict:
    """Assess a design: P1/P10 duplex free energies, stability labels,
    composition-controlled specificity margins, and any warnings.

    All numeric values are None without ViennaRNA.
    """
    p1_dg = duplex_dg(p1_exon_segment, new_igs)
    p10_dg = duplex_dg(p10_exon_segment, new_p1ex_p10)
    p1_margin, p1_frac = specificity_margin(p1_exon_segment, new_igs)
    p10_margin, p10_frac = specificity_margin(p10_exon_segment, new_p1ex_p10)

    warnings: list[str] = []
    if not _HAVE_VIENNA:
        warnings.append(
            "ViennaRNA not installed; thermodynamic scores skipped "
            "(pip install ViennaRNA to enable)."
        )
    else:
        if p1_dg is not None and p1_dg > -4.0:
            warnings.append(
                f"Weak predicted P1 helix (dG={p1_dg} kcal/mol); AU-rich 5' exon "
                "targets are harder to splice and rely more on library selection."
            )
        if p1_margin is not None and p1_margin < 2.0:
            warnings.append(
                f"Low P1 exon-specificity margin ({p1_margin} kcal/mol vs random): "
                "the guide pairs this target little better than composition-matched "
                "sequence; off-target pairing is more likely."
            )
    return {
        "viennarna": _HAVE_VIENNA,
        "p1_dg": p1_dg,
        "p1_stability": classify_stability(p1_dg),
        "p1_specificity_margin": p1_margin,
        "p1_frac_random_stronger": p1_frac,
        "p10_dg": p10_dg,
        "p10_stability": classify_stability(p10_dg),
        "p10_specificity_margin": p10_margin,
        "p10_frac_random_stronger": p10_frac,
        "warnings": warnings,
    }
