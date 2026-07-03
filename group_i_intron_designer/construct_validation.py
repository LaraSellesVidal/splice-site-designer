"""Re-validate a retargeted intron against the covariance models.

Rebuilds the designed intron (defined mutations applied) and re-runs cmscan to
confirm it still classifies as the expected group I subtype. The IGS lies in the
P1 guide, outside the intron's conserved core, so retargeting should preserve the
catalytic architecture; this check verifies that and flags any design whose CM
match degrades. In library mode the check reflects the rational template, not
individual variants. Reuses the Infernal helpers from intron_analyzer.
"""

from __future__ import annotations

from dataclasses import dataclass

from .intron_analyzer import _check_infernal, _find_cm_database, _run_cmscan, _write_temp_fasta
from .models import Mutation
from .structural_parser import parse_cmscan_tblout


@dataclass
class ConstructValidation:
    """Result of re-scanning the retargeted construct against the CM database."""

    performed: bool
    still_group_i: bool
    subtype: str | None
    score: float | None          # bits
    e_value: float | None
    same_subtype: bool
    delta_bits: float | None     # construct minus native bits
    note: str


def apply_mutations(intron_seq: str, mutations: list[Mutation]) -> str:
    """Return the intron sequence with all DEFINED (non-'N') mutations applied."""
    seq = list(intron_seq.upper())
    for m in mutations:
        if m.new_nt.upper() != "N" and 0 <= m.intron_position < len(seq):
            seq[m.intron_position] = m.new_nt.upper()
    return "".join(seq)


def validate_construct(
    intron_seq: str,
    mutations: list[Mutation],
    native_subtype: str | None,
    native_score: float | None = None,
    library: bool = False,
) -> ConstructValidation:
    """Rebuild the retargeted intron and re-scan it against the CM database.

    Returns a ConstructValidation. If Infernal is unavailable the result has
    performed=False and the pipeline is otherwise unaffected.
    """
    try:
        cmscan_bin = _check_infernal()
        cm_dir = _find_cm_database()
    except Exception as exc:  # Infernal not on PATH
        return ConstructValidation(False, False, None, None, None, False, None,
                                   f"CM re-validation skipped: {exc}")

    designed = apply_mutations(intron_seq, mutations)
    fasta = _write_temp_fasta(designed, name="construct")
    try:
        tblout, _ = _run_cmscan(cm_dir, fasta, cmscan_bin)
    except Exception as exc:  # pragma: no cover
        return ConstructValidation(False, False, None, None, None, False, None,
                                   f"CM re-validation failed: {exc}")
    finally:
        import os
        try:
            os.unlink(fasta)
        except OSError:
            pass

    hits = parse_cmscan_tblout(tblout)
    if not hits:
        return ConstructValidation(True, False, None, None, None, False, None,
                                   "Retargeted construct no longer matches any group I CM "
                                   "(the edits may have disrupted the catalytic core).")
    best = hits[0]
    subtype = best["target_name"]
    score = best["score"]
    same = (native_subtype is not None and subtype == native_subtype)
    delta = round(score - native_score, 1) if native_score is not None else None
    lib = " (rational template; individual library variants not scanned)" if library else ""
    note = (f"Retargeted construct still classifies as {subtype} "
            f"({score:.1f} bits, E={best['e_value']:.1e})"
            + ("" if same else f"; note: native was {native_subtype}") + lib + ".")
    return ConstructValidation(True, True, subtype, score, best["e_value"], same, delta, note)
