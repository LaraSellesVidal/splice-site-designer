"""Input validation for the Group I Intron Primer Designer."""

from __future__ import annotations

import re

from .models import LibraryMode

_DNA_RE = re.compile(r"^[ATGCatgc]+$")
_DNA_WITH_N_RE = re.compile(r"^[ATGCNatgcn]+$")


def validate_dna(seq: str, name: str, allow_n: bool = False) -> list[str]:
    """Validate that *seq* contains only DNA characters.

    Returns a list of error messages (empty if valid).
    """
    errors: list[str] = []
    if not seq:
        errors.append(f"{name}: sequence is empty.")
        return errors

    pattern = _DNA_WITH_N_RE if allow_n else _DNA_RE
    if not pattern.match(seq):
        bad = set(seq.upper()) - set("ATGCN" if allow_n else "ATGC")
        errors.append(
            f"{name}: contains non-DNA characters: {', '.join(sorted(bad))}"
        )
    return errors


def validate_inputs(
    intron_seq: str,
    five_prime_exon: str,
    three_prime_exon: str,
    library_mode: LibraryMode,
    p1_length: int | None,
    p10_length: int | None,
    p10_offset: int,
    target_tm: float,
) -> list[str]:
    """Validate all user-supplied inputs.

    Returns a list of error messages (empty if everything is valid).
    """
    errors: list[str] = []

    # --- Sequence validation ---
    errors.extend(validate_dna(intron_seq, "Intron"))
    errors.extend(validate_dna(five_prime_exon, "5' exon"))
    errors.extend(validate_dna(three_prime_exon, "3' exon"))

    # --- Length constraints ---
    if p1_length is not None:
        if not (4 <= p1_length <= 6):
            errors.append(
                f"P1 length must be 4–6 (got {p1_length})."
            )

    if p10_length is not None:
        if not (2 <= p10_length <= 5):
            errors.append(
                f"P10 length must be 2–5 (got {p10_length})."
            )

    if p10_offset < 0:
        errors.append(f"P10 offset must be ≥ 0 (got {p10_offset}).")

    # --- Exon length checks ---
    effective_p1_length = p1_length if p1_length is not None else 6
    # Need p1_length + 1 bases: p1_length for IGS pairing + 1 wobble base
    if len(five_prime_exon) < effective_p1_length + 1:
        errors.append(
            f"5' exon ({len(five_prime_exon)} nt) is shorter than "
            f"P1 length + wobble base ({effective_p1_length + 1})."
        )

    effective_p10_length = p10_length if p10_length is not None else 3
    required_3prime = p10_offset + effective_p10_length
    if len(three_prime_exon) < required_3prime:
        errors.append(
            f"3' exon ({len(three_prime_exon)} nt) is too short for "
            f"P10 offset ({p10_offset}) + P10 length ({effective_p10_length}) "
            f"= {required_3prime} nt required."
        )

    # --- Tm ---
    if not (40.0 <= target_tm <= 80.0):
        errors.append(
            f"Target Tm should be 40–80 °C (got {target_tm})."
        )

    # --- Intron minimum length ---
    if intron_seq and len(intron_seq) < 50:
        errors.append(
            f"Intron ({len(intron_seq)} nt) seems too short to be a "
            f"group I intron (typical minimum ~200 nt)."
        )

    return errors
