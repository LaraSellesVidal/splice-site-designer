"""Primer construction for Q5 site-directed mutagenesis.

Builds forward and reverse primers incorporating the designed mutations.
The forward primer spans the P1ex and IGS regions, with degenerate (N)
positions inserted according to the library mode.  The reverse primer
anneals upstream of the mutation region.
"""

from __future__ import annotations

from .models import LibraryMode, P1exStemInfo, Primer
from .sequence_utils import gc_content, reverse_complement, tm_nearest_neighbour


def _extend_annealing_to_tm(
    seq: str,
    start: int,
    target_tm: float,
    min_length: int = 18,
    max_length: int = 35,
) -> str:
    """Extend an annealing region from *start* until Tm ≥ target_tm.

    Reads downstream (5'→3') from position *start* in *seq*.
    Returns the annealing segment.
    """
    seq_len = len(seq)
    if start >= seq_len:
        return ""

    # Start with minimum length
    end = min(start + min_length, seq_len)
    segment = seq[start:end]

    # Extend until we reach target Tm or max length
    while end < seq_len and (end - start) < max_length:
        try:
            tm = tm_nearest_neighbour(segment)
        except ValueError:
            # Sequence too short or contains non-standard bases
            end += 1
            segment = seq[start:end]
            continue

        if tm >= target_tm:
            break
        end += 1
        segment = seq[start:end]

    return segment


def _extend_annealing_upstream(
    seq: str,
    end_pos: int,
    target_tm: float,
    min_length: int = 18,
    max_length: int = 35,
) -> str:
    """Extend an annealing region upstream (towards 5') from *end_pos*.

    Reads upstream from *end_pos* (exclusive) in *seq*.
    Returns the annealing segment (5'→3' orientation).
    """
    if end_pos <= 0:
        return ""

    start = max(end_pos - min_length, 0)
    segment = seq[start:end_pos]

    while start > 0 and (end_pos - start) < max_length:
        try:
            tm = tm_nearest_neighbour(segment)
        except ValueError:
            start -= 1
            segment = seq[start:end_pos]
            continue

        if tm >= target_tm:
            break
        start -= 1
        segment = seq[start:end_pos]

    return segment


def build_forward_primer(
    intron_seq: str,
    new_p1ex_p10: str,
    new_igs: str,
    p1ex_p10_positions: list[int],
    stem: P1exStemInfo,
    igs_positions: list[int],
    wobble_g_position: int,
    library_mode: LibraryMode,
    target_tm: float,
    min_annealing: int = 18,
    max_primer_length: int = 60,
) -> Primer:
    """Build the forward mutagenesis primer.

    The primer spans (5'→3'):
        [P1ex-P10][P1ex-stem-5'][P1ex-stem-3'][G][IGS][annealing tail]
         designed   mode-dep.     mode-dep.   fix  designed  unchanged

    Parameters
    ----------
    intron_seq : str
        Full intron sequence.
    new_p1ex_p10 : str
        Designed P1ex-P10 sequence (reverse complement of 3' exon segment).
    new_igs : str
        Designed IGS sequence (reverse complement of 5' exon segment).
    p1ex_p10_positions : list[int]
        Intron positions for P1ex-P10.
    stem : P1exStemInfo
        Stem strand positions and native sequences.
    igs_positions : list[int]
        Intron positions for the IGS.
    wobble_g_position : int
        Position of the conserved wobble G.
    library_mode : LibraryMode
        Which stem strands to randomise.
    target_tm : float
        Target Tm for the annealing region.
    min_annealing : int
        Minimum annealing region length.
    max_primer_length : int
        Maximum total primer length.
    """
    # Work on a mutable copy
    intron_list = list(intron_seq)

    # Apply P1ex-P10 mutations
    sorted_p10_pos = sorted(p1ex_p10_positions)
    for i, pos in enumerate(sorted_p10_pos):
        if i < len(new_p1ex_p10):
            intron_list[pos] = new_p1ex_p10[i]

    # Apply P1ex-stem 5' strand
    sorted_stem_5p = sorted(stem.five_prime_positions)
    if library_mode in (LibraryMode.FIVE_PRIME, LibraryMode.BOTH):
        for pos in sorted_stem_5p:
            intron_list[pos] = "N"

    # Apply P1ex-stem 3' strand
    sorted_stem_3p = sorted(stem.three_prime_positions)
    if library_mode in (LibraryMode.THREE_PRIME, LibraryMode.BOTH):
        for pos in sorted_stem_3p:
            intron_list[pos] = "N"

    # Wobble G: enforce
    intron_list[wobble_g_position] = "G"

    # Apply IGS mutations
    sorted_igs_pos = sorted(igs_positions)
    for i, pos in enumerate(sorted_igs_pos):
        if i < len(new_igs):
            intron_list[pos] = new_igs[i]

    # Determine the mutation region span
    all_modified = sorted(
        set(
            sorted_p10_pos
            + sorted_stem_5p
            + sorted_stem_3p
            + [wobble_g_position]
            + sorted_igs_pos
        )
    )
    first_pos = min(all_modified)
    last_pos = max(all_modified)

    mutation_region = "".join(intron_list[first_pos : last_pos + 1])

    # Annealing tail — extends downstream of the last modified position
    annealing_start = last_pos + 1
    max_anneal = max_primer_length - len(mutation_region)
    if max_anneal < min_annealing:
        max_anneal = min_annealing

    annealing = _extend_annealing_to_tm(
        intron_seq,
        annealing_start,
        target_tm,
        min_annealing,
        max_anneal,
    )

    primer_seq = mutation_region + annealing
    n_count = primer_seq.count("N")

    # Compute Tm of the annealing region only (no degenerate bases)
    try:
        ann_tm = tm_nearest_neighbour(annealing)
    except ValueError:
        ann_tm = 0.0

    return Primer(
        name="Forward",
        sequence=primer_seq,
        length=len(primer_seq),
        mutation_region=mutation_region,
        annealing_region=annealing,
        annealing_tm=ann_tm,
        annealing_gc_percent=round(gc_content(annealing) * 100, 1),
        contains_degenerate=n_count > 0,
        degenerate_count=n_count,
    )


def build_reverse_primer(
    intron_seq: str,
    five_prime_exon: str,
    first_modified_pos: int,
    target_tm: float,
    min_length: int = 18,
    max_length: int = 35,
) -> Primer:
    """Build the reverse primer.

    The reverse primer anneals to the opposite strand, upstream of the
    forward primer's 5' end.  It never contains degenerate positions.

    For non-overlapping Q5-style mutagenesis the forward primer starts at
    intron position 0 (covering the full P1ex/IGS mutation region).  The
    reverse primer therefore anneals to the **5' exon** sequence immediately
    upstream of the intron, reading 5'→3' on the antisense strand.

    When the first mutation is further into the intron (e.g. rational mode
    where only P10/IGS change), the annealing region may include native
    intron sequence upstream of the first mutation as well.
    """
    # Build upstream context: 5' exon + native intron before first mutation
    upstream = five_prime_exon + intron_seq[:first_modified_pos]

    annealing_sense = _extend_annealing_upstream(
        upstream,
        len(upstream),
        target_tm,
        min_length,
        max_length,
    )

    if not annealing_sense:
        # Fallback: use whatever upstream context we have
        take = min(min_length, len(upstream))
        annealing_sense = upstream[-take:]

    # Reverse complement for the actual primer
    primer_seq = reverse_complement(annealing_sense)

    try:
        ann_tm = tm_nearest_neighbour(annealing_sense)
    except ValueError:
        ann_tm = 0.0

    return Primer(
        name="Reverse",
        sequence=primer_seq,
        length=len(primer_seq),
        mutation_region="",
        annealing_region=primer_seq,
        annealing_tm=ann_tm,
        annealing_gc_percent=round(gc_content(primer_seq) * 100, 1),
        contains_degenerate=False,
        degenerate_count=0,
    )
