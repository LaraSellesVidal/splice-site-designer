"""Parse Infernal Stockholm alignment output to identify P1/P1ex structural elements.

This module interprets the WUSS (Washington University Secondary Structure)
notation produced by ``cmalign`` to locate:

- P1ex-P10 positions (at the 5' tip of the intron)
- P1ex-stem 5' and 3' strand positions (internal stem within P1ex)
- IGS / P1 positions (pair with 5' exon)
- Wobble G position (conserved G between P1ex-stem and IGS)
- Omega G position (3' terminal G of intron)

**Key insight about CM alignments of group I introns:**

The CM models the P1/P1ex helix as ``<<<...>>>`` at the 5' end.  The
5' side (``<<<``) is inside the intron; the 3' side (``>>>``) is the
exon (absent from the intron-only alignment → gap columns).  The CM
captures the *conserved core* of P1/P1ex; the actual IGS may extend
into insert states because it varies across introns inserted in
different genes.  The parser therefore uses a hybrid CM + heuristic
strategy.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class StructuralElements:
    """Positions of key structural elements in the intron sequence.

    All positions are 0-based indices into the *ungapped* intron sequence.
    """

    p1ex_p10_positions: list[int]
    p1ex_stem_5prime_positions: list[int]
    p1ex_stem_3prime_positions: list[int]
    igs_positions: list[int]
    wobble_g_position: int
    omega_g_position: int
    p1_region_end: int  # first position of P2 (boundary of P1 region)


def _ungapped_index_map(aligned_seq: str) -> dict[int, int]:
    """Map alignment column indices to ungapped sequence positions.

    Returns a dict ``{alignment_col: seq_pos}`` for every non-gap column.
    Gap characters are ``'.'`` and ``'-'``.
    """
    mapping: dict[int, int] = {}
    seq_pos = 0
    for col, ch in enumerate(aligned_seq):
        if ch not in (".", "-"):
            mapping[col] = seq_pos
            seq_pos += 1
    return mapping


def _parse_wuss_pairs(ss_string: str) -> dict[int, int]:
    """Parse WUSS secondary structure to find base-pair partners.

    Handles nested pairs with ``()``, ``<>``, ``{}``, ``[]``.
    Returns a dict mapping each paired column to its partner column.
    """
    pairs: dict[int, int] = {}
    stacks: dict[str, list[int]] = {
        "(": [],
        "<": [],
        "[": [],
        "{": [],
    }
    close_to_open = {")": "(", ">": "<", "]": "[", "}": "{"}

    for i, ch in enumerate(ss_string):
        if ch in stacks:
            stacks[ch].append(i)
        elif ch in close_to_open:
            opener = close_to_open[ch]
            if stacks[opener]:
                j = stacks[opener].pop()
                pairs[j] = i
                pairs[i] = j
    return pairs


def _parse_pseudoknot_pairs(ss_string: str) -> dict[int, int]:
    """Parse WUSS pseudoknot annotation.

    Pseudoknot pairs are marked with uppercase letters (AA, BB, etc.)
    or sometimes lowercase (aa, bb) in WUSS.  We handle both.
    """
    pairs: dict[int, int] = {}
    letter_positions: dict[str, list[int]] = {}
    for i, ch in enumerate(ss_string):
        if ch.isalpha():
            key = ch.upper()
            if key not in letter_positions:
                letter_positions[key] = []
            letter_positions[key].append(i)

    for positions in letter_positions.values():
        n = len(positions)
        if n % 2 != 0:
            continue
        half = n // 2
        for k in range(half):
            left = positions[k]
            right = positions[n - 1 - k]
            pairs[left] = right
            pairs[right] = left

    return pairs


def parse_stockholm_alignment(
    stockholm_text: str,
    intron_seq: str,
) -> StructuralElements:
    """Extract structural elements from an Infernal cmalign Stockholm output.

    Strategy
    --------
    1. Find the ``<`` paired block at the 5' end — this is the CM's
       P1/P1ex core.  Partners are gap columns (exon, not in alignment).
    2. Find where the next helix (P2) starts — this bounds the P1 region.
    3. Within the P1 region, identify the wobble G and assign elements
       using the known group I intron 5' architecture.
    """
    aligned_seq, ss_cons = _extract_alignment_data(stockholm_text)

    col_to_seq = _ungapped_index_map(aligned_seq)
    nested_pairs = _parse_wuss_pairs(ss_cons)
    intron_len = len(intron_seq)

    sorted_cols = sorted(col_to_seq.keys())

    # --- Step 1: Find the P1/P1ex core (first < block at 5' end) ---
    p1_core_cols: list[int] = []
    for col in sorted_cols:
        ch = ss_cons[col] if col < len(ss_cons) else "."
        if ch == "<":
            p1_core_cols.append(col)
        elif ch in (",", "(", "[", "{"):
            break  # Junction or different helix — end of P1 region
        elif p1_core_cols:
            break  # Any non-< after the block ends it

    p1_core_positions = sorted(
        col_to_seq[c] for c in p1_core_cols if c in col_to_seq
    )

    # --- Step 2: Find the P2 helix start (next < or ( after P1 region) ---
    p1_region_end_pos = _find_p2_boundary(
        ss_cons, col_to_seq, sorted_cols, p1_core_cols
    )

    # --- Step 3: Find internal pairs within the P1 region (P1ex-stem) ---
    stem_5prime: list[int] = []
    stem_3prime: list[int] = []

    for col in p1_core_cols:
        if col in nested_pairs:
            partner_col = nested_pairs[col]
            if partner_col in col_to_seq:
                # Partner is also in the intron → internal pair (P1ex-stem)
                seq_pos = col_to_seq[col]
                partner_pos = col_to_seq[partner_col]
                if seq_pos < partner_pos:
                    stem_5prime.append(seq_pos)
                    stem_3prime.append(partner_pos)
                else:
                    stem_5prime.append(partner_pos)
                    stem_3prime.append(seq_pos)
            # else: partner is a gap column → pairs with exon (P1ex-P10 or IGS)

    stem_5prime = sorted(set(stem_5prime))
    stem_3prime = sorted(set(stem_3prime))

    # --- Step 4: Identify wobble G and IGS ---
    wobble_g_pos, igs_positions = _find_wobble_g_and_igs(
        intron_seq=intron_seq,
        p1_core_positions=p1_core_positions,
        p1_region_end=p1_region_end_pos,
    )

    # --- Step 5: Assign P1ex-P10 positions ---
    # P1ex-P10 = positions at the 5' tip, before the stem and wobble G
    # These are the first N positions of the P1 core that are NOT part
    # of the stem or IGS.
    all_non_p10 = set(stem_5prime) | set(stem_3prime) | set(igs_positions)
    if wobble_g_pos is not None:
        all_non_p10.add(wobble_g_pos)

    p1ex_p10_positions = sorted(
        p for p in p1_core_positions if p not in all_non_p10
    )

    # If the CM didn't resolve internal pairs (no stem), use the heuristic
    # to split pre-wobble positions into P1ex-P10 and stem strands.
    if not stem_5prime:
        p1ex_p10_positions, stem_5prime, stem_3prime, wobble_g_pos, igs_positions = (
            _heuristic_full_assignment(
                intron_seq, p1_core_positions, p1_region_end_pos
            )
        )

    omega_g_pos = intron_len - 1

    return StructuralElements(
        p1ex_p10_positions=sorted(p1ex_p10_positions),
        p1ex_stem_5prime_positions=sorted(stem_5prime),
        p1ex_stem_3prime_positions=sorted(stem_3prime),
        igs_positions=sorted(igs_positions),
        wobble_g_position=wobble_g_pos,
        omega_g_position=omega_g_pos,
        p1_region_end=p1_region_end_pos,
    )


def _extract_alignment_data(stockholm_text: str) -> tuple[str, str]:
    """Extract aligned sequence and SS_cons from Stockholm text."""
    aligned_seq = ""
    ss_cons = ""
    for line in stockholm_text.splitlines():
        line = line.strip()
        if not line or line.startswith("#=GF") or line.startswith("//"):
            continue
        if line.startswith("#=GC SS_cons"):
            ss_cons += line.split(None, 2)[2]
        elif line.startswith("#=GC") or line.startswith("#"):
            continue
        else:
            parts = line.split()
            if len(parts) == 2 and not parts[0].startswith("#"):
                aligned_seq += parts[1]

    if not aligned_seq or not ss_cons:
        raise ValueError(
            "Could not extract aligned sequence and structure from "
            "Stockholm output"
        )
    return aligned_seq, ss_cons


def _find_p2_boundary(
    ss_cons: str,
    col_to_seq: dict[int, int],
    sorted_cols: list[int],
    p1_core_cols: list[int],
) -> int:
    """Find the first position of the P2 helix (end of P1 region).

    Scans past the P1 ``<>`` block, through junction markers (``,``),
    and finds the next helix opening.
    """
    if not p1_core_cols:
        return 0

    # Find the closing > block that matches the opening < block
    max_p1_col = max(p1_core_cols)

    for col in sorted_cols:
        if col <= max_p1_col:
            continue

        ch = ss_cons[col] if col < len(ss_cons) else "."

        # Skip through unpaired, inserts, junction markers, and closing >
        if ch in (">", "_", ".", "-", "~", ":", ","):
            continue

        # First opening bracket after the P1 close → P2 start
        if ch in ("<", "(", "[", "{") and col in col_to_seq:
            return col_to_seq[col]

    # Fallback: if no P2 found, use a reasonable default
    if p1_core_cols and p1_core_cols[-1] in col_to_seq:
        last_p1_pos = col_to_seq[p1_core_cols[-1]]
        # Estimate P1 region as ~3x the P1 core length
        return min(last_p1_pos * 3 + 10, len(ss_cons))

    return 30  # conservative default


def _find_wobble_g_and_igs(
    intron_seq: str,
    p1_core_positions: list[int],
    p1_region_end: int,
) -> tuple[int, list[int]]:
    """Identify the wobble G and IGS positions within the P1 region.

    Strategy:
    1. If the last position of the P1 core block is G → that's the wobble G.
       IGS = positions immediately after it.
    2. Otherwise, scan from the end of the P1 core to P1 region end for
       the first G → wobble G candidate.  IGS follows.

    The IGS length is estimated from the P1 core size (number of CM-modeled
    paired positions).
    """
    if not p1_core_positions:
        # No CM core found — pure heuristic
        wobble = _find_wobble_g_heuristic(intron_seq)
        return wobble, []

    last_core_pos = max(p1_core_positions)

    # Estimate expected IGS length from P1 core.
    # The CM core represents the full P1/P1ex helix.  The IGS portion
    # is typically 4-6 nt.  If the core has N positions, the IGS
    # is approximately min(N, 6) nt.
    core_size = len(p1_core_positions)
    estimated_igs_len = min(core_size, 6)

    # Strategy 1: last core position is G → wobble G
    if intron_seq[last_core_pos] == "G":
        wobble_g = last_core_pos
        # IGS starts right after the wobble G
        igs_start = wobble_g + 1
        igs_end = min(igs_start + estimated_igs_len, len(intron_seq))
        igs_positions = list(range(igs_start, igs_end))
        return wobble_g, igs_positions

    # Strategy 2: scan forward from the P1 core end for a G
    scan_limit = min(last_core_pos + 20, p1_region_end, len(intron_seq))
    for pos in range(last_core_pos + 1, scan_limit):
        if intron_seq[pos] == "G":
            wobble_g = pos
            igs_start = wobble_g + 1
            igs_end = min(igs_start + estimated_igs_len, len(intron_seq))
            igs_positions = list(range(igs_start, igs_end))
            return wobble_g, igs_positions

    # Strategy 3: scan within the core for G (might be interior)
    for pos in reversed(p1_core_positions):
        if intron_seq[pos] == "G":
            wobble_g = pos
            igs_start = wobble_g + 1
            igs_end = min(igs_start + estimated_igs_len, len(intron_seq))
            igs_positions = list(range(igs_start, igs_end))
            return wobble_g, igs_positions

    # Fallback
    wobble_g = _find_wobble_g_heuristic(intron_seq)
    return wobble_g, []


def _heuristic_full_assignment(
    intron_seq: str,
    p1_core_positions: list[int],
    p1_region_end: int,
) -> tuple[list[int], list[int], list[int], int, list[int]]:
    """Full heuristic assignment of P1/P1ex elements.

    Used when the CM alignment doesn't provide enough structural detail
    (e.g., no internal pairs, no pseudoknot annotation).

    Uses the real group I intron 5' architecture::

        [stem-5'][internal loop with P10][stem-3'][wobble G][IGS][junction] → P2

    The algorithm scans for a conserved G (wobble G) and checks that the
    regions flanking the internal loop form an antiparallel stem with
    Watson–Crick complementarity.

    Returns (p1ex_p10, stem_5p, stem_3p, wobble_g, igs).
    """
    if not p1_core_positions:
        wobble_g = _find_wobble_g_heuristic(intron_seq)
        return [], [], [], wobble_g, []

    _WC = {"A": "T", "T": "A", "G": "C", "C": "G"}

    max_stem = min(len(p1_core_positions), 6)
    best: dict | None = None

    for threshold in (0.75, 0.60):
        for stem_len in range(max_stem, 1, -1):
            stem_5p_seq = intron_seq[:stem_len]

            # Wobble G must leave room for stem_3p before it and IGS after.
            min_g = 2 * stem_len          # stem_5p + stem_3p minimum
            max_g = p1_region_end - 3     # leave >= 3 nt for IGS
            max_g = min(max_g, len(intron_seq))

            for g_pos in range(min_g, max_g):
                if intron_seq[g_pos] != "G":
                    continue

                s3_start = g_pos - stem_len
                if s3_start < stem_len:
                    continue  # stem_3p would overlap stem_5p

                stem_3p_seq = intron_seq[s3_start:g_pos]

                # Antiparallel Watson–Crick complementarity check
                matches = sum(
                    1
                    for i in range(stem_len)
                    if _WC.get(stem_5p_seq[i]) == stem_3p_seq[stem_len - 1 - i]
                )

                if matches / stem_len >= threshold:
                    # Internal loop between the two stem strands
                    loop_positions = list(range(stem_len, s3_start))

                    # P10 = last 3 nt of internal loop (just before stem-3')
                    p10_len = min(3, len(loop_positions))
                    p10_pos = loop_positions[-p10_len:] if p10_len > 0 else []

                    # IGS length: typically 4-6 nt
                    available = p1_region_end - g_pos - 1
                    igs_len = min(5, max(3, available - 2))
                    igs_pos = list(range(g_pos + 1, g_pos + 1 + igs_len))

                    if best is None or stem_len > best["stem_len"]:
                        best = {
                            "stem_len": stem_len,
                            "p10": p10_pos,
                            "s5p": list(range(stem_len)),
                            "s3p": list(range(s3_start, g_pos)),
                            "g": g_pos,
                            "igs": igs_pos,
                        }

            # Already found best at this stem_len → skip shorter stems
            if best is not None and best["stem_len"] == stem_len:
                break

        if best is not None:
            break  # found at this threshold, skip relaxed threshold

    if best is None:
        wobble_g = _find_wobble_g_heuristic(intron_seq)
        return [], [], [], wobble_g, []

    return best["p10"], best["s5p"], best["s3p"], best["g"], best["igs"]


def _find_wobble_g_heuristic(intron_seq: str) -> int:
    """Find the wobble G by scanning the 5' end for a G in the expected range."""
    for i in range(3, min(20, len(intron_seq))):
        if intron_seq[i] == "G":
            return i
    return min(7, len(intron_seq) - 1)


def parse_cmscan_tblout(tblout_text: str) -> list[dict]:
    """Parse cmscan ``--tblout`` output.

    Returns a list of hit dicts sorted by E-value, each containing:
    - target_name, query_name, score, e_value, bias
    """
    hits: list[dict] = []
    for line in tblout_text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        fields = line.split()
        if len(fields) < 18:
            continue
        hits.append(
            {
                "target_name": fields[0],
                "query_name": fields[2],
                "score": float(fields[14]),
                "e_value": float(fields[15]),
                "bias": float(fields[13]),
                "seq_from": int(fields[7]),
                "seq_to": int(fields[8]),
                "mdl_from": int(fields[5]),
                "mdl_to": int(fields[6]),
            }
        )
    hits.sort(key=lambda h: h["e_value"])
    return hits
