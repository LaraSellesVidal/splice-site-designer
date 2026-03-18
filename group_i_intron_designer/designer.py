"""Core design logic for P1 and P1ex primer design.

Implements the biological design rules:
- Rule 1: P1 (IGS) design from 5' exon reverse complement with forced G at +1
- Rule 2: P1ex-P10 design from 3' exon reverse complement
- Rule 3: P1ex-stem randomisation according to library mode
- Rule 4: Conserved boundaries (wobble G, omega G)
"""

from __future__ import annotations

from .helix_drawer import draw_p1_helix, draw_p10_helix
from .models import (
    DesignReport,
    IntronAnalysis,
    LibraryInfo,
    LibraryMode,
    Mutation,
    P1Design,
    P1exDesign,
)
from .primer_designer import build_forward_primer, build_reverse_primer
from .sequence_utils import reverse_complement


class IntronPrimerDesigner:
    """Design mutagenesis primers for group I intron P1/P1ex retargeting.

    Parameters
    ----------
    intron_sequence : str
        The full intron DNA sequence.
    five_prime_exon : str
        Sequence of the 5' exon (at least p1_length nt from the splice site).
    three_prime_exon : str
        Sequence of the 3' exon (at least p10_offset + p10_length nt).
    library_mode : LibraryMode
        Which P1ex-stem strands to randomise.
    p1_length : int or None
        Number of base pairs in the P1 helix (4–6). None = use native length.
    p10_length : int or None
        Number of base pairs in P10 (2–5). None = use native length.
    p10_offset : int
        Nucleotides to skip in the 3' exon before the P10 pairing region.
    target_tm : float
        Target annealing Tm for primers (°C).
    analysis : IntronAnalysis or None
        Pre-computed intron analysis (if None, must call analyse_intron first).
    """

    def __init__(
        self,
        intron_sequence: str,
        five_prime_exon: str,
        three_prime_exon: str,
        library_mode: LibraryMode = LibraryMode.NONE,
        p1_length: int | None = None,
        p10_length: int | None = None,
        p10_offset: int = 4,
        target_tm: float = 60.0,
        analysis: IntronAnalysis | None = None,
    ):
        self.intron_seq = intron_sequence.upper()
        self.five_prime_exon = five_prime_exon.upper()
        self.three_prime_exon = three_prime_exon.upper()
        self.library_mode = library_mode
        self.p1_length = p1_length
        self.p10_length = p10_length
        self.p10_offset = p10_offset
        self.target_tm = target_tm
        self.analysis = analysis

    def design_p1(self, analysis: IntronAnalysis) -> P1Design:
        """Design the new P1 helix (IGS).

        Rule 1: Take last p1_length nt of 5' exon → reverse complement →
        force G at position +1 (the wobble position on the exon strand
        must pair with G on the intron).

        The wobble G on the intron is always maintained.  On the exon side,
        the last nt of the 5' exon ideally should be U (T in DNA) for a
        canonical G·U wobble pair.
        """
        p1_len = self.p1_length if self.p1_length is not None else analysis.native_p1_length
        if p1_len < 1:
            p1_len = analysis.native_p1_length

        # The last nt of the 5' exon is the wobble base (pairs with the
        # conserved wobble G on the intron).  The IGS pairs with the
        # p1_len exon bases *upstream* of the wobble base.
        last_exon_nt = self.five_prime_exon[-1]
        wobble_valid = last_exon_nt == "T"  # G·U wobble pair requires T(U)

        exon_segment = self.five_prime_exon[-(p1_len + 1):-1]

        # Reverse complement to get the new IGS
        new_igs = reverse_complement(exon_segment)

        helix = draw_p1_helix(exon_segment, new_igs, wobble_valid)

        return P1Design(
            new_igs=new_igs,
            exon_segment=exon_segment,
            wobble_pair=("G", last_exon_nt),
            wobble_valid=wobble_valid,
            helix_diagram=helix,
            p1_length=p1_len,
        )

    def design_p1ex(self, analysis: IntronAnalysis) -> P1exDesign:
        """Design the P1ex region (P10 pairing + stem).

        Rule 2: Skip p10_offset nt of 3' exon → take p10_length nt →
        reverse complement to get P1ex-P10.

        Rule 3: Stem strands set to N or native depending on library_mode.
        """
        p10_len = (
            self.p10_length
            if self.p10_length is not None
            else analysis.native_p1ex_p10_length
        )
        if p10_len < 1:
            p10_len = analysis.native_p1ex_p10_length

        # Extract the 3' exon segment for P10
        exon_start = self.p10_offset
        exon_end = self.p10_offset + p10_len
        exon_segment = self.three_prime_exon[exon_start:exon_end]

        # Reverse complement → new P1ex-P10
        new_p1ex_p10 = reverse_complement(exon_segment)

        # Stem sequences based on library mode
        stem = analysis.p1ex_stem
        if self.library_mode in (LibraryMode.FIVE_PRIME, LibraryMode.BOTH):
            stem_5p = "N" * len(stem.five_prime_positions)
        else:
            stem_5p = stem.five_prime_native

        if self.library_mode in (LibraryMode.THREE_PRIME, LibraryMode.BOTH):
            stem_3p = "N" * len(stem.three_prime_positions)
        else:
            stem_3p = stem.three_prime_native

        helix = draw_p10_helix(
            exon_segment, new_p1ex_p10, stem_5p, stem_3p, self.library_mode
        )

        return P1exDesign(
            new_p1ex_p10=new_p1ex_p10,
            stem_5prime_sequence=stem_5p,
            stem_3prime_sequence=stem_3p,
            library_mode=self.library_mode,
            helix_diagram=helix,
            p10_offset=self.p10_offset,
            p10_length=p10_len,
            exon_segment=exon_segment,
        )

    def compute_mutations(
        self,
        analysis: IntronAnalysis,
        p1: P1Design,
        p1ex: P1exDesign,
    ) -> list[Mutation]:
        """Identify all nucleotide changes required in the intron.

        Compares the designed sequences against the native intron to
        enumerate each position that must be mutated.
        """
        mutations: list[Mutation] = []

        # P1ex-P10 mutations
        p10_positions = sorted(analysis.p1ex_p10_positions)
        for i, pos in enumerate(p10_positions):
            if i >= len(p1ex.new_p1ex_p10):
                break
            original = self.intron_seq[pos]
            new = p1ex.new_p1ex_p10[i]
            if original != new:
                mutations.append(Mutation(pos, original, new, "P1ex-P10"))

        # P1ex-stem 5' strand
        stem = analysis.p1ex_stem
        for i, pos in enumerate(sorted(stem.five_prime_positions)):
            original = self.intron_seq[pos]
            new = p1ex.stem_5prime_sequence[i] if i < len(p1ex.stem_5prime_sequence) else original
            if original != new:
                mutations.append(
                    Mutation(pos, original, new, "P1ex-stem-5prime")
                )

        # P1ex-stem 3' strand
        for i, pos in enumerate(sorted(stem.three_prime_positions)):
            original = self.intron_seq[pos]
            new = p1ex.stem_3prime_sequence[i] if i < len(p1ex.stem_3prime_sequence) else original
            if original != new:
                mutations.append(
                    Mutation(pos, original, new, "P1ex-stem-3prime")
                )

        # IGS mutations
        igs_positions = sorted(analysis.igs_positions)
        for i, pos in enumerate(igs_positions):
            if i >= len(p1.new_igs):
                break
            original = self.intron_seq[pos]
            new = p1.new_igs[i]
            if original != new:
                mutations.append(Mutation(pos, original, new, "IGS"))

        mutations.sort(key=lambda m: m.intron_position)
        return mutations

    def compute_library_info(
        self,
        analysis: IntronAnalysis,
        p1ex: P1exDesign,
    ) -> LibraryInfo | None:
        """Compute library statistics. Returns None for rational mode."""
        if self.library_mode == LibraryMode.NONE:
            return None

        stem = analysis.p1ex_stem
        five_randomised = self.library_mode in (
            LibraryMode.FIVE_PRIME,
            LibraryMode.BOTH,
        )
        three_randomised = self.library_mode in (
            LibraryMode.THREE_PRIME,
            LibraryMode.BOTH,
        )
        n5 = len(stem.five_prime_positions) if five_randomised else 0
        n3 = len(stem.three_prime_positions) if three_randomised else 0
        total = n5 + n3
        complexity = 4**total

        warnings: list[str] = []
        if complexity > 1_000_000:
            warnings.append(
                f"Large library ({complexity:,} variants). "
                f"Ensure sufficient transformation efficiency."
            )
        if (five_randomised and not three_randomised) or (
            three_randomised and not five_randomised
        ):
            warnings.append(
                "The non-randomised strand retains native sequence. "
                "If the new stem context is very different, consider "
                "--library both."
            )

        return LibraryInfo(
            mode=self.library_mode,
            stem_5prime_randomised=five_randomised,
            stem_3prime_randomised=three_randomised,
            n_positions_5prime=n5,
            n_positions_3prime=n3,
            total_degenerate=total,
            complexity=complexity,
            warnings=warnings,
        )

    def _collect_warnings(
        self,
        analysis: IntronAnalysis,
        p1: P1Design,
        p1ex: P1exDesign,
        mutations: list[Mutation],
        fwd_primer: "Primer",  # noqa: F821
        rev_primer: "Primer",  # noqa: F821
        library: LibraryInfo | None,
    ) -> list[str]:
        """Collect all design warnings."""
        warnings: list[str] = []

        if not p1.wobble_valid:
            warnings.append(
                "G·U wobble pair disrupted: 5' exon last nt is "
                f"'{p1.wobble_pair[1]}', not T(U)."
            )

        if p1ex.p10_length > 3:
            warnings.append(
                f"Strong P10 ({p1ex.p10_length} bp) may impair exon release."
            )

        if not mutations:
            warnings.append(
                "Native intron already matches target exonic context — "
                "no mutations needed."
            )

        if fwd_primer.length > 60:
            warnings.append(
                f"Forward primer is {fwd_primer.length} nt (>60). "
                f"Consider synthesis constraints."
            )

        if fwd_primer.annealing_tm < 55.0:
            warnings.append(
                f"Low forward primer annealing Tm ({fwd_primer.annealing_tm:.1f} °C)."
            )

        if rev_primer.length < 18:
            warnings.append(
                f"Reverse primer is only {rev_primer.length} nt (<18). "
                f"Provide a longer 5' exon sequence (>=20 nt) for a proper reverse primer."
            )

        if rev_primer.annealing_tm < 55.0:
            warnings.append(
                f"Low reverse primer Tm ({rev_primer.annealing_tm:.1f} °C). "
                f"Provide a longer 5' exon sequence (>=20 nt) to improve Tm."
            )

        if analysis.score < 20.0:
            warnings.append(
                f"Low Infernal score ({analysis.score:.1f} bits). "
                f"Sequence may not be a canonical group I intron."
            )

        if library:
            warnings.extend(library.warnings)

        return warnings

    def run(self) -> DesignReport:
        """Execute the full design pipeline.

        Returns
        -------
        DesignReport
            Complete design output including primers, mutations, and
            library information.

        Raises
        ------
        ValueError
            If the intron analysis has not been provided or computed.
        """
        if self.analysis is None:
            raise ValueError(
                "IntronAnalysis must be provided.  Run analyse_intron() first."
            )

        analysis = self.analysis

        # Use native lengths if not specified
        if self.p1_length is None:
            self.p1_length = analysis.native_p1_length
        if self.p10_length is None:
            self.p10_length = analysis.native_p1ex_p10_length

        # Design P1 and P1ex
        p1 = self.design_p1(analysis)
        p1ex = self.design_p1ex(analysis)

        # Compute mutations
        mutations = self.compute_mutations(analysis, p1, p1ex)

        # Build primers
        fwd_primer = build_forward_primer(
            intron_seq=self.intron_seq,
            new_p1ex_p10=p1ex.new_p1ex_p10,
            new_igs=p1.new_igs,
            p1ex_p10_positions=analysis.p1ex_p10_positions,
            stem=analysis.p1ex_stem,
            igs_positions=analysis.igs_positions,
            wobble_g_position=analysis.wobble_g_position,
            library_mode=self.library_mode,
            target_tm=self.target_tm,
        )

        # First modified position = start of the forward primer's mutation region
        all_structural_positions = sorted(
            set(
                analysis.p1ex_p10_positions
                + analysis.p1ex_stem.five_prime_positions
                + analysis.p1ex_stem.three_prime_positions
                + [analysis.wobble_g_position]
                + analysis.igs_positions
            )
        )
        first_modified_pos = min(all_structural_positions) if all_structural_positions else 0

        rev_primer = build_reverse_primer(
            intron_seq=self.intron_seq,
            five_prime_exon=self.five_prime_exon,
            first_modified_pos=first_modified_pos,
            target_tm=self.target_tm,
        )

        # Library info
        library = self.compute_library_info(analysis, p1ex)

        # Warnings
        warnings = self._collect_warnings(
            analysis, p1, p1ex, mutations, fwd_primer, rev_primer, library
        )

        return DesignReport(
            intron_analysis=analysis,
            p1=p1,
            p1ex=p1ex,
            mutations=mutations,
            forward_primer=fwd_primer,
            reverse_primer=rev_primer,
            library=library,
            warnings=warnings,
        )
