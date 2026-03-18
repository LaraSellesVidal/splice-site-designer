"""Tests for the core IntronPrimerDesigner class."""

from __future__ import annotations

import pytest

from group_i_intron_designer.designer import IntronPrimerDesigner
from group_i_intron_designer.models import (
    IntronAnalysis,
    LibraryMode,
    P1exStemInfo,
)


def _make_analysis(
    intron_seq: str = "ACGACCTGAGGTAACGT" + "A" * 383 + "G",
    igs_positions: list[int] | None = None,
    p1ex_p10_positions: list[int] | None = None,
    stem_5p_positions: list[int] | None = None,
    stem_3p_positions: list[int] | None = None,
    wobble_g: int = 7,
    omega_g: int | None = None,
) -> IntronAnalysis:
    """Helper to construct a mock IntronAnalysis."""
    if igs_positions is None:
        igs_positions = [8, 9, 10, 11, 12]
    if p1ex_p10_positions is None:
        p1ex_p10_positions = [0, 1, 2]
    if stem_5p_positions is None:
        stem_5p_positions = [3, 4]
    if stem_3p_positions is None:
        stem_3p_positions = [5, 6]
    if omega_g is None:
        omega_g = len(intron_seq) - 1

    native_igs = "".join(intron_seq[p] for p in igs_positions)
    native_p10 = "".join(intron_seq[p] for p in p1ex_p10_positions)
    stem_5p_native = "".join(intron_seq[p] for p in stem_5p_positions)
    stem_3p_native = "".join(intron_seq[p] for p in stem_3p_positions)

    return IntronAnalysis(
        is_group_i=True,
        subtype="IC1",
        score=89.3,
        e_value=2.1e-24,
        igs_positions=igs_positions,
        p1ex_p10_positions=p1ex_p10_positions,
        p1ex_stem=P1exStemInfo(
            five_prime_positions=stem_5p_positions,
            three_prime_positions=stem_3p_positions,
            five_prime_native=stem_5p_native,
            three_prime_native=stem_3p_native,
            stem_length=min(len(stem_5p_positions), len(stem_3p_positions)),
        ),
        wobble_g_position=wobble_g,
        omega_g_position=omega_g,
        native_igs=native_igs,
        native_p1ex_p10=native_p10,
        native_p1_length=len(igs_positions),
        native_p1ex_p10_length=len(p1ex_p10_positions),
        intron_length=len(intron_seq),
    )


class TestDesignP1:
    """Tests for P1 (IGS) design — Rule 1."""

    def test_basic_reverse_complement(self):
        """IGS should be the RC of the exon segment upstream of the wobble base.

        Mock analysis has native_p1_length=5 (positions 8-12), so
        designer takes 5 nt upstream of wobble base from 'GTCCGT':
        exon[-6:-1] = 'GTCCG' → RC = 'CGGAC'.
        """
        intron = "ACGACCTGAGGTAACGT" + "A" * 383 + "G"
        analysis = _make_analysis(intron)
        designer = IntronPrimerDesigner(
            intron_sequence=intron,
            five_prime_exon="GTCCGT",
            three_prime_exon="GACAGGCATG",
            analysis=analysis,
        )
        p1 = designer.design_p1(analysis)
        # native_p1_length=5, exon[-6:-1] = GTCCG → RC = CGGAC
        assert p1.new_igs == "CGGAC"
        assert p1.p1_length == 5

    def test_p1_length_matches_request(self):
        intron = "ACGACCTGAGGTAACGT" + "A" * 383 + "G"
        analysis = _make_analysis(intron)
        designer = IntronPrimerDesigner(
            intron_sequence=intron,
            five_prime_exon="AGTCCGT",
            three_prime_exon="GACAGGCATG",
            p1_length=4,
            analysis=analysis,
        )
        p1 = designer.design_p1(analysis)
        # p1_length=4, exon[-5:-1] = TCCG → RC = CGGA
        assert p1.new_igs == "CGGA"
        assert p1.p1_length == 4

    def test_wobble_valid_with_t(self):
        """Last nt of 5' exon = T → valid G·U wobble."""
        intron = "ACGACCTGAGGTAACGT" + "A" * 383 + "G"
        analysis = _make_analysis(intron)
        designer = IntronPrimerDesigner(
            intron_sequence=intron,
            five_prime_exon="GTCCGT",
            three_prime_exon="GACAGGCATG",
            analysis=analysis,
        )
        p1 = designer.design_p1(analysis)
        assert p1.wobble_valid is True
        assert p1.wobble_pair == ("G", "T")

    def test_wobble_disrupted_with_non_t(self):
        """Last nt of 5' exon ≠ T → wobble disrupted."""
        intron = "ACGACCTGAGGTAACGT" + "A" * 383 + "G"
        analysis = _make_analysis(intron)
        designer = IntronPrimerDesigner(
            intron_sequence=intron,
            five_prime_exon="GTCCGA",
            three_prime_exon="GACAGGCATG",
            analysis=analysis,
        )
        p1 = designer.design_p1(analysis)
        assert p1.wobble_valid is False


class TestDesignP1ex:
    """Tests for P1ex (P10 + stem) design — Rules 2 & 3."""

    def test_p10_reverse_complement(self):
        """P1ex-P10 = RC of 3' exon segment."""
        intron = "ACGACCTGAGGTAACGT" + "A" * 383 + "G"
        analysis = _make_analysis(intron)
        # 3' exon = GACAGGCATG, offset=4, length=3 → segment = GGC → RC = GCC
        designer = IntronPrimerDesigner(
            intron_sequence=intron,
            five_prime_exon="GTCCGT",
            three_prime_exon="GACAGGCATG",
            p10_offset=4,
            p10_length=3,
            analysis=analysis,
        )
        p1ex = designer.design_p1ex(analysis)
        assert p1ex.new_p1ex_p10 == "GCC"
        assert p1ex.exon_segment == "GGC"

    def test_stem_native_in_rational_mode(self):
        """In NONE mode, stem sequences should be native."""
        intron = "ACGACCTGAGGTAACGT" + "A" * 383 + "G"
        analysis = _make_analysis(intron)
        designer = IntronPrimerDesigner(
            intron_sequence=intron,
            five_prime_exon="GTCCGT",
            three_prime_exon="GACAGGCATG",
            library_mode=LibraryMode.NONE,
            analysis=analysis,
        )
        p1ex = designer.design_p1ex(analysis)
        assert "N" not in p1ex.stem_5prime_sequence
        assert "N" not in p1ex.stem_3prime_sequence

    def test_stem_5prime_randomised(self):
        """In FIVE_PRIME mode, 5' strand → N, 3' strand → native."""
        intron = "ACGACCTGAGGTAACGT" + "A" * 383 + "G"
        analysis = _make_analysis(intron)
        designer = IntronPrimerDesigner(
            intron_sequence=intron,
            five_prime_exon="GTCCGT",
            three_prime_exon="GACAGGCATG",
            library_mode=LibraryMode.FIVE_PRIME,
            analysis=analysis,
        )
        p1ex = designer.design_p1ex(analysis)
        assert all(c == "N" for c in p1ex.stem_5prime_sequence)
        assert "N" not in p1ex.stem_3prime_sequence

    def test_stem_3prime_randomised(self):
        """In THREE_PRIME mode, 3' strand → N, 5' strand → native."""
        intron = "ACGACCTGAGGTAACGT" + "A" * 383 + "G"
        analysis = _make_analysis(intron)
        designer = IntronPrimerDesigner(
            intron_sequence=intron,
            five_prime_exon="GTCCGT",
            three_prime_exon="GACAGGCATG",
            library_mode=LibraryMode.THREE_PRIME,
            analysis=analysis,
        )
        p1ex = designer.design_p1ex(analysis)
        assert "N" not in p1ex.stem_5prime_sequence
        assert all(c == "N" for c in p1ex.stem_3prime_sequence)

    def test_stem_both_randomised(self):
        """In BOTH mode, both strands → N."""
        intron = "ACGACCTGAGGTAACGT" + "A" * 383 + "G"
        analysis = _make_analysis(intron)
        designer = IntronPrimerDesigner(
            intron_sequence=intron,
            five_prime_exon="GTCCGT",
            three_prime_exon="GACAGGCATG",
            library_mode=LibraryMode.BOTH,
            analysis=analysis,
        )
        p1ex = designer.design_p1ex(analysis)
        assert all(c == "N" for c in p1ex.stem_5prime_sequence)
        assert all(c == "N" for c in p1ex.stem_3prime_sequence)


class TestComputeMutations:
    def test_no_mutations_when_native_matches(self):
        """If native intron already matches target, no mutations."""
        # Construct intron where positions match what design would produce
        # 5' exon GTCCGT, exon[-6:-1] = GTCCG → IGS = RC(GTCCG) = CGGAC (positions 8-12)
        # 3' exon segment at offset 4, len 3 from GACAGGCATG = GGC → P10 = GCC (positions 0-2)
        intron = "GCCACCTGCGGAC" + "A" * 387 + "G"
        analysis = _make_analysis(intron)
        designer = IntronPrimerDesigner(
            intron_sequence=intron,
            five_prime_exon="GTCCGT",
            three_prime_exon="GACAGGCATG",
            p10_offset=4,
            p10_length=3,
            library_mode=LibraryMode.NONE,
            analysis=analysis,
        )
        p1 = designer.design_p1(analysis)
        p1ex = designer.design_p1ex(analysis)
        mutations = designer.compute_mutations(analysis, p1, p1ex)
        # Only count non-stem, non-degenerate mutations
        deterministic_muts = [m for m in mutations if m.new_nt != "N"]
        # Check that no IGS or P10 mutations are needed
        igs_muts = [m for m in deterministic_muts if m.region == "IGS"]
        p10_muts = [m for m in deterministic_muts if m.region == "P1ex-P10"]
        assert len(igs_muts) == 0
        assert len(p10_muts) == 0

    def test_mutations_sorted_by_position(self):
        intron = "ACGACCTGAGGTAACGT" + "A" * 383 + "G"
        analysis = _make_analysis(intron)
        designer = IntronPrimerDesigner(
            intron_sequence=intron,
            five_prime_exon="GTCCGT",
            three_prime_exon="GACAGGCATG",
            p10_offset=4,
            p10_length=3,
            analysis=analysis,
        )
        p1 = designer.design_p1(analysis)
        p1ex = designer.design_p1ex(analysis)
        mutations = designer.compute_mutations(analysis, p1, p1ex)
        positions = [m.intron_position for m in mutations]
        assert positions == sorted(positions)

    def test_degenerate_mutations_in_library_mode(self):
        """Library mode should produce N mutations at stem positions."""
        intron = "ACGACCTGAGGTAACGT" + "A" * 383 + "G"
        analysis = _make_analysis(intron)
        designer = IntronPrimerDesigner(
            intron_sequence=intron,
            five_prime_exon="GTCCGT",
            three_prime_exon="GACAGGCATG",
            library_mode=LibraryMode.BOTH,
            analysis=analysis,
        )
        p1 = designer.design_p1(analysis)
        p1ex = designer.design_p1ex(analysis)
        mutations = designer.compute_mutations(analysis, p1, p1ex)
        n_mutations = [m for m in mutations if m.new_nt == "N"]
        assert len(n_mutations) > 0


class TestComputeLibraryInfo:
    def test_none_returns_none(self):
        intron = "ACGACCTGAGGTAACGT" + "A" * 383 + "G"
        analysis = _make_analysis(intron)
        designer = IntronPrimerDesigner(
            intron_sequence=intron,
            five_prime_exon="GTCCGT",
            three_prime_exon="GACAGGCATG",
            library_mode=LibraryMode.NONE,
            analysis=analysis,
        )
        p1ex = designer.design_p1ex(analysis)
        lib = designer.compute_library_info(analysis, p1ex)
        assert lib is None

    def test_5prime_complexity(self):
        """Complexity = 4^n₅ for 5prime mode."""
        intron = "ACGACCTGAGGTAACGT" + "A" * 383 + "G"
        analysis = _make_analysis(intron)
        designer = IntronPrimerDesigner(
            intron_sequence=intron,
            five_prime_exon="GTCCGT",
            three_prime_exon="GACAGGCATG",
            library_mode=LibraryMode.FIVE_PRIME,
            analysis=analysis,
        )
        p1ex = designer.design_p1ex(analysis)
        lib = designer.compute_library_info(analysis, p1ex)
        n5 = len(analysis.p1ex_stem.five_prime_positions)
        assert lib is not None
        assert lib.complexity == 4**n5
        assert lib.stem_5prime_randomised is True
        assert lib.stem_3prime_randomised is False

    def test_3prime_complexity(self):
        """Complexity = 4^n₃ for 3prime mode."""
        intron = "ACGACCTGAGGTAACGT" + "A" * 383 + "G"
        analysis = _make_analysis(intron)
        designer = IntronPrimerDesigner(
            intron_sequence=intron,
            five_prime_exon="GTCCGT",
            three_prime_exon="GACAGGCATG",
            library_mode=LibraryMode.THREE_PRIME,
            analysis=analysis,
        )
        p1ex = designer.design_p1ex(analysis)
        lib = designer.compute_library_info(analysis, p1ex)
        n3 = len(analysis.p1ex_stem.three_prime_positions)
        assert lib is not None
        assert lib.complexity == 4**n3

    def test_both_complexity(self):
        """Complexity = 4^(n₅+n₃) for both mode."""
        intron = "ACGACCTGAGGTAACGT" + "A" * 383 + "G"
        analysis = _make_analysis(intron)
        designer = IntronPrimerDesigner(
            intron_sequence=intron,
            five_prime_exon="GTCCGT",
            three_prime_exon="GACAGGCATG",
            library_mode=LibraryMode.BOTH,
            analysis=analysis,
        )
        p1ex = designer.design_p1ex(analysis)
        lib = designer.compute_library_info(analysis, p1ex)
        n5 = len(analysis.p1ex_stem.five_prime_positions)
        n3 = len(analysis.p1ex_stem.three_prime_positions)
        assert lib is not None
        assert lib.complexity == 4 ** (n5 + n3)
        assert lib.stem_5prime_randomised is True
        assert lib.stem_3prime_randomised is True


class TestRunPipeline:
    def test_run_rational_mode(self):
        """Full pipeline in rational mode should produce a complete report."""
        intron = "ACGACCTGAGGTAACGT" + "A" * 383 + "G"
        analysis = _make_analysis(intron)
        designer = IntronPrimerDesigner(
            intron_sequence=intron,
            five_prime_exon="GTCCGT",
            three_prime_exon="GACAGGCATG",
            p10_offset=4,
            p10_length=3,
            library_mode=LibraryMode.NONE,
            analysis=analysis,
        )
        report = designer.run()
        assert report.library is None
        assert report.forward_primer.contains_degenerate is False
        assert report.reverse_primer.contains_degenerate is False

    def test_run_library_mode(self):
        """Full pipeline in library mode should produce primers with Ns."""
        intron = "ACGACCTGAGGTAACGT" + "A" * 383 + "G"
        analysis = _make_analysis(intron)
        designer = IntronPrimerDesigner(
            intron_sequence=intron,
            five_prime_exon="GTCCGT",
            three_prime_exon="GACAGGCATG",
            p10_offset=4,
            p10_length=3,
            library_mode=LibraryMode.BOTH,
            analysis=analysis,
        )
        report = designer.run()
        assert report.library is not None
        assert report.forward_primer.contains_degenerate is True
        assert report.reverse_primer.contains_degenerate is False

    def test_run_without_analysis_raises(self):
        """Should raise if no analysis is provided."""
        designer = IntronPrimerDesigner(
            intron_sequence="A" * 400 + "G",
            five_prime_exon="GTCCGT",
            three_prime_exon="GACAGGCATG",
        )
        with pytest.raises(ValueError, match="IntronAnalysis must be provided"):
            designer.run()
