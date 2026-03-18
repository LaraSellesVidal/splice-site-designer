"""Tests for library mode behaviour across all four modes.

These tests verify the key invariants specified in the plan:
- P1ex-P10 is NEVER randomised (always rationally designed)
- IGS is NEVER randomised (always rationally designed)
- Wobble G is always preserved
- N positions in the primer match the selected stem strands
- Reverse primer never has degenerate positions
"""

import pytest

from group_i_intron_designer.designer import IntronPrimerDesigner
from group_i_intron_designer.models import (
    IntronAnalysis,
    LibraryMode,
    P1exStemInfo,
)


@pytest.fixture
def mock_intron() -> str:
    """A mock intron sequence with known positions."""
    # Positions:
    # 0-2: P1ex-P10 = ACG
    # 3-4: P1ex-stem 5' = AC
    # 5-6: P1ex-stem 3' = CT
    # 7: Wobble G
    # 8-12: IGS = GAGGT
    # 13-16: annealing region
    # ... rest of intron ... last position = G (omega)
    return "ACGACCTGAGGTAACG" + "T" * 384 + "G"


@pytest.fixture
def mock_analysis(mock_intron: str) -> IntronAnalysis:
    return IntronAnalysis(
        is_group_i=True,
        subtype="IC1",
        score=89.3,
        e_value=2.1e-24,
        igs_positions=[8, 9, 10, 11, 12],
        p1ex_p10_positions=[0, 1, 2],
        p1ex_stem=P1exStemInfo(
            five_prime_positions=[3, 4],
            three_prime_positions=[5, 6],
            five_prime_native="AC",
            three_prime_native="CT",
            stem_length=2,
        ),
        wobble_g_position=7,
        omega_g_position=len(mock_intron) - 1,
        native_igs="GAGGT",
        native_p1ex_p10="ACG",
        native_p1_length=5,
        native_p1ex_p10_length=3,
        intron_length=len(mock_intron),
    )


class TestLibraryNoneNoN:
    """test_library_none_no_N: No degenerate positions in rational mode."""

    def test_no_n_in_forward_primer(self, mock_intron, mock_analysis):
        designer = IntronPrimerDesigner(
            intron_sequence=mock_intron,
            five_prime_exon="GTCCGT",
            three_prime_exon="GACAGGCATG",
            library_mode=LibraryMode.NONE,
            analysis=mock_analysis,
        )
        report = designer.run()
        assert "N" not in report.forward_primer.sequence
        assert report.forward_primer.degenerate_count == 0

    def test_no_n_in_reverse_primer(self, mock_intron, mock_analysis):
        designer = IntronPrimerDesigner(
            intron_sequence=mock_intron,
            five_prime_exon="GTCCGT",
            three_prime_exon="GACAGGCATG",
            library_mode=LibraryMode.NONE,
            analysis=mock_analysis,
        )
        report = designer.run()
        assert "N" not in report.reverse_primer.sequence


class TestLibrary5PrimeOnly:
    """test_library_5prime_only: Only 5' strand positions become N."""

    def test_5prime_strand_randomised(self, mock_intron, mock_analysis):
        designer = IntronPrimerDesigner(
            intron_sequence=mock_intron,
            five_prime_exon="GTCCGT",
            three_prime_exon="GACAGGCATG",
            library_mode=LibraryMode.FIVE_PRIME,
            analysis=mock_analysis,
        )
        report = designer.run()
        assert report.p1ex.stem_5prime_sequence == "NN"
        assert "N" not in report.p1ex.stem_3prime_sequence

    def test_forward_primer_has_correct_n_count(self, mock_intron, mock_analysis):
        designer = IntronPrimerDesigner(
            intron_sequence=mock_intron,
            five_prime_exon="GTCCGT",
            three_prime_exon="GACAGGCATG",
            library_mode=LibraryMode.FIVE_PRIME,
            analysis=mock_analysis,
        )
        report = designer.run()
        assert report.forward_primer.degenerate_count == 2


class TestLibrary3PrimeOnly:
    """test_library_3prime_only: Only 3' strand positions become N."""

    def test_3prime_strand_randomised(self, mock_intron, mock_analysis):
        designer = IntronPrimerDesigner(
            intron_sequence=mock_intron,
            five_prime_exon="GTCCGT",
            three_prime_exon="GACAGGCATG",
            library_mode=LibraryMode.THREE_PRIME,
            analysis=mock_analysis,
        )
        report = designer.run()
        assert "N" not in report.p1ex.stem_5prime_sequence
        assert report.p1ex.stem_3prime_sequence == "NN"

    def test_forward_primer_has_correct_n_count(self, mock_intron, mock_analysis):
        designer = IntronPrimerDesigner(
            intron_sequence=mock_intron,
            five_prime_exon="GTCCGT",
            three_prime_exon="GACAGGCATG",
            library_mode=LibraryMode.THREE_PRIME,
            analysis=mock_analysis,
        )
        report = designer.run()
        assert report.forward_primer.degenerate_count == 2


class TestLibraryBoth:
    """test_library_both: Both strands become N."""

    def test_both_strands_randomised(self, mock_intron, mock_analysis):
        designer = IntronPrimerDesigner(
            intron_sequence=mock_intron,
            five_prime_exon="GTCCGT",
            three_prime_exon="GACAGGCATG",
            library_mode=LibraryMode.BOTH,
            analysis=mock_analysis,
        )
        report = designer.run()
        assert report.p1ex.stem_5prime_sequence == "NN"
        assert report.p1ex.stem_3prime_sequence == "NN"

    def test_forward_primer_has_correct_n_count(self, mock_intron, mock_analysis):
        designer = IntronPrimerDesigner(
            intron_sequence=mock_intron,
            five_prime_exon="GTCCGT",
            three_prime_exon="GACAGGCATG",
            library_mode=LibraryMode.BOTH,
            analysis=mock_analysis,
        )
        report = designer.run()
        assert report.forward_primer.degenerate_count == 4


class TestP1exP10NeverRandomised:
    """test_p1ex_p10_never_randomised: P1ex-P10 always designed, all modes."""

    @pytest.mark.parametrize(
        "mode",
        [LibraryMode.NONE, LibraryMode.FIVE_PRIME, LibraryMode.THREE_PRIME, LibraryMode.BOTH],
    )
    def test_p10_is_designed_not_random(self, mock_intron, mock_analysis, mode):
        designer = IntronPrimerDesigner(
            intron_sequence=mock_intron,
            five_prime_exon="GTCCGT",
            three_prime_exon="GACAGGCATG",
            library_mode=mode,
            analysis=mock_analysis,
        )
        report = designer.run()
        # P1ex-P10 should be the reverse complement of 3' exon segment
        assert "N" not in report.p1ex.new_p1ex_p10


class TestIgsNeverRandomised:
    """test_igs_never_randomised: IGS always designed, all modes."""

    @pytest.mark.parametrize(
        "mode",
        [LibraryMode.NONE, LibraryMode.FIVE_PRIME, LibraryMode.THREE_PRIME, LibraryMode.BOTH],
    )
    def test_igs_is_designed_not_random(self, mock_intron, mock_analysis, mode):
        designer = IntronPrimerDesigner(
            intron_sequence=mock_intron,
            five_prime_exon="GTCCGT",
            three_prime_exon="GACAGGCATG",
            library_mode=mode,
            analysis=mock_analysis,
        )
        report = designer.run()
        assert "N" not in report.p1.new_igs


class TestWobbleGPreserved:
    """test_wobble_g_preserved: Wobble G fixed in all modes."""

    @pytest.mark.parametrize(
        "mode",
        [LibraryMode.NONE, LibraryMode.FIVE_PRIME, LibraryMode.THREE_PRIME, LibraryMode.BOTH],
    )
    def test_wobble_g_in_primer(self, mock_intron, mock_analysis, mode):
        designer = IntronPrimerDesigner(
            intron_sequence=mock_intron,
            five_prime_exon="GTCCGT",
            three_prime_exon="GACAGGCATG",
            library_mode=mode,
            analysis=mock_analysis,
        )
        report = designer.run()
        # The wobble G position in the mutation region should be G
        fwd_seq = report.forward_primer.sequence
        # The wobble G is between stem and IGS in the primer
        # It should always be G, never N
        mutation_region = report.forward_primer.mutation_region
        # Position 7 (wobble G) is the 8th position from start, index 7 in mutation region
        wobble_idx = mock_analysis.wobble_g_position - min(mock_analysis.p1ex_p10_positions)
        if wobble_idx < len(mutation_region):
            assert mutation_region[wobble_idx] == "G"


class TestReversePrimerNoN:
    """test_reverse_primer_no_N: Reverse primer never has degenerate positions."""

    @pytest.mark.parametrize(
        "mode",
        [LibraryMode.NONE, LibraryMode.FIVE_PRIME, LibraryMode.THREE_PRIME, LibraryMode.BOTH],
    )
    def test_reverse_no_degenerate(self, mock_intron, mock_analysis, mode):
        designer = IntronPrimerDesigner(
            intron_sequence=mock_intron,
            five_prime_exon="GTCCGT",
            three_prime_exon="GACAGGCATG",
            library_mode=mode,
            analysis=mock_analysis,
        )
        report = designer.run()
        assert "N" not in report.reverse_primer.sequence
        assert report.reverse_primer.contains_degenerate is False
        assert report.reverse_primer.degenerate_count == 0


class TestPrimerNPositionsMatchMode:
    """test_primer_N_positions_match_mode: Ns in primer match selected strands."""

    def test_5prime_n_positions(self, mock_intron, mock_analysis):
        designer = IntronPrimerDesigner(
            intron_sequence=mock_intron,
            five_prime_exon="GTCCGT",
            three_prime_exon="GACAGGCATG",
            library_mode=LibraryMode.FIVE_PRIME,
            analysis=mock_analysis,
        )
        report = designer.run()
        # N count should equal number of 5' strand positions
        assert report.forward_primer.degenerate_count == len(
            mock_analysis.p1ex_stem.five_prime_positions
        )

    def test_3prime_n_positions(self, mock_intron, mock_analysis):
        designer = IntronPrimerDesigner(
            intron_sequence=mock_intron,
            five_prime_exon="GTCCGT",
            three_prime_exon="GACAGGCATG",
            library_mode=LibraryMode.THREE_PRIME,
            analysis=mock_analysis,
        )
        report = designer.run()
        assert report.forward_primer.degenerate_count == len(
            mock_analysis.p1ex_stem.three_prime_positions
        )

    def test_both_n_positions(self, mock_intron, mock_analysis):
        designer = IntronPrimerDesigner(
            intron_sequence=mock_intron,
            five_prime_exon="GTCCGT",
            three_prime_exon="GACAGGCATG",
            library_mode=LibraryMode.BOTH,
            analysis=mock_analysis,
        )
        report = designer.run()
        expected = len(mock_analysis.p1ex_stem.five_prime_positions) + len(
            mock_analysis.p1ex_stem.three_prime_positions
        )
        assert report.forward_primer.degenerate_count == expected
