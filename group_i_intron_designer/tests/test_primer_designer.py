"""Tests for the primer_designer module."""

import pytest

from group_i_intron_designer.models import LibraryMode, P1exStemInfo
from group_i_intron_designer.primer_designer import (
    build_forward_primer,
    build_reverse_primer,
)


@pytest.fixture
def intron_seq() -> str:
    """A mock intron sequence."""
    return "ACGACCTGAGGTAACG" + "T" * 384 + "G"


@pytest.fixture
def stem_info() -> P1exStemInfo:
    return P1exStemInfo(
        five_prime_positions=[3, 4],
        three_prime_positions=[5, 6],
        five_prime_native="AC",
        three_prime_native="CT",
        stem_length=2,
    )


class TestBuildForwardPrimer:
    def test_includes_mutation_and_annealing_regions(
        self, intron_seq, stem_info
    ):
        primer = build_forward_primer(
            intron_seq=intron_seq,
            new_p1ex_p10="GCC",
            new_igs="ACGGAC",
            p1ex_p10_positions=[0, 1, 2],
            stem=stem_info,
            igs_positions=[8, 9, 10, 11, 12],
            wobble_g_position=7,
            library_mode=LibraryMode.NONE,
            target_tm=60.0,
        )
        assert len(primer.mutation_region) > 0
        assert len(primer.annealing_region) > 0
        assert primer.sequence == primer.mutation_region + primer.annealing_region

    def test_no_n_in_rational_mode(self, intron_seq, stem_info):
        primer = build_forward_primer(
            intron_seq=intron_seq,
            new_p1ex_p10="GCC",
            new_igs="ACGGAC",
            p1ex_p10_positions=[0, 1, 2],
            stem=stem_info,
            igs_positions=[8, 9, 10, 11, 12],
            wobble_g_position=7,
            library_mode=LibraryMode.NONE,
            target_tm=60.0,
        )
        assert "N" not in primer.sequence
        assert primer.contains_degenerate is False

    def test_n_in_5prime_mode(self, intron_seq, stem_info):
        primer = build_forward_primer(
            intron_seq=intron_seq,
            new_p1ex_p10="GCC",
            new_igs="ACGGAC",
            p1ex_p10_positions=[0, 1, 2],
            stem=stem_info,
            igs_positions=[8, 9, 10, 11, 12],
            wobble_g_position=7,
            library_mode=LibraryMode.FIVE_PRIME,
            target_tm=60.0,
        )
        assert primer.contains_degenerate is True
        assert primer.degenerate_count == 2

    def test_n_in_3prime_mode(self, intron_seq, stem_info):
        primer = build_forward_primer(
            intron_seq=intron_seq,
            new_p1ex_p10="GCC",
            new_igs="ACGGAC",
            p1ex_p10_positions=[0, 1, 2],
            stem=stem_info,
            igs_positions=[8, 9, 10, 11, 12],
            wobble_g_position=7,
            library_mode=LibraryMode.THREE_PRIME,
            target_tm=60.0,
        )
        assert primer.contains_degenerate is True
        assert primer.degenerate_count == 2

    def test_n_in_both_mode(self, intron_seq, stem_info):
        primer = build_forward_primer(
            intron_seq=intron_seq,
            new_p1ex_p10="GCC",
            new_igs="ACGGAC",
            p1ex_p10_positions=[0, 1, 2],
            stem=stem_info,
            igs_positions=[8, 9, 10, 11, 12],
            wobble_g_position=7,
            library_mode=LibraryMode.BOTH,
            target_tm=60.0,
        )
        assert primer.contains_degenerate is True
        assert primer.degenerate_count == 4

    def test_wobble_g_always_present(self, intron_seq, stem_info):
        """Wobble G at position 7 should always be G in the primer."""
        for mode in LibraryMode:
            primer = build_forward_primer(
                intron_seq=intron_seq,
                new_p1ex_p10="GCC",
                new_igs="ACGGAC",
                p1ex_p10_positions=[0, 1, 2],
                stem=stem_info,
                igs_positions=[8, 9, 10, 11, 12],
                wobble_g_position=7,
                library_mode=mode,
                target_tm=60.0,
            )
            # Position 7 in the intron = index 7 in the mutation region
            # (since mutation region starts at position 0)
            mut = primer.mutation_region
            assert mut[7] == "G", f"Wobble G missing in mode {mode}"

    def test_igs_correctly_placed(self, intron_seq, stem_info):
        """IGS should be correctly placed in the mutation region."""
        primer = build_forward_primer(
            intron_seq=intron_seq,
            new_p1ex_p10="GCC",
            new_igs="ACGGAC",
            p1ex_p10_positions=[0, 1, 2],
            stem=stem_info,
            igs_positions=[8, 9, 10, 11, 12],
            wobble_g_position=7,
            library_mode=LibraryMode.NONE,
            target_tm=60.0,
        )
        mut = primer.mutation_region
        # IGS positions 8-12, first position in mutation region is 0
        # So IGS starts at index 8
        igs_in_primer = mut[8:13]
        assert igs_in_primer == "ACGGA"  # First 5 chars of IGS at those positions

    def test_p10_correctly_placed(self, intron_seq, stem_info):
        """P1ex-P10 should be at the start of the mutation region."""
        primer = build_forward_primer(
            intron_seq=intron_seq,
            new_p1ex_p10="GCC",
            new_igs="ACGGAC",
            p1ex_p10_positions=[0, 1, 2],
            stem=stem_info,
            igs_positions=[8, 9, 10, 11, 12],
            wobble_g_position=7,
            library_mode=LibraryMode.NONE,
            target_tm=60.0,
        )
        mut = primer.mutation_region
        assert mut[:3] == "GCC"


class TestBuildReversePrimer:
    def test_no_degenerate_bases(self, intron_seq):
        primer = build_reverse_primer(
            intron_seq=intron_seq,
            five_prime_exon="GAAGTCAAGTCCGTAGCA",
            first_modified_pos=0,
            target_tm=60.0,
        )
        assert "N" not in primer.sequence
        assert primer.contains_degenerate is False
        assert primer.degenerate_count == 0

    def test_name_is_reverse(self, intron_seq):
        primer = build_reverse_primer(
            intron_seq=intron_seq,
            five_prime_exon="GAAGTCAAGTCCGTAGCA",
            first_modified_pos=0,
            target_tm=60.0,
        )
        assert primer.name == "Reverse"

    def test_minimum_length(self, intron_seq):
        primer = build_reverse_primer(
            intron_seq=intron_seq,
            five_prime_exon="GAAGTCAAGTCCGTAGCATTGCAG",
            first_modified_pos=0,
            target_tm=60.0,
        )
        assert primer.length >= 18

    def test_anneals_to_exon(self, intron_seq):
        """Reverse primer should anneal to the 5' exon when mutations
        start at position 0."""
        exon = "GAAGTCAAGTCCGTAGCATTGCAG"
        primer = build_reverse_primer(
            intron_seq=intron_seq,
            five_prime_exon=exon,
            first_modified_pos=0,
            target_tm=60.0,
        )
        from group_i_intron_designer.sequence_utils import reverse_complement
        sense = reverse_complement(primer.sequence)
        assert exon.endswith(sense)
