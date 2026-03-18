"""Tests for sequence_utils module."""

import math

import pytest

from group_i_intron_designer.sequence_utils import (
    complement,
    gc_content,
    parse_fasta,
    reverse_complement,
    tm_nearest_neighbour,
)


class TestComplement:
    def test_basic(self):
        assert complement("ATGC") == "TACG"

    def test_single_base(self):
        assert complement("A") == "T"
        assert complement("T") == "A"
        assert complement("G") == "C"
        assert complement("C") == "G"

    def test_lowercase_input(self):
        assert complement("atgc") == "TACG"

    def test_degenerate_n(self):
        assert complement("ANG") == "TNC"

    def test_invalid_base(self):
        with pytest.raises(ValueError, match="Non-DNA character"):
            complement("ATXGC")

    def test_empty(self):
        assert complement("") == ""


class TestReverseComplement:
    def test_basic(self):
        assert reverse_complement("ATGC") == "GCAT"

    def test_palindrome(self):
        assert reverse_complement("ATAT") == "ATAT"

    def test_single_base(self):
        assert reverse_complement("A") == "T"

    def test_self_inverse(self):
        """RC of RC should be the original sequence."""
        seq = "GCTAGCTA"
        assert reverse_complement(reverse_complement(seq)) == seq

    def test_biological_example_p1(self):
        """5' exon GTCCGT → RC = ACGGAC (expected IGS)."""
        assert reverse_complement("GTCCGT") == "ACGGAC"

    def test_biological_example_p10(self):
        """3' exon segment GCA → RC = TGC (expected P1ex-P10)."""
        assert reverse_complement("GCA") == "TGC"


class TestGcContent:
    def test_all_gc(self):
        assert gc_content("GGCC") == 1.0

    def test_all_at(self):
        assert gc_content("AATT") == 0.0

    def test_half(self):
        assert gc_content("ATGC") == 0.5

    def test_empty(self):
        assert gc_content("") == 0.0

    def test_with_n(self):
        """N bases should be excluded from both numerator and denominator."""
        assert gc_content("GN") == 1.0

    def test_case_insensitive(self):
        assert gc_content("atgc") == 0.5


class TestTmNearestNeighbour:
    def test_short_oligo(self):
        """A typical 18-mer should have a reasonable Tm."""
        tm = tm_nearest_neighbour("GCTAGCTAGCTAGCTAGC")
        assert 40 < tm < 80

    def test_gc_rich_higher_tm(self):
        """GC-rich oligos should have higher Tm than AT-rich."""
        tm_gc = tm_nearest_neighbour("GCGCGCGCGCGCGCGCGC")
        tm_at = tm_nearest_neighbour("ATATATATATATATATATA")
        assert tm_gc > tm_at

    def test_too_short(self):
        with pytest.raises(ValueError, match="at least 2 nt"):
            tm_nearest_neighbour("A")

    def test_returns_float(self):
        tm = tm_nearest_neighbour("ATGCATGCATGCATGCATGC")
        assert isinstance(tm, float)

    def test_concentration_effect(self):
        """Higher oligo concentration should increase Tm."""
        tm_low = tm_nearest_neighbour("GCTAGCTAGCTAGCTAGC", oligo_conc_nm=100)
        tm_high = tm_nearest_neighbour("GCTAGCTAGCTAGCTAGC", oligo_conc_nm=1000)
        assert tm_high > tm_low


class TestParseFasta:
    def test_simple(self):
        header, seq = parse_fasta(">test\nATGCATGC\n")
        assert header == "test"
        assert seq == "ATGCATGC"

    def test_multiline_sequence(self):
        header, seq = parse_fasta(">test\nATGC\nATGC\n")
        assert seq == "ATGCATGC"

    def test_no_header(self):
        with pytest.raises(ValueError, match="no '>' header"):
            parse_fasta("ATGCATGC")

    def test_empty_sequence(self):
        with pytest.raises(ValueError, match="no sequence data"):
            parse_fasta(">test\n")

    def test_lowercase_uppercased(self):
        _, seq = parse_fasta(">test\natgcatgc\n")
        assert seq == "ATGCATGC"

    def test_multiple_records_takes_first(self):
        header, seq = parse_fasta(">first\nATGC\n>second\nGGGG\n")
        assert header == "first"
        assert seq == "ATGC"

    def test_spaces_removed(self):
        _, seq = parse_fasta(">test\nATGC ATGC\n")
        assert seq == "ATGCATGC"
