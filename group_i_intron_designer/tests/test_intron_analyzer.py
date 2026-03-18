"""Tests for the intron_analyzer module.

These tests primarily test the structural parser since the Infernal
integration requires external tools.
"""

import pytest

from group_i_intron_designer.structural_parser import (
    _parse_wuss_pairs,
    _parse_pseudoknot_pairs,
    _ungapped_index_map,
    parse_cmscan_tblout,
)


class TestUngappedIndexMap:
    def test_no_gaps(self):
        mapping = _ungapped_index_map("ACGT")
        assert mapping == {0: 0, 1: 1, 2: 2, 3: 3}

    def test_with_gaps(self):
        mapping = _ungapped_index_map("AC.GT")
        assert mapping == {0: 0, 1: 1, 3: 2, 4: 3}

    def test_leading_gaps(self):
        mapping = _ungapped_index_map("..ACG")
        assert mapping == {2: 0, 3: 1, 4: 2}

    def test_dash_gaps(self):
        mapping = _ungapped_index_map("A-C-G")
        assert mapping == {0: 0, 2: 1, 4: 2}

    def test_empty(self):
        assert _ungapped_index_map("") == {}


class TestParseWussPairs:
    def test_simple_stem(self):
        pairs = _parse_wuss_pairs("((..))")
        assert pairs[0] == 5
        assert pairs[1] == 4
        assert pairs[4] == 1
        assert pairs[5] == 0

    def test_nested(self):
        pairs = _parse_wuss_pairs("(((...)))")
        assert pairs[0] == 8
        assert pairs[1] == 7
        assert pairs[2] == 6

    def test_angle_brackets(self):
        pairs = _parse_wuss_pairs("<<..>>")
        assert pairs[0] == 5
        assert pairs[1] == 4

    def test_mixed_brackets(self):
        pairs = _parse_wuss_pairs("(<..>)")
        assert pairs[0] == 5  # ( pairs with )
        assert pairs[1] == 4  # < pairs with >

    def test_unpaired_not_in_pairs(self):
        pairs = _parse_wuss_pairs("(....)")
        assert 1 not in pairs
        assert 2 not in pairs
        assert 3 not in pairs
        assert 4 not in pairs

    def test_empty(self):
        pairs = _parse_wuss_pairs("")
        assert pairs == {}


class TestParsePseudoknotPairs:
    def test_simple_pk(self):
        pairs = _parse_pseudoknot_pairs("AA....AA")
        assert pairs[0] == 7
        assert pairs[1] == 6

    def test_different_letters(self):
        pairs = _parse_pseudoknot_pairs("AAA...BBB...AAA...BBB")
        # Each letter group's first half pairs with second half
        assert 0 in pairs  # A group

    def test_lowercase(self):
        pairs = _parse_pseudoknot_pairs("aa....aa")
        assert pairs[0] == 7
        assert pairs[1] == 6

    def test_no_letters(self):
        pairs = _parse_pseudoknot_pairs("((..))")
        assert pairs == {}


class TestParseCmscanTblout:
    def test_empty(self):
        assert parse_cmscan_tblout("") == []

    def test_comment_lines_skipped(self):
        text = "# comment line\n# another\n"
        assert parse_cmscan_tblout(text) == []

    def test_sorted_by_evalue(self):
        # Create mock tblout lines (18+ fields)
        line1 = "IC1 - query - cm 1 100 1 100 + no 1 - 0.0 1.0 50.0 1e-10 - -"
        line2 = "IC3 - query - cm 1 100 1 100 + no 1 - 0.0 2.0 40.0 1e-5 - -"
        text = f"{line1}\n{line2}\n"
        hits = parse_cmscan_tblout(text)
        assert len(hits) == 2
        assert hits[0]["e_value"] <= hits[1]["e_value"]
