"""Tests for the helix_drawer module."""

from group_i_intron_designer.helix_drawer import draw_p1_helix, draw_p10_helix
from group_i_intron_designer.models import LibraryMode


class TestDrawP1Helix:
    def test_contains_exon_and_igs(self):
        result = draw_p1_helix("GTCCGT", "ACGGAC", wobble_valid=True)
        assert "G T C C G T" in result
        assert "C A G G C A" in result

    def test_wobble_valid_label(self):
        result = draw_p1_helix("GTCCGT", "ACGGAC", wobble_valid=True)
        assert "G·U wobble" in result

    def test_wobble_disrupted_label(self):
        result = draw_p1_helix("GTCCGA", "TCGGAC", wobble_valid=False)
        assert "disrupted" in result

    def test_has_pairing_symbols(self):
        result = draw_p1_helix("GTCCGT", "ACGGAC", wobble_valid=True)
        assert "|" in result  # Watson-Crick pairs


class TestDrawP10Helix:
    def test_contains_exon_and_p10(self):
        result = draw_p10_helix("GGC", "GCC", "AC", "CT", LibraryMode.NONE)
        assert "G G C" in result
        assert "C C G" in result

    def test_stem_native_label(self):
        result = draw_p10_helix("GGC", "GCC", "AC", "CT", LibraryMode.NONE)
        assert "native" in result

    def test_stem_randomised_label(self):
        result = draw_p10_helix("GGC", "GCC", "NN", "CT", LibraryMode.FIVE_PRIME)
        assert "randomised" in result

    def test_stem_both_randomised(self):
        result = draw_p10_helix("GGC", "GCC", "NN", "NN", LibraryMode.BOTH)
        assert result.count("randomised") >= 2
