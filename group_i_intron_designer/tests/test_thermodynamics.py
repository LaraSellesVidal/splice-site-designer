"""Tests for the thermodynamic assessment module and its integration.

Exercises the ViennaRNA-based P1/P10 scoring when ViennaRNA is installed, checks
the None-score behaviour when it is not, and confirms that enabling scoring does
not alter the design output.
"""

from __future__ import annotations

import pytest

from group_i_intron_designer import thermodynamics as thm
from group_i_intron_designer.designer import IntronPrimerDesigner
from group_i_intron_designer.models import (
    IntronAnalysis,
    LibraryMode,
    P1exStemInfo,
    ThermodynamicAssessment,
)
from group_i_intron_designer.report import format_json_report, format_text_report

HAVE_VIENNA = thm.viennarna_available()
vienna_only = pytest.mark.skipif(not HAVE_VIENNA, reason="ViennaRNA not installed")


def _make_analysis() -> IntronAnalysis:
    intron = "ACGACCTGAGGTAACGT" + "A" * 383 + "G"
    igs_positions = [8, 9, 10, 11, 12]
    p10_positions = [0, 1, 2]
    stem_5p = [3, 4]
    stem_3p = [5, 6]
    return IntronAnalysis(
        is_group_i=True,
        subtype="IC1",
        score=89.3,
        e_value=2.1e-24,
        igs_positions=igs_positions,
        p1ex_p10_positions=p10_positions,
        p1ex_stem=P1exStemInfo(stem_5p, stem_3p,
                               "".join(intron[p] for p in stem_5p),
                               "".join(intron[p] for p in stem_3p), 2),
        wobble_g_position=7,
        omega_g_position=len(intron) - 1,
        native_igs="".join(intron[p] for p in igs_positions),
        native_p1ex_p10="".join(intron[p] for p in p10_positions),
        native_p1_length=5,
        native_p1ex_p10_length=3,
        intron_length=len(intron),
    )


def _designer(mode=LibraryMode.NONE) -> IntronPrimerDesigner:
    intron = "ACGACCTGAGGTAACGT" + "A" * 383 + "G"
    return IntronPrimerDesigner(
        intron_sequence=intron,
        five_prime_exon="GGGACGCGACTGAATGAAATGGTGAAGGACGGGTCCAGTAGTT",  # ends in T
        three_prime_exon="GAGCTCCGTAACTAGTCGCGTC",
        library_mode=mode,
        analysis=_make_analysis(),
    )


class TestModuleContract:
    def test_available_flag_is_bool(self):
        assert isinstance(thm.viennarna_available(), bool)

    def test_classify_stability_thresholds(self):
        assert thm.classify_stability(None) == "unknown"
        assert thm.classify_stability(-10.0) == "strong"
        assert thm.classify_stability(-5.0) == "moderate"
        assert thm.classify_stability(-1.0) == "weak"

    def test_assess_returns_all_keys(self):
        a = thm.assess("GCGCGC", "GCGCGC", "AUAU", "AUAU")
        for k in ("viennarna", "p1_dg", "p1_stability", "p1_specificity_margin",
                  "p10_dg", "p10_stability", "warnings"):
            assert k in a


class TestGracefulDegradation:
    def test_duplex_dg_short_returns_none(self):
        assert thm.duplex_dg("A", "U") is None

    def test_specificity_none_without_vienna(self, monkeypatch):
        monkeypatch.setattr(thm, "_HAVE_VIENNA", False)
        assert thm.specificity_margin("GCGCGC", "GCGCGC") == (None, None)

    def test_assess_without_vienna_warns(self, monkeypatch):
        monkeypatch.setattr(thm, "_HAVE_VIENNA", False)
        a = thm.assess("GCGCGC", "GCGCGC", "AUAU", "AUAU")
        assert a["p1_dg"] is None
        assert any("ViennaRNA" in w for w in a["warnings"])


@vienna_only
class TestScoring:
    def test_gc_rich_is_strong(self):
        dg = thm.duplex_dg("GCGCGC", "GCGCGC")
        assert dg is not None and dg < -8.0

    def test_complex_exon_is_specific(self):
        # A complex exon should pair its exact RC better than composition shuffles.
        margin, frac = thm.specificity_margin("GAUCGA", "UCGAUC")
        assert margin is not None and margin > 1.0
        assert frac is not None and frac < 0.3

    def test_homopolymer_has_no_specificity(self):
        margin, frac = thm.specificity_margin("AAAAAA", "UUUUUU")
        # Every shuffle of a homopolymer is identical → zero margin, all tie.
        assert margin == 0.0
        assert frac == 1.0


class TestIntegrationPreservesDesign:
    """Enabling thermodynamics/validation must not change the rational design."""

    def test_design_identical_with_and_without_thermo(self):
        r_off = _designer().run(assess_thermodynamics=False)
        r_on = _designer().run(assess_thermodynamics=True)
        assert r_on.p1.new_igs == r_off.p1.new_igs
        assert r_on.p1ex.new_p1ex_p10 == r_off.p1ex.new_p1ex_p10
        assert r_on.forward_primer.sequence == r_off.forward_primer.sequence
        assert r_on.reverse_primer.sequence == r_off.reverse_primer.sequence
        assert [m.new_nt for m in r_on.mutations] == [m.new_nt for m in r_off.mutations]

    def test_thermo_absent_when_disabled(self):
        r = _designer().run(assess_thermodynamics=False)
        assert r.thermodynamics is None

    def test_thermo_present_when_enabled(self):
        r = _designer().run(assess_thermodynamics=True)
        assert isinstance(r.thermodynamics, ThermodynamicAssessment)
        assert r.thermodynamics.viennarna_available == HAVE_VIENNA

    @vienna_only
    def test_thermo_scores_populated(self):
        r = _designer().run(assess_thermodynamics=True)
        assert r.thermodynamics.p1_dg is not None
        assert r.thermodynamics.p1_stability in ("strong", "moderate", "weak")

    def test_reports_render_with_thermo(self):
        r = _designer().run(assess_thermodynamics=True)
        text = format_text_report(r)
        assert "Thermodynamic Assessment" in text
        js = format_json_report(r)
        assert "thermodynamics" in js

    def test_construct_validation_default_off(self):
        r = _designer().run()
        assert r.construct_validation is None
