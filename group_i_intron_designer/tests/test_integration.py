"""Integration tests for the full design pipeline.

These tests use mock IntronAnalysis objects to test the end-to-end
flow without requiring Infernal to be installed.

For tests with actual Infernal, set the environment variable
GROUP_I_INTRON_TEST_INFERNAL=1 (skipped by default).
"""

from __future__ import annotations

import json
import os

import pytest

from group_i_intron_designer.designer import IntronPrimerDesigner
from group_i_intron_designer.models import (
    DesignReport,
    IntronAnalysis,
    LibraryMode,
    P1exStemInfo,
)
from group_i_intron_designer.report import format_json_report, format_text_report
from group_i_intron_designer.validators import validate_inputs


# ---------------------------------------------------------------------------
# Fixtures: biologically inspired mock introns
# ---------------------------------------------------------------------------

# Tetrahymena thermophila-like group I intron (IC1 subtype)
# The 5' end architecture:
# pos 0-2: P1ex-P10 (pairs with 3' exon)
# pos 3-6: P1ex-stem 5' strand
# pos 7-10: P1ex-stem 3' strand
# pos 11: wobble G
# pos 12-17: IGS (6 nt)
# pos 18+: intron core
TETRAHYMENA_LIKE = (
    "CUGAAUUGAUGGCUACAAU"  # 5' end (P1ex + stem + G + IGS)
    + "GAAAGUUAACUAUAAACUAA"  # intron core (simplified)
    + "ACCCUGGUUUAUAGGUAUGG"
    + "UUAAGCUUGGUAAUCCGUUA"
    + "CCGUGAAGUCAUUCCAUAAA"
    + "GGCUUAGCCUUGCGAAGAUU"
    + "GAUAGUGCCAUACCAGCAUGG"
    + "GCUGUCCCCAUAAAUAGUUU"
    + "AAAUCUUAGAUAAAAGCUGA"
    + "CCUGUCUUAUAAAUGGAUGG"
    + "GCCUCUAAACGGGUUUUGAG"
    + "GGAGCUUAAGAGUUUAAUAA"
    + "GGUAAAAAAUAACUUUACCC"
    + "UAUAGUUAAAUUUCCUAUAG"
    + "UAAAGAUCAAGGAUAGUUCC"
    + "UUGGGGUCCUUUGUUACUGG"
    + "CCUCUUUAAACCUAAUAUAA"
    + "AAAUUAAUGCUUAUAUAUAG"
    + "UCUAAAUUUAAUCUAUAGGG"
    + "CCUCUUUAAAGUAGAUCAAG"
    + "G"  # omega G
)

# Convert to DNA for the tool
TETRAHYMENA_DNA = TETRAHYMENA_LIKE.replace("U", "T")


def _make_tetrahymena_analysis() -> IntronAnalysis:
    """Mock analysis for the Tetrahymena-like intron."""
    seq = TETRAHYMENA_DNA
    return IntronAnalysis(
        is_group_i=True,
        subtype="IC1",
        score=89.3,
        e_value=2.1e-24,
        igs_positions=[12, 13, 14, 15, 16, 17],
        p1ex_p10_positions=[0, 1, 2],
        p1ex_stem=P1exStemInfo(
            five_prime_positions=[3, 4, 5, 6],
            three_prime_positions=[7, 8, 9, 10],
            five_prime_native=seq[3:7],
            three_prime_native=seq[7:11],
            stem_length=4,
        ),
        wobble_g_position=11,
        omega_g_position=len(seq) - 1,
        native_igs=seq[12:18],
        native_p1ex_p10=seq[0:3],
        native_p1_length=6,
        native_p1ex_p10_length=3,
        intron_length=len(seq),
    )


# Azoarcus-like group I intron (IC3 subtype) — shorter intron
AZOARCUS_DNA = (
    "GCCGATGCTAAGGCTA"  # 5' end
    + "GATGAGTCCGTGAGGACGAAACGGTACCCGGTACCGTC"
    + "AGCTACGATCAATCGGTAACGTAGCTGATCGATCGTAC"
    + "GCTAAATCCAGCTTGGCTAACCTAAATCGATCGGGACC"
    + "TTAAGGCCTAGCTGACTAAATCGTAACCCGATCGTTAA"
    + "CCTTGGATCGATCGAATCGGCTTAACCTGGA"
    + "G"  # omega G
)


def _make_azoarcus_analysis() -> IntronAnalysis:
    """Mock analysis for the Azoarcus-like intron."""
    seq = AZOARCUS_DNA
    return IntronAnalysis(
        is_group_i=True,
        subtype="IC3",
        score=72.1,
        e_value=3.5e-18,
        igs_positions=[10, 11, 12, 13, 14],
        p1ex_p10_positions=[0, 1, 2],
        p1ex_stem=P1exStemInfo(
            five_prime_positions=[3, 4, 5],
            three_prime_positions=[6, 7, 8],
            five_prime_native=seq[3:6],
            three_prime_native=seq[6:9],
            stem_length=3,
        ),
        wobble_g_position=9,
        omega_g_position=len(seq) - 1,
        native_igs=seq[10:15],
        native_p1ex_p10=seq[0:3],
        native_p1_length=5,
        native_p1ex_p10_length=3,
        intron_length=len(seq),
    )


# ---------------------------------------------------------------------------
# Integration tests
# ---------------------------------------------------------------------------

class TestTetrahymenaAllModes:
    """test_tetrahymena_all_modes: Full pipeline × 4 modes."""

    @pytest.mark.parametrize("mode", list(LibraryMode))
    def test_full_pipeline(self, mode):
        analysis = _make_tetrahymena_analysis()
        designer = IntronPrimerDesigner(
            intron_sequence=TETRAHYMENA_DNA,
            five_prime_exon="AGTCCGT",
            three_prime_exon="GACAGGCATGCC",
            p10_offset=4,
            p10_length=3,
            library_mode=mode,
            target_tm=60.0,
            analysis=analysis,
        )
        report = designer.run()

        # Basic invariants
        assert report.p1.new_igs is not None
        assert report.p1ex.new_p1ex_p10 is not None
        assert "N" not in report.p1.new_igs
        assert "N" not in report.p1ex.new_p1ex_p10
        assert report.forward_primer.length > 0
        assert report.reverse_primer.length > 0
        assert "N" not in report.reverse_primer.sequence

        # Library-specific checks
        if mode == LibraryMode.NONE:
            assert report.library is None
            assert "N" not in report.forward_primer.sequence
        else:
            assert report.library is not None
            assert report.library.complexity >= 1

    @pytest.mark.parametrize("mode", list(LibraryMode))
    def test_text_report_generation(self, mode):
        analysis = _make_tetrahymena_analysis()
        designer = IntronPrimerDesigner(
            intron_sequence=TETRAHYMENA_DNA,
            five_prime_exon="AGTCCGT",
            three_prime_exon="GACAGGCATGCC",
            library_mode=mode,
            analysis=analysis,
        )
        report = designer.run()
        text = format_text_report(report)
        assert "DESIGN REPORT" in text
        assert "Forward" in text
        assert "Reverse" in text

    @pytest.mark.parametrize("mode", list(LibraryMode))
    def test_json_report_generation(self, mode):
        analysis = _make_tetrahymena_analysis()
        designer = IntronPrimerDesigner(
            intron_sequence=TETRAHYMENA_DNA,
            five_prime_exon="AGTCCGT",
            three_prime_exon="GACAGGCATGCC",
            library_mode=mode,
            analysis=analysis,
        )
        report = designer.run()
        json_str = format_json_report(report)
        data = json.loads(json_str)
        assert "forward_primer" in data
        assert "reverse_primer" in data
        assert "intron_analysis" in data


class TestAzoarcusAllModes:
    """test_azoarcus_all_modes: Different intron × 4 modes."""

    @pytest.mark.parametrize("mode", list(LibraryMode))
    def test_full_pipeline(self, mode):
        analysis = _make_azoarcus_analysis()
        designer = IntronPrimerDesigner(
            intron_sequence=AZOARCUS_DNA,
            five_prime_exon="TGCCGT",
            three_prime_exon="CACAGGCATG",
            p10_offset=3,
            p10_length=3,
            library_mode=mode,
            target_tm=60.0,
            analysis=analysis,
        )
        report = designer.run()

        assert report.p1.new_igs is not None
        assert report.p1ex.new_p1ex_p10 is not None
        assert report.forward_primer.length > 0
        assert report.reverse_primer.length > 0

        if mode == LibraryMode.NONE:
            assert report.library is None
        else:
            assert report.library is not None


class TestNativeContextNoMutations:
    """test_native_context_no_mutations: Native exons → no IGS/P10 changes."""

    def test_native_exons_produce_no_igs_p10_mutations(self):
        """When exon sequences match the native pairing, IGS and P10
        should not require mutations (in rational mode)."""
        seq = TETRAHYMENA_DNA
        analysis = _make_tetrahymena_analysis()

        # Construct exon sequences that would produce the native IGS and P10
        from group_i_intron_designer.sequence_utils import reverse_complement

        native_igs = analysis.native_igs
        native_p10 = analysis.native_p1ex_p10

        # The IGS pairs with exon bases upstream of the wobble base.
        # exon = RC(native_igs) + wobble_T
        matching_5prime_exon = reverse_complement(native_igs) + "T"
        # 3' exon at offset 4, len 3 → RC = native P10 → segment = RC(native_p10)
        matching_3prime_segment = reverse_complement(native_p10)
        matching_3prime_exon = "AAAA" + matching_3prime_segment + "AAAA"

        designer = IntronPrimerDesigner(
            intron_sequence=seq,
            five_prime_exon=matching_5prime_exon,
            three_prime_exon=matching_3prime_exon,
            p10_offset=4,
            p10_length=3,
            library_mode=LibraryMode.NONE,
            analysis=analysis,
        )
        report = designer.run()

        # No IGS or P10 mutations should be needed
        igs_mutations = [m for m in report.mutations if m.region == "IGS"]
        p10_mutations = [m for m in report.mutations if m.region == "P1ex-P10"]
        assert len(igs_mutations) == 0, f"Unexpected IGS mutations: {igs_mutations}"
        assert len(p10_mutations) == 0, f"Unexpected P10 mutations: {p10_mutations}"


class TestValidationIntegration:
    """Test that validation catches bad inputs before design."""

    def test_short_exon_rejected(self):
        errors = validate_inputs(
            intron_seq="A" * 200,
            five_prime_exon="AT",  # too short for p1_length=6
            three_prime_exon="GACAGGCATG",
            library_mode=LibraryMode.NONE,
            p1_length=6,
            p10_length=3,
            p10_offset=4,
            target_tm=60.0,
        )
        assert any("5' exon" in e for e in errors)

    def test_short_3prime_exon_rejected(self):
        errors = validate_inputs(
            intron_seq="A" * 200,
            five_prime_exon="AGTCCGT",
            three_prime_exon="GAC",  # too short for offset=4 + length=3
            library_mode=LibraryMode.NONE,
            p1_length=5,
            p10_length=3,
            p10_offset=4,
            target_tm=60.0,
        )
        assert any("3' exon" in e for e in errors)

    def test_valid_inputs_pass(self):
        errors = validate_inputs(
            intron_seq="A" * 200,
            five_prime_exon="AGTCCGT",
            three_prime_exon="GACAGGCATG",
            library_mode=LibraryMode.NONE,
            p1_length=5,
            p10_length=3,
            p10_offset=4,
            target_tm=60.0,
        )
        assert errors == []


class TestWarningsIntegration:
    """Test that appropriate warnings are generated."""

    def test_wobble_disruption_warning(self):
        analysis = _make_tetrahymena_analysis()
        designer = IntronPrimerDesigner(
            intron_sequence=TETRAHYMENA_DNA,
            five_prime_exon="AGTCCGA",  # ends in A, not T
            three_prime_exon="GACAGGCATGCC",
            library_mode=LibraryMode.NONE,
            analysis=analysis,
        )
        report = designer.run()
        assert any("wobble" in w.lower() for w in report.warnings)

    def test_long_p10_warning(self):
        analysis = _make_tetrahymena_analysis()
        designer = IntronPrimerDesigner(
            intron_sequence=TETRAHYMENA_DNA,
            five_prime_exon="AGTCCGT",
            three_prime_exon="GACAGGCATGCC",
            p10_length=5,  # > 3 bp
            p10_offset=2,
            library_mode=LibraryMode.NONE,
            analysis=analysis,
        )
        report = designer.run()
        assert any("P10" in w for w in report.warnings)


class TestReportFormats:
    """Test both output formats produce valid output."""

    def test_text_report_structure(self):
        analysis = _make_tetrahymena_analysis()
        designer = IntronPrimerDesigner(
            intron_sequence=TETRAHYMENA_DNA,
            five_prime_exon="AGTCCGT",
            three_prime_exon="GACAGGCATGCC",
            analysis=analysis,
        )
        report = designer.run()
        text = format_text_report(report)

        # Check for expected sections
        assert "Intron Analysis" in text
        assert "New P1 Helix" in text
        assert "P1ex Design" in text
        assert "Mutations" in text
        assert "Primers" in text

    def test_json_report_roundtrip(self):
        analysis = _make_tetrahymena_analysis()
        designer = IntronPrimerDesigner(
            intron_sequence=TETRAHYMENA_DNA,
            five_prime_exon="AGTCCGT",
            three_prime_exon="GACAGGCATGCC",
            library_mode=LibraryMode.BOTH,
            analysis=analysis,
        )
        report = designer.run()
        json_str = format_json_report(report)
        data = json.loads(json_str)

        # Verify key fields exist and have correct types
        assert isinstance(data["forward_primer"]["sequence"], str)
        assert isinstance(data["forward_primer"]["length"], int)
        assert isinstance(data["intron_analysis"]["score"], float)
        assert data["p1ex"]["library_mode"] == "both"
        assert data["library"]["complexity"] > 1
