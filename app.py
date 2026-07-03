"""Streamlit web app for the Group I Intron P1/P1ex Primer Designer."""

from __future__ import annotations

import re

import streamlit as st

from group_i_intron_designer.designer import IntronPrimerDesigner
from group_i_intron_designer.intron_analyzer import analyse_intron
from group_i_intron_designer.models import LibraryMode
from group_i_intron_designer.report import format_json_report, format_text_report
from group_i_intron_designer.sequence_utils import parse_fasta
from group_i_intron_designer.validators import validate_inputs

# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------

st.set_page_config(
    page_title="Group I Intron P1/P1ex Primer Designer",
    page_icon=":dna:",
    layout="wide",
)

st.title("Group I Intron P1/P1ex Primer Designer")
st.caption("Design site-directed mutagenesis primers for intron retargeting")

# ---------------------------------------------------------------------------
# Sidebar — inputs
# ---------------------------------------------------------------------------

with st.sidebar:
    st.header("Inputs")

    # Intron sequence: paste or upload
    intron_input_method = st.radio(
        "Intron sequence input",
        ["Paste sequence", "Upload FASTA"],
        horizontal=True,
    )

    intron_raw = ""
    if intron_input_method == "Paste sequence":
        intron_raw = st.text_area(
            "Intron DNA sequence",
            height=120,
            placeholder="Paste intron sequence here (plain DNA or FASTA)...",
        )
    else:
        uploaded = st.file_uploader("Upload FASTA file", type=["fa", "fasta", "fna", "txt"])
        if uploaded is not None:
            intron_raw = uploaded.read().decode("utf-8")

    # Exon sequences
    five_prime_exon = st.text_input(
        "5' exon sequence",
        placeholder="e.g. GAAGTCAAGT",
        help="Provide >=20 nt of the 5' exon ending at the splice site. Must end in T for a valid G-U wobble pair. The reverse primer anneals here.",
    )

    three_prime_exon = st.text_input(
        "3' exon sequence",
        placeholder="e.g. TTGAAGGTG",
        help="First ~7-12 nt of the 3' exon from the splice site.",
    )

    # Library mode
    library_mode_label = st.radio(
        "Design mode",
        [
            "Rational (single construct)",
            "Library — randomise 5' stem",
            "Library — randomise 3' stem",
            "Library — randomise both stems",
        ],
    )
    _MODE_MAP = {
        "Rational (single construct)": LibraryMode.NONE,
        "Library — randomise 5' stem": LibraryMode.FIVE_PRIME,
        "Library — randomise 3' stem": LibraryMode.THREE_PRIME,
        "Library — randomise both stems": LibraryMode.BOTH,
    }
    library_mode = _MODE_MAP[library_mode_label]

    # Advanced parameters
    with st.expander("Advanced parameters"):
        use_custom_p1 = st.checkbox("Custom P1 length", value=False)
        p1_length: int | None = None
        if use_custom_p1:
            p1_length = st.slider("P1 helix length (bp)", 4, 6, 5)

        use_custom_p10 = st.checkbox("Custom P10 length", value=False)
        p10_length: int | None = None
        if use_custom_p10:
            p10_length = st.slider("P10 helix length (bp)", 2, 5, 3)

        p10_offset = st.number_input("P10 offset (nt)", min_value=0, max_value=20, value=4)
        target_tm = st.number_input("Target annealing Tm (C)", min_value=40.0, max_value=80.0, value=60.0, step=1.0)

        st.markdown("---")
        assess_thermo = st.checkbox(
            "Thermodynamic P1/P10 scoring (ViennaRNA)",
            value=True,
            help=(
                "Score the designed P1/P10 helices for stability and "
                "composition-controlled exon specificity. Skipped if ViennaRNA "
                "is not installed."
            ),
        )
        validate_construct = st.checkbox(
            "Re-validate retargeted construct (cmscan)",
            value=False,
            help=(
                "Rebuild the retargeted intron and re-run cmscan to confirm it "
                "still classifies as the expected subtype (slower; extra CM scan)."
            ),
        )

    run_clicked = st.button("Run Design", type="primary", use_container_width=True)


# ---------------------------------------------------------------------------
# Helper: parse intron input
# ---------------------------------------------------------------------------

def _parse_intron(raw: str) -> str:
    """Extract a clean DNA sequence from raw input (plain or FASTA)."""
    text = raw.strip()
    if not text:
        return ""
    if text.startswith(">"):
        _, seq = parse_fasta(text)
        return seq
    # Plain sequence — strip whitespace and digits
    seq = re.sub(r"[\s\d]+", "", text).upper()
    return seq


# ---------------------------------------------------------------------------
# Cached Infernal analysis
# ---------------------------------------------------------------------------

@st.cache_data(show_spinner=False)
def _cached_analysis(intron_seq: str):
    """Run Infernal analysis, cached by intron sequence."""
    return analyse_intron(intron_seq)


# ---------------------------------------------------------------------------
# Main area
# ---------------------------------------------------------------------------

if not run_clicked:
    st.info(
        "Enter your intron and exon sequences in the sidebar, then click **Run Design**."
    )
    st.stop()

# Parse intron
intron_seq = _parse_intron(intron_raw)
five_prime = five_prime_exon.strip().upper()
three_prime = three_prime_exon.strip().upper()

# Validate
if not intron_seq:
    st.error("Please provide an intron sequence.")
    st.stop()
if not five_prime:
    st.error("Please provide a 5' exon sequence.")
    st.stop()
if not three_prime:
    st.error("Please provide a 3' exon sequence.")
    st.stop()

errors = validate_inputs(
    intron_seq=intron_seq,
    five_prime_exon=five_prime,
    three_prime_exon=three_prime,
    library_mode=library_mode,
    p1_length=p1_length,
    p10_length=p10_length,
    p10_offset=p10_offset,
    target_tm=target_tm,
)
if errors:
    for err in errors:
        st.error(err)
    st.stop()

# Run analysis
try:
    with st.spinner("Analysing intron structure with Infernal..."):
        analysis = _cached_analysis(intron_seq)
except FileNotFoundError as exc:
    st.error(f"Infernal not found: {exc}")
    st.stop()
except (RuntimeError, ValueError) as exc:
    st.error(f"Analysis error: {exc}")
    st.stop()

if not analysis.is_group_i:
    st.warning("Sequence may not be a group I intron.")

# Run design
try:
    designer = IntronPrimerDesigner(
        intron_sequence=intron_seq,
        five_prime_exon=five_prime,
        three_prime_exon=three_prime,
        library_mode=library_mode,
        p1_length=p1_length,
        p10_length=p10_length,
        p10_offset=p10_offset,
        target_tm=target_tm,
        analysis=analysis,
    )
    if validate_construct:
        with st.spinner("Re-validating retargeted construct with cmscan..."):
            report = designer.run(
                assess_thermodynamics=assess_thermo,
                validate_construct=True,
            )
    else:
        report = designer.run(
            assess_thermodynamics=assess_thermo,
            validate_construct=False,
        )
except Exception as exc:
    st.error(f"Design error: {exc}")
    st.stop()

# Warnings
for w in report.warnings:
    st.warning(w)

# Results in tabs
tab_report, tab_primers, tab_mutations, tab_json = st.tabs(
    ["Report", "Primers", "Mutations", "JSON"]
)

text_report = format_text_report(report)
json_report = format_json_report(report)

with tab_report:
    st.code(text_report, language=None)
    st.download_button(
        "Download text report",
        data=text_report,
        file_name="design_report.txt",
        mime="text/plain",
    )

with tab_primers:
    fwd = report.forward_primer
    rev = report.reverse_primer

    st.subheader("Forward primer")
    st.code(f"5'- {fwd.sequence} -3'", language=None)
    col1, col2, col3 = st.columns(3)
    col1.metric("Length", f"{fwd.length} nt")
    col2.metric("Annealing Tm", f"{fwd.annealing_tm:.1f} C")
    col3.metric("GC%", f"{fwd.annealing_gc_percent:.0f}%")
    if fwd.contains_degenerate:
        st.info(f"Contains {fwd.degenerate_count} degenerate (N) positions. Order with mixed bases.")

    st.subheader("Reverse primer")
    st.code(f"5'- {rev.sequence} -3'", language=None)
    col1, col2, col3 = st.columns(3)
    col1.metric("Length", f"{rev.length} nt")
    col2.metric("Annealing Tm", f"{rev.annealing_tm:.1f} C")
    col3.metric("GC%", f"{rev.annealing_gc_percent:.0f}%")

with tab_mutations:
    if report.mutations:
        rows = [
            {
                "Position": m.intron_position,
                "Original": m.original_nt,
                "New": m.new_nt,
                "Region": m.region,
            }
            for m in report.mutations
        ]
        st.dataframe(rows, use_container_width=True, hide_index=True)
    else:
        st.success("No mutations required — native intron matches the target exonic context.")

with tab_json:
    st.code(json_report, language="json")
    st.download_button(
        "Download JSON report",
        data=json_report,
        file_name="design_report.json",
        mime="application/json",
    )
