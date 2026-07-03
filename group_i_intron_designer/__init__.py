"""Group I Intron P1/P1ex Primer Designer.

A tool for designing site-directed mutagenesis primers to modify the IGS (P1)
and P1ex region of group I introns for splicing in user-specified exonic contexts.
"""

__version__ = "0.2.0"

from .designer import IntronPrimerDesigner
from .intron_analyzer import analyse_intron
from .models import (
    DesignReport,
    LibraryMode,
    ThermodynamicAssessment,
)
from .report import format_json_report, format_text_report
from .thermodynamics import viennarna_available

__all__ = [
    "IntronPrimerDesigner",
    "analyse_intron",
    "DesignReport",
    "LibraryMode",
    "ThermodynamicAssessment",
    "format_json_report",
    "format_text_report",
    "viennarna_available",
]
