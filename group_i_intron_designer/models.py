"""Data models for the Group I Intron P1/P1ex Primer Designer."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum


class LibraryMode(Enum):
    """Which strands of the P1ex stem to randomise."""

    NONE = "none"
    FIVE_PRIME = "5prime"
    THREE_PRIME = "3prime"
    BOTH = "both"


@dataclass
class P1exStemInfo:
    """Positions and sequences for the two strands of the P1ex internal stem.

    In the linear intron sequence the 5' strand lies between the P1ex-P10
    portion and the 3' strand, which in turn precedes the wobble G.
    """

    five_prime_positions: list[int]
    three_prime_positions: list[int]
    five_prime_native: str
    three_prime_native: str
    stem_length: int


@dataclass
class IntronAnalysis:
    """Result of structural analysis of a group I intron via Infernal."""

    is_group_i: bool
    subtype: str | None
    score: float
    e_value: float
    igs_positions: list[int]
    p1ex_p10_positions: list[int]
    p1ex_stem: P1exStemInfo
    wobble_g_position: int
    omega_g_position: int
    native_igs: str
    native_p1ex_p10: str
    native_p1_length: int
    native_p1ex_p10_length: int
    intron_length: int


@dataclass
class P1Design:
    """Design of the new P1 helix (IGS pairing with the 5' exon)."""

    new_igs: str
    exon_segment: str
    wobble_pair: tuple[str, str]
    wobble_valid: bool
    helix_diagram: str
    p1_length: int


@dataclass
class P1exDesign:
    """Design of the new P1ex region (P10 + stem)."""

    new_p1ex_p10: str
    stem_5prime_sequence: str
    stem_3prime_sequence: str
    library_mode: LibraryMode
    helix_diagram: str
    p10_offset: int
    p10_length: int
    exon_segment: str


@dataclass
class Mutation:
    """A single nucleotide change in the intron."""

    intron_position: int
    original_nt: str
    new_nt: str  # 'N' for degenerate positions
    region: str  # "IGS", "P1ex-P10", "P1ex-stem-5prime", "P1ex-stem-3prime"


@dataclass
class Primer:
    """A mutagenesis primer (forward or reverse)."""

    name: str
    sequence: str  # may contain N for degenerate positions
    length: int
    mutation_region: str
    annealing_region: str
    annealing_tm: float
    annealing_gc_percent: float
    contains_degenerate: bool
    degenerate_count: int


@dataclass
class LibraryInfo:
    """Summary of the combinatorial library characteristics."""

    mode: LibraryMode
    stem_5prime_randomised: bool
    stem_3prime_randomised: bool
    n_positions_5prime: int
    n_positions_3prime: int
    total_degenerate: int
    complexity: int  # 4^total_degenerate
    warnings: list[str] = field(default_factory=list)


@dataclass
class ThermodynamicAssessment:
    """Thermodynamic scoring of the designed P1/P10 helices.

    Duplex free energies, stability labels, and composition-controlled specificity
    margins (see thermodynamics.py). Numeric fields are None if ViennaRNA is not
    installed.
    """

    viennarna_available: bool
    p1_dg: float | None
    p1_stability: str
    p1_specificity_margin: float | None
    p1_frac_random_stronger: float | None
    p10_dg: float | None
    p10_stability: str
    p10_specificity_margin: float | None
    p10_frac_random_stronger: float | None
    warnings: list[str] = field(default_factory=list)


@dataclass
class DesignReport:
    """Complete output of a design run."""

    intron_analysis: IntronAnalysis
    p1: P1Design
    p1ex: P1exDesign
    mutations: list[Mutation]
    forward_primer: Primer
    reverse_primer: Primer
    library: LibraryInfo | None  # None when mode = none
    warnings: list[str] = field(default_factory=list)
    # Optional, additive assessments (None when the feature is disabled/unavailable).
    thermodynamics: ThermodynamicAssessment | None = None
    construct_validation: "ConstructValidation | None" = None
