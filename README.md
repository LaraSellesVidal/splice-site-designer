# Group I Intron P1/P1ex Primer Designer

A computational tool for designing site-directed mutagenesis primers to retarget group I introns to new exonic contexts. Given an intron sequence and target exon sequences, the tool redesigns the IGS (P1) and P1ex-P10 regions and outputs Q5-compatible mutagenesis primers.

## What it does

Group I introns recognise their splice sites through two short RNA helices:
- **P1 (IGS)**: pairs with the 5' exon
- **P10**: pairs with the 3' exon

To make an intron splice in a different gene, these pairings must be redesigned. This tool automates that process:

1. Classifies the intron and maps its P1/P1ex structure using [Infernal](http://eddylab.org/infernal/) covariance models
2. Designs new IGS and P10 sequences complementary to your target exons
3. Optionally randomises the P1ex internal stem for combinatorial library screening
4. **Scores the designed P1/P10 helices** for thermodynamic stability and composition-controlled exon specificity ([ViennaRNA](https://www.tbi.univie.ac.at/RNA/), optional)
5. **Re-validates the retargeted construct** against the covariance models to confirm it is still a recognisable group I intron (optional)
6. Outputs ready-to-order mutagenesis primers with melting temperatures

### Scope

The reverse-complement rule guarantees that the redesigned IGS/P10 can base-pair the new exons, and the CM re-validation confirms the catalytic core is preserved. It does not guarantee efficient splicing: the P1 helix is short and wobble-permissive, and its exon pairing does not covary above phylogenetic expectation across natural introns, so splicing competence cannot be predicted computationally. The thermodynamic scores below prioritise designs (a stronger, more exon-specific P1 is a better starting point); the functional test remains the genetic selection of the P1ex library. Rational (`--library none`) designs are best used as templates for that library rather than as guaranteed single constructs.

## Installation

### Prerequisites

- Python >= 3.10
- [Infernal](http://eddylab.org/infernal/) >= 1.1.4 (`cmscan`, `cmalign`, `cmpress` must be on PATH)
- [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) >= 2.5 (optional; enables thermodynamic P1/P10 scoring. The tool runs without it and omits those scores)

### Install the package

```bash
git clone https://github.com/YOUR_USERNAME/intron-p1ex-designer.git
cd intron-p1ex-designer
pip install -e .
```

Optional extras:
```bash
pip install -e ".[webapp]"   # Streamlit web app
pip install -e ".[thermo]"   # ViennaRNA thermodynamic scoring
pip install -e ".[all]"      # both
```

## Usage

### Web App (recommended)

```bash
streamlit run app.py
```

Open `http://localhost:8501` in your browser. Enter your intron sequence, target exon sequences, and design parameters. Results are displayed with downloadable reports.

### Command Line

```bash
# Rational design (single construct)
group-i-intron-designer \
    --intron intron.fa \
    --five-prime-exon GCTTCAGATCGCCATCGTAGCTTGAAGTCAAGT \
    --three-prime-exon TTGAAGGTG \
    --p10-offset 4 --p10-length 3

# Library mode (randomise both P1ex-stem strands)
group-i-intron-designer \
    --intron intron.fa \
    --five-prime-exon GCTTCAGATCGCCATCGTAGCTTGAAGTCAAGT \
    --three-prime-exon TTGAAGGTG \
    --library both

# With CM re-validation of the retargeted construct
group-i-intron-designer \
    --intron intron.fa \
    --five-prime-exon GCTTCAGATCGCCATCGTAGCTTGAAGTCAAGT \
    --three-prime-exon TTGAAGGTG \
    --validate-construct

# JSON output (thermodynamic scoring is on by default; --no-thermo to disable)
group-i-intron-designer \
    --intron intron.fa \
    --five-prime-exon GCTTCAGATCGCCATCGTAGCTTGAAGTCAAGT \
    --three-prime-exon TTGAAGGTG \
    --format json --output report.json
```

Additional flags:

| Flag | Default | Description |
|------|---------|-------------|
| `--no-thermo` | scoring on | Disable ViennaRNA P1/P10 thermodynamic scoring |
| `--validate-construct` | off | Rebuild the retargeted intron and re-run `cmscan` to confirm it still classifies as the expected subtype (extra CM scan) |

### Input guidelines

| Input | Description | Recommendation |
|-------|-------------|----------------|
| Intron | Full intron DNA sequence (FASTA or plain) | Must be a group I intron |
| 5' exon | Sequence ending at the splice site | **>= 20 nt** (needed for reverse primer) |
| 3' exon | Sequence starting at the splice site | >= p10_offset + p10_length nt |

The last nucleotide of the 5' exon should be **T** for a canonical G-U wobble pair at the splice site.

### Parameters

| Parameter | Default | Range | Description |
|-----------|---------|-------|-------------|
| P1 length | native | 4-6 | Base pairs in the P1 (IGS-exon) helix |
| P10 length | native | 2-5 | Base pairs in the P10 (P1ex-exon) helix |
| P10 offset | 4 | >= 0 | Nucleotides to skip in 3' exon before P10 |
| Target Tm | 60 C | 40-80 | Target annealing temperature for primers |
| Library mode | none | none/5prime/3prime/both | Which P1ex-stem strands to randomise |

### Library modes

| Mode | Description | Use case |
|------|-------------|----------|
| `none` | Rational design, single construct | Known functional context |
| `5prime` | Randomise 5' stem strand | Screen one strand |
| `3prime` | Randomise 3' stem strand | Screen complementary strand |
| `both` | Randomise both stem strands | Full combinatorial library |

## Thermodynamic assessment (optional, ViennaRNA)

When ViennaRNA is installed, each design is scored on two axes that the
reverse-complement rule alone does not capture:

- **Stability**: the minimum free energy of the intended `exon:guide` duplex
  (ViennaRNA `duplexfold`), classified `strong` / `moderate` / `weak`. GC-rich
  target exons give strong helices; AU-rich targets give weak ones that are
  harder to splice and lean more heavily on library selection.
- **Exon specificity (composition-controlled)**: the designed guide is compared
  against ~60 composition-matched shuffles of the target exon paired to the same
  guide. The reported margin is `median(random ΔG) − real ΔG` (positive means
  the real exon pairs the guide better than a composition-matched random
  sequence), and `P(random ≤ real)` is the fraction of shuffles at least as
  stable as the real pair (small means specific). This control matters because a
  short, wobble-permissive P1 can be satisfied by base composition alone for
  low-complexity exons.

As noted under [Scope](#scope), these scores prioritise designs but do not predict
splicing; the P1ex library selection is the functional test.

## Construct re-validation (optional, Infernal)

With `--validate-construct` the retargeted intron is rebuilt (mutations applied)
and re-scanned with `cmscan`. The IGS lies in the exon-pairing guide, outside the
intron's conserved covarying core, so retargeting should leave the catalytic
architecture intact. This check confirms it and flags any design whose CM match
degrades or changes subtype.

## Output

The tool produces a design report containing:
- Intron structural analysis (subtype, score, native element positions)
- New P1 and P10 helix designs with ASCII diagrams
- List of all required mutations
- Forward primer (mutagenic, may contain N for library) and reverse primer
- Primer properties (length, Tm, GC%)
- Library statistics (complexity, warnings)
- Thermodynamic assessment (P1/P10 ΔG, stability class, composition-controlled specificity margin), when ViennaRNA is available
- Construct re-validation (best CM match, bit-score, Δbits vs native), when `--validate-construct` is set
- Design warnings

## Testing

```bash
pip install pytest
pytest
```

---

The sections below document the biology, the design rules, the tool architecture,
and the validation behind the pipeline.

## Biological background

### Group I intron splicing

Group I introns are self-splicing ribozymes found across all domains of life. Their catalytic activity depends on a conserved RNA tertiary structure, but the specificity of splice-site recognition is encoded in a small region at the intron's 5' end: the **P1 helix** and its extension, the **P1ex region**.

Splicing proceeds through two transesterification steps:

1. An exogenous guanosine (exoG) attacks the 5' splice site, cleaving the 5' exon from the intron.
2. The free 3' hydroxyl of the 5' exon attacks the 3' splice site, ligating the exons and releasing the intron.

Both steps require the 5' splice site to be correctly positioned by the P1 helix and the 3' exon to be held by the P10 interaction.

### The P1/P1ex architecture

The 5' end of a group I intron forms a characteristic nested structure:

```
5' of intron, read 5' → 3':

[stem-5'] [internal loop containing P10] [stem-3'] [wobble G] [IGS] [junction] → P2
    ↑              ↑                                             ↑
    └── P1ex-stem ─┘                                        pairs with
    (internal                                               5' exon (P1)
     base-paired
     stem)

P10: pairs with 3' exon (via reverse complementarity)
IGS: pairs with 5' exon (via reverse complementarity)
```

**P1 helix (IGS):** The Internal Guide Sequence base-pairs with the last few nucleotides of the 5' exon in an antiparallel duplex. This interaction defines where the intron recognises its 5' splice site. Typical P1 length is 4-6 bp.

**P10 helix:** Located within the P1ex region, P10 base-pairs with the first few nucleotides of the 3' exon (after a short offset). This interaction positions the 3' splice site for the second step of splicing. Typical P10 length is 2-5 bp.

**P1ex internal stem:** An antiparallel stem formed between two segments within the internal loop of the P1ex region. This stem is structural rather than specificity-determining, making it a candidate for randomisation in library approaches.

**Wobble G:** A conserved guanosine on the intron side of the P1 helix that forms a non-Watson-Crick G-U wobble pair with the last nucleotide (uridine) of the 5' exon. This wobble pair is functionally essential and is maintained in all designs.

**Omega G:** The terminal guanosine at the 3' end of the intron, conserved across group I introns.

### Retargeting rationale

To make a group I intron splice in a different gene (or at a different position within the same gene), one must redesign the P1 and P10 pairings so that:
- The new IGS is complementary to the new 5' exon.
- The new P1ex-P10 is complementary to the new 3' exon.
- The conserved wobble G and omega G are preserved.
- The P1ex internal stem remains structurally viable (or is randomised to screen for functional variants).

## Design rules

The tool implements four core biological design rules.

### Rule 1: P1 (IGS) design from the 5' exon

The IGS is redesigned as the reverse complement of the target 5' exon sequence. Crucially, the **last nucleotide of the 5' exon is the wobble base** — it pairs with the conserved wobble G on the intron, not with the IGS. Therefore:

```
exon_segment = five_prime_exon[-(p1_length + 1) : -1]
new_igs = reverse_complement(exon_segment)
```

The tool checks whether the exon's terminal nucleotide is T (encoding U in RNA), which is required for a canonical G-U wobble pair. If not, a warning is issued.

### Rule 2: P1ex-P10 design from the 3' exon

P10 is redesigned as the reverse complement of a segment of the target 3' exon, offset from the splice junction:

```
exon_segment = three_prime_exon[p10_offset : p10_offset + p10_length]
new_p1ex_p10 = reverse_complement(exon_segment)
```

The default P10 offset is 4 nucleotides (matching the typical distance in natural group I introns). P10 lengths greater than 3 bp trigger a warning, as excessively strong P10 interactions can impair exon release after splicing.

### Rule 3: P1ex-stem randomisation (library modes)

The P1ex internal stem can be handled in four ways:

| Mode | 5' strand | 3' strand | Use case |
|------|-----------|-----------|----------|
| `none` (rational) | Native | Native | Single optimised construct |
| `5prime` | NNNN... | Native | Screen one strand |
| `3prime` | Native | NNNN... | Screen complementary strand |
| `both` | NNNN... | NNNN... | Full combinatorial library |

Library complexity = 4^N where N is the total number of degenerate positions (typically 4-8 positions, giving 256-65,536 variants). Primers with degenerate positions are ordered with mixed bases at N positions.

### Rule 4: Conserved boundaries

Two positions are never altered:
- **Wobble G:** Always set to G in the designed intron. The exon side should be T(U).
- **Omega G:** The terminal intron nucleotide remains G.

## Architecture

### Pipeline

```
User Input
    │
    ▼
┌─────────────────────┐
│  Input Validation    │  validators.py
│  (DNA, lengths, Tm)  │
└─────────┬───────────┘
          │
          ▼
┌─────────────────────────────────────┐
│  Intron Structural Analysis         │  intron_analyzer.py
│                                     │
│  1. cmscan → classify subtype       │  (Infernal external tool)
│  2. cmalign → Stockholm alignment   │
│  3. Parse WUSS structure            │  structural_parser.py
│     → P1ex-P10, stem, wobble G,     │
│       IGS positions                 │
└─────────┬───────────────────────────┘
          │
          ▼
┌─────────────────────────────────────┐
│  Design Engine                      │  designer.py
│                                     │
│  1. P1 design (IGS from 5' exon)    │
│  2. P1ex design (P10 from 3' exon)  │
│  3. Stem randomisation              │
│  4. Mutation enumeration            │
│  5. Library statistics              │
└─────────┬───────────────────────────┘
          │
          ▼
┌─────────────────────────────────────┐
│  Primer Construction                │  primer_designer.py
│                                     │
│  Forward: [mutations][annealing]    │
│  Reverse: RC(5' exon context)       │
│  Tm via SantaLucia 1998 NN method   │
└─────────┬───────────────────────────┘
          │
          ▼
┌─────────────────────────────────────┐
│  Report Generation                  │  report.py
│                                     │
│  - Text report with helix diagrams  │
│  - JSON for programmatic use        │
│  - Warnings and quality checks      │
└─────────────────────────────────────┘
```

### Modules

| Module | Purpose |
|--------|---------|
| `models.py` | Data classes: `IntronAnalysis`, `P1Design`, `P1exDesign`, `Mutation`, `Primer`, `LibraryInfo`, `ThermodynamicAssessment`, `DesignReport` |
| `validators.py` | Input validation (DNA characters, parameter ranges, sequence lengths) |
| `intron_analyzer.py` | Orchestrates Infernal `cmscan` + `cmalign`, returns `IntronAnalysis` |
| `structural_parser.py` | Parses Stockholm alignments and WUSS secondary structure to extract P1/P1ex element positions |
| `designer.py` | Core `IntronPrimerDesigner` class implementing all four design rules |
| `primer_designer.py` | Constructs forward and reverse Q5 SDM primers with Tm optimisation |
| `thermodynamics.py` | ViennaRNA-based P1/P10 duplex free energies + composition-controlled exon-specificity (optional; skipped if ViennaRNA absent) |
| `construct_validation.py` | Rebuilds the retargeted intron and re-runs `cmscan` to confirm it remains a valid group I intron of the expected subtype |
| `helix_drawer.py` | ASCII diagrams of P1 and P10 helices with Watson-Crick and wobble pair annotation |
| `sequence_utils.py` | DNA utilities: complement, reverse complement, GC content, nearest-neighbour Tm, FASTA parsing |
| `report.py` | Text and JSON report formatters |
| `cli.py` | Command-line interface (argparse) |
| `app.py` | Streamlit web application |

## Structural analysis via Infernal

### Classification (cmscan)

The tool ships with covariance models (CMs) for group I intron subtypes (e.g., IA1, IA2, IB, IC1, IC3, etc.). Infernal's `cmscan` is run against all available CMs to classify the input intron by subtype. The best hit (lowest E-value) determines the subtype assignment. The bit score reflects how well the sequence matches the structural profile of that subtype.

### Structural alignment (cmalign)

The input intron is then aligned to the best-matching CM using `cmalign --notrunc -g` (global, no-truncation mode). This produces a Stockholm-format alignment containing the consensus secondary structure (`SS_cons` line) in WUSS notation.

### Structural parsing

The parser (`structural_parser.py`) extracts P1/P1ex element positions through a hybrid approach:

**Step 1 — CM core identification:** The parser identifies the `<<<<<<` block at the 5' end of the alignment's SS_cons line. These positions represent the conserved outer helix recognised by the covariance model. However, the CM typically captures only a subset of the full P1/P1ex structure (often ~6 paired positions out of ~25 total).

**Step 2 — P2 boundary detection:** The parser locates the start of the next helix (P2), which defines the 3' boundary of the P1 region.

**Step 3 — Heuristic full assignment:** Because the CM core alone does not resolve the full `[stem-5'][loop+P10][stem-3'][wobbleG][IGS]` architecture, a heuristic algorithm fills in the details:

```
Algorithm: _heuristic_full_assignment

For each stem_len from max_possible down to 2:
    stem_5p = intron_seq[0 : stem_len]

    For each candidate G position in [2*stem_len, p1_region_end - 3):
        if intron[g_pos] != G: skip

        stem_3p = intron_seq[g_pos - stem_len : g_pos]

        Check antiparallel Watson-Crick complementarity:
            stem_5p[i] pairs with stem_3p[stem_len - 1 - i]

        Accept if ≥ 75% match (first pass) or ≥ 60% (fallback)

    Assign elements:
        stem_5p_positions = [0, ..., stem_len-1]
        internal_loop = [stem_len, ..., g_pos - stem_len - 1]
        P10 = last 3 positions of internal loop
        stem_3p_positions = [g_pos - stem_len, ..., g_pos - 1]
        wobble_g = g_pos
        IGS = [g_pos + 1, ..., g_pos + igs_len]
```

This algorithm was validated against the *Tetrahymena thermophila* group I intron (413 nt), where it correctly identifies:
- stem-5' at positions 0-3 (AAAT, complementary to stem-3')
- P10 at positions 13-15
- stem-3' at positions 16-19
- wobble G at position 20
- IGS at positions 21-25

## Primer design

### Forward primer (mutagenic)

The forward primer is constructed by concatenating the mutation region with a downstream annealing tail:

```
5'─[P1ex-P10][stem-5'][stem-3'][G][IGS][annealing tail]─3'
    ├──────── mutation region ────────┤├── native ─────┤
```

In the mutation region:
- P1ex-P10 positions carry the new (designed) P10 sequence.
- Stem positions carry native sequence (rational mode) or N (library modes).
- Wobble G is always G.
- IGS positions carry the new (designed) IGS sequence.

The annealing tail extends downstream into native intron sequence until the nearest-neighbour Tm reaches the target (default 60 C). Minimum annealing length is 18 nt; maximum total primer length is 60 nt.

In library modes, degenerate (N) positions are included at the appropriate stem positions. The resulting primer must be ordered with mixed bases at these positions.

### Reverse primer (non-mutagenic)

The reverse primer anneals to the opposite strand, immediately upstream of the forward primer. Since the forward primer starts at intron position 0, the reverse primer anneals to the **5' exon** sequence:

```
upstream_context = five_prime_exon + intron_seq[:first_modified_pos]
annealing_sense = extend_upstream(upstream_context, target_tm, min=18nt, max=35nt)
reverse_primer = reverse_complement(annealing_sense)
```

This produces a non-overlapping, back-to-back primer pair suitable for Q5 site-directed mutagenesis (inverse PCR). The reverse primer never contains degenerate positions.

### Melting temperature calculation

All Tm calculations use the SantaLucia (1998) unified nearest-neighbour model with:
- 16 dinucleotide enthalpy/entropy parameters
- Initiation correction for terminal A/T vs G/C
- Salt correction: Tm_corrected = Tm + 16.6 * log10([Na+])
- Default conditions: 250 nM oligo, 50 mM Na+

## Validation rules and warnings

### Input validation

| Parameter | Constraint | Rationale |
|-----------|-----------|-----------|
| All sequences | A/T/G/C only | Standard DNA alphabet |
| P1 length | 4-6 bp | Biological range for functional P1 helices |
| P10 length | 2-5 bp | Too short = no pairing; too long = impaired exon release |
| P10 offset | >= 0 | Nucleotides skipped in 3' exon |
| 5' exon | >= p1_length + 1 nt | Need p1_length for IGS + 1 wobble base |
| 3' exon | >= p10_offset + p10_length nt | Must cover the P10 pairing region |
| Target Tm | 40-80 C | Practical PCR range |
| Intron | >= 50 nt | Group I introns are typically 200+ nt |

### Design warnings

| Condition | Warning |
|-----------|---------|
| 5' exon last nt != T | G-U wobble pair disrupted |
| P10 length > 3 bp | Strong P10 may impair exon release |
| No mutations needed | Native intron already matches target context |
| Forward primer > 60 nt | Synthesis constraints |
| Forward primer Tm < 55 C | Low annealing temperature |
| Reverse primer < 18 nt | Provide longer 5' exon (>= 20 nt) |
| Reverse primer Tm < 55 C | Provide longer 5' exon (>= 20 nt) |
| Infernal score < 20 bits | May not be a canonical group I intron |
| Library complexity > 1M | Ensure sufficient transformation efficiency |
| Single-strand randomisation | Consider randomising both strands |

## Validation against Tetrahymena

The tool includes 153 unit tests covering:
- Sequence utilities (complement, reverse complement, GC content, Tm, FASTA parsing)
- Input validation (DNA characters, parameter ranges, edge cases)
- Structural parser (Stockholm parsing, WUSS notation, heuristic assignment)
- Designer (P1 design, P1ex design, mutations, library modes, wobble pair handling)
- Primer construction (forward/reverse, degenerate positions, annealing to exon)
- Integration tests (full pipeline from analysis mock to report)
- Report formatting (text and JSON)
- Thermodynamic scoring, and its behaviour without ViennaRNA (15 tests), including a check
  that enabling scoring does not change the design output

All designs were verified against manually computed results for the *Tetrahymena thermophila* group I intron (413 nt, subtype IC1) across three independent exonic contexts:

| 5' Exon | 3' Exon | Expected IGS | Tool IGS | Match |
|---------|---------|-------------|----------|-------|
| GAAGTCAAGT | TTGAAGGTG | CTTGA | CTTGA | Yes |
| GTAGAGTGT | GAGCTCCGT | CACTC | CACTC | Yes |
| CATCTTACGGAT | ATGACAGTAAGA | TCCGT | TCCGT | Yes |

In each case, the full mutant intron sequence produced by the tool was character-by-character identical to the user's independently computed expected sequence.

## License

MIT
