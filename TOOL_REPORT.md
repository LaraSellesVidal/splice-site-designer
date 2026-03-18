# Group I Intron P1/P1ex Primer Designer — Technical Report

## 1. Overview

The Group I Intron P1/P1ex Primer Designer is a computational tool for designing site-directed mutagenesis primers that retarget a group I intron to splice in a new exonic context. Given an intron sequence, a target 5' exon, and a target 3' exon, the tool:

1. Classifies the intron and maps its 5' structural elements using Infernal covariance models.
2. Designs new IGS (Internal Guide Sequence) and P1ex-P10 sequences that pair with the user-specified exons.
3. Optionally randomises the P1ex internal stem to generate a combinatorial library.
4. Outputs Q5-style mutagenesis primers with computed melting temperatures.

The tool is available as a command-line interface (`group-i-intron-designer`) and as a Streamlit web application (`app.py`).

---

## 2. Biological Background

### 2.1 Group I Intron Splicing

Group I introns are self-splicing ribozymes found across all domains of life. Their catalytic activity depends on a conserved RNA tertiary structure, but the specificity of splice-site recognition is encoded in a small region at the intron's 5' end: the **P1 helix** and its extension, the **P1ex region**.

Splicing proceeds through two transesterification steps:

1. An exogenous guanosine (exoG) attacks the 5' splice site, cleaving the 5' exon from the intron.
2. The free 3' hydroxyl of the 5' exon attacks the 3' splice site, ligating the exons and releasing the intron.

Both steps require the 5' splice site to be correctly positioned by the P1 helix and the 3' exon to be held by the P10 interaction.

### 2.2 The P1/P1ex Architecture

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

### 2.3 Retargeting Rationale

To make a group I intron splice in a different gene (or at a different position within the same gene), one must redesign the P1 and P10 pairings so that:
- The new IGS is complementary to the new 5' exon.
- The new P1ex-P10 is complementary to the new 3' exon.
- The conserved wobble G and omega G are preserved.
- The P1ex internal stem remains structurally viable (or is randomised to screen for functional variants).

---

## 3. Design Rules

The tool implements four core biological design rules:

### Rule 1: P1 (IGS) Design from the 5' Exon

The IGS is redesigned as the reverse complement of the target 5' exon sequence. Crucially, the **last nucleotide of the 5' exon is the wobble base** — it pairs with the conserved wobble G on the intron, not with the IGS. Therefore:

```
exon_segment = five_prime_exon[-(p1_length + 1) : -1]
new_igs = reverse_complement(exon_segment)
```

The tool checks whether the exon's terminal nucleotide is T (encoding U in RNA), which is required for a canonical G-U wobble pair. If not, a warning is issued.

### Rule 2: P1ex-P10 Design from the 3' Exon

P10 is redesigned as the reverse complement of a segment of the target 3' exon, offset from the splice junction:

```
exon_segment = three_prime_exon[p10_offset : p10_offset + p10_length]
new_p1ex_p10 = reverse_complement(exon_segment)
```

The default P10 offset is 4 nucleotides (matching the typical distance in natural group I introns). P10 lengths greater than 3 bp trigger a warning, as excessively strong P10 interactions can impair exon release after splicing.

### Rule 3: P1ex-Stem Randomisation (Library Modes)

The P1ex internal stem can be handled in four ways:

| Mode | 5' strand | 3' strand | Use case |
|------|-----------|-----------|----------|
| `none` (rational) | Native | Native | Single optimised construct |
| `5prime` | NNNN... | Native | Screen one strand |
| `3prime` | Native | NNNN... | Screen complementary strand |
| `both` | NNNN... | NNNN... | Full combinatorial library |

Library complexity = 4^N where N is the total number of degenerate positions (typically 4-8 positions, giving 256-65,536 variants). Primers with degenerate positions are ordered with mixed bases at N positions.

### Rule 4: Conserved Boundaries

Two positions are never altered:
- **Wobble G:** Always set to G in the designed intron. The exon side should be T(U).
- **Omega G:** The terminal intron nucleotide remains G.

---

## 4. Tool Architecture

### 4.1 Pipeline

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

### 4.2 Module Descriptions

| Module | Purpose |
|--------|---------|
| `models.py` | Data classes: `IntronAnalysis`, `P1Design`, `P1exDesign`, `Mutation`, `Primer`, `LibraryInfo`, `DesignReport` |
| `validators.py` | Input validation (DNA characters, parameter ranges, sequence lengths) |
| `intron_analyzer.py` | Orchestrates Infernal `cmscan` + `cmalign`, returns `IntronAnalysis` |
| `structural_parser.py` | Parses Stockholm alignments and WUSS secondary structure to extract P1/P1ex element positions |
| `designer.py` | Core `IntronPrimerDesigner` class implementing all four design rules |
| `primer_designer.py` | Constructs forward and reverse Q5 SDM primers with Tm optimisation |
| `helix_drawer.py` | ASCII diagrams of P1 and P10 helices with Watson-Crick and wobble pair annotation |
| `sequence_utils.py` | DNA utilities: complement, reverse complement, GC content, nearest-neighbour Tm, FASTA parsing |
| `report.py` | Text and JSON report formatters |
| `cli.py` | Command-line interface (argparse) |
| `app.py` | Streamlit web application |

---

## 5. Structural Analysis via Infernal

### 5.1 Classification (cmscan)

The tool ships with covariance models (CMs) for group I intron subtypes (e.g., IA1, IA2, IB, IC1, IC3, etc.). Infernal's `cmscan` is run against all available CMs to classify the input intron by subtype. The best hit (lowest E-value) determines the subtype assignment. The bit score reflects how well the sequence matches the structural profile of that subtype.

### 5.2 Structural Alignment (cmalign)

The input intron is then aligned to the best-matching CM using `cmalign --notrunc -g` (global, no-truncation mode). This produces a Stockholm-format alignment containing the consensus secondary structure (`SS_cons` line) in WUSS notation.

### 5.3 Structural Parsing

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

---

## 6. Primer Design

### 6.1 Forward Primer (Mutagenic)

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

### 6.2 Reverse Primer (Non-mutagenic)

The reverse primer anneals to the opposite strand, immediately upstream of the forward primer. Since the forward primer starts at intron position 0, the reverse primer anneals to the **5' exon** sequence:

```
upstream_context = five_prime_exon + intron_seq[:first_modified_pos]
annealing_sense = extend_upstream(upstream_context, target_tm, min=18nt, max=35nt)
reverse_primer = reverse_complement(annealing_sense)
```

This produces a non-overlapping, back-to-back primer pair suitable for Q5 site-directed mutagenesis (inverse PCR). The reverse primer never contains degenerate positions.

### 6.3 Melting Temperature Calculation

All Tm calculations use the SantaLucia (1998) unified nearest-neighbour model with:
- 16 dinucleotide enthalpy/entropy parameters
- Initiation correction for terminal A/T vs G/C
- Salt correction: Tm_corrected = Tm + 16.6 * log10([Na+])
- Default conditions: 250 nM oligo, 50 mM Na+

---

## 7. Validation and Warnings

### 7.1 Input Validation

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

### 7.2 Design Warnings

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

---

## 8. Output

### 8.1 Text Report

The text report contains:

1. **Intron Analysis:** Subtype classification, bit score, E-value, intron length, and a breakdown of the native 5' structure showing positions for P1ex-P10, P1ex-stem (5' and 3' strands), wobble G, and IGS.

2. **New P1 Helix:** The target 5' exon segment, new IGS (reverse complement), wobble pair status, and an ASCII helix diagram showing Watson-Crick and wobble base pairs.

3. **P1ex Design:** The target 3' exon segment (with offset), new P1ex-P10 (reverse complement), stem strand sequences (native or randomised), and an ASCII helix diagram.

4. **Library Statistics** (if applicable): Randomisation mode, number of degenerate positions per strand, total degenerate count, library complexity (4^N variants).

5. **Mutations:** A position-by-position list of all nucleotide changes required, with the region label (IGS, P1ex-P10, P1ex-stem-5prime, P1ex-stem-3prime).

6. **Primers:** Full sequences for forward and reverse primers (5' to 3'), with length, annealing Tm, GC%, and degenerate position count.

7. **Warnings:** Any design concerns (see Section 7.2).

### 8.2 JSON Report

A machine-readable JSON serialisation of the full `DesignReport` dataclass, suitable for integration into automated pipelines.

---

## 9. Interfaces

### 9.1 Command-Line Interface

```bash
group-i-intron-designer \
    --intron intron.fa \
    --five-prime-exon GCTTCAGATCGCCATCGTAGCTTGAAGTCAAGT \
    --three-prime-exon TTGAAGGTG \
    --library both \
    --p10-offset 4 \
    --p10-length 3 \
    --target-tm 60 \
    --format text \
    --output report.txt
```

### 9.2 Web Application (Streamlit)

```bash
streamlit run app.py
```

The web interface provides:
- **Sidebar inputs:** Intron sequence (paste or FASTA upload), exon sequences, library mode selector, advanced parameter controls.
- **Tabbed results:** Report (monospace text), Primers (with metrics cards), Mutations (interactive table), JSON.
- **Download buttons** for text and JSON reports.
- **Caching:** Infernal analysis is cached by intron sequence, so changing exon parameters re-runs the design instantly without re-running structural analysis.

---

## 10. Dependencies

| Dependency | Version | Purpose |
|------------|---------|---------|
| Python | >= 3.10 | Runtime |
| Infernal | >= 1.1.4 | Intron classification and structural alignment (cmscan, cmalign, cmpress) |
| Streamlit | >= 1.30 | Web application (optional, in `[project.optional-dependencies]`) |

No Python library dependencies are required beyond the standard library for the core tool. Streamlit is an optional dependency for the web interface.

---

## 11. Validation

The tool includes 138 unit tests covering:
- Sequence utilities (complement, reverse complement, GC content, Tm, FASTA parsing)
- Input validation (DNA characters, parameter ranges, edge cases)
- Structural parser (Stockholm parsing, WUSS notation, heuristic assignment)
- Designer (P1 design, P1ex design, mutations, library modes, wobble pair handling)
- Primer construction (forward/reverse, degenerate positions, annealing to exon)
- Integration tests (full pipeline from analysis mock to report)
- Report formatting (text and JSON)

All designs were verified against manually computed results for the *Tetrahymena thermophila* group I intron (413 nt, subtype IC1) across three independent exonic contexts:

| 5' Exon | 3' Exon | Expected IGS | Tool IGS | Match |
|---------|---------|-------------|----------|-------|
| GAAGTCAAGT | TTGAAGGTG | CTTGA | CTTGA | Yes |
| GTAGAGTGT | GAGCTCCGT | CACTC | CACTC | Yes |
| CATCTTACGGAT | ATGACAGTAAGA | TCCGT | TCCGT | Yes |

In each case, the full mutant intron sequence produced by the tool was character-by-character identical to the user's independently computed expected sequence.
