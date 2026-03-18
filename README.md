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
4. Outputs ready-to-order mutagenesis primers with melting temperatures

## Installation

### Prerequisites

- Python >= 3.10
- [Infernal](http://eddylab.org/infernal/) >= 1.1.4 (`cmscan`, `cmalign`, `cmpress` must be on PATH)

### Install the package

```bash
git clone https://github.com/YOUR_USERNAME/intron-p1ex-designer.git
cd intron-p1ex-designer
pip install -e .
```

For the web app:
```bash
pip install -e ".[webapp]"
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

# JSON output
group-i-intron-designer \
    --intron intron.fa \
    --five-prime-exon GCTTCAGATCGCCATCGTAGCTTGAAGTCAAGT \
    --three-prime-exon TTGAAGGTG \
    --format json --output report.json
```

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

## Design rules

1. **P1 (IGS)**: Reverse complement of the 5' exon segment upstream of the wobble base
2. **P10**: Reverse complement of the 3' exon segment (after offset)
3. **Wobble G**: Conserved G on the intron always maintained (pairs with exon U)
4. **P1ex-stem**: Native (rational) or degenerate N (library modes)

## Output

The tool produces a design report containing:
- Intron structural analysis (subtype, score, native element positions)
- New P1 and P10 helix designs with ASCII diagrams
- List of all required mutations
- Forward primer (mutagenic, may contain N for library) and reverse primer
- Primer properties (length, Tm, GC%)
- Library statistics (complexity, warnings)
- Design warnings

## Testing

```bash
pip install pytest
pytest
```

## License

MIT
