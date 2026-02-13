# ref-solver

[![Crates.io](https://img.shields.io/crates/v/ref-solver.svg)](https://crates.io/crates/ref-solver)
[![Documentation](https://docs.rs/ref-solver/badge.svg)](https://docs.rs/ref-solver)
[![CI](https://github.com/fulcrumgenomics/ref-solver/workflows/CI/badge.svg)](https://github.com/fulcrumgenomics/ref-solver/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Rust](https://img.shields.io/badge/rust-1.75%2B-blue.svg)](https://www.rust-lang.org)

Identify which human reference genome was used to align a BAM/SAM/CRAM file.

<p>
<a href="https://fulcrumgenomics.com"><img src=".github/logos/fulcrumgenomics.svg" alt="Fulcrum Genomics" height="100"/></a>
</p>

[Visit us at Fulcrum Genomics](https://www.fulcrumgenomics.com) to learn more about how we can power your Bioinformatics with ref-solver and beyond.

<a href="mailto:contact@fulcrumgenomics.com?subject=[GitHub inquiry]"><img src="https://img.shields.io/badge/Email_us-brightgreen.svg?&style=for-the-badge&logo=gmail&logoColor=white"/></a>
<a href="https://www.fulcrumgenomics.com"><img src="https://img.shields.io/badge/Visit_Us-blue.svg?&style=for-the-badge&logo=wordpress&logoColor=white"/></a>

## The Problem

When working with BAM files from external sources (collaborators, public repositories, sequencing vendors), it's often unclear exactly which reference genome was used for alignment. While the reference might be labeled "GRCh38" or "hg19", there are dozens of variations:

- **Naming conventions**: `chr1` vs `1` vs `NC_000001.11`
- **Contig sets**: With or without ALT contigs, decoys, HLA alleles
- **Mitochondrial sequences**: rCRS (16,569 bp) vs old Cambridge (16,571 bp)
- **Sources**: UCSC, NCBI, Broad, Illumina DRAGEN, 1000 Genomes

`ref-solver` solves this by matching the sequence dictionary from your BAM file against a catalog of known human reference genomes, providing:

- Exact matches when possible
- Detailed diagnostics when differences exist
- Actionable suggestions for fixing mismatches

## Quickstart

### Installation

```bash
# From crates.io (when published)
cargo install ref-solver

# From source
git clone https://github.com/fulcrumgenomics/ref-solver
cd ref-solver
cargo build --release
```

### Basic Usage

```bash
# Identify reference from a BAM file
ref-solver identify sample.bam

# From stdin (pipe samtools header)
samtools view -H sample.bam | ref-solver identify -

# JSON output for scripting
ref-solver identify sample.bam --format json

# Compare two files/references
ref-solver compare sample.bam grch38_ncbi --reference

# Score one file against another directly
ref-solver score query.bam reference.dict

# List all known references
ref-solver catalog list

# Start interactive web UI
ref-solver serve --port 8080 --open
```

### Example Output

```
#1 hg38 (UCSC) (EXACT)
   ID: hg38_ucsc
   Assembly: GRCh38
   Source: UCSC
   Match Type: Exact
   Score: 100.0%

   Contigs: 25 exact, 0 renamed, 0 by name+length, 0 unmatched, 0 conflicts

   Suggestions:
   - Safe to use as-is

   Download: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
```

## Features

- **Multiple input formats**: BAM, SAM, CRAM, Picard `.dict`, TSV
- **MD5-based matching**: Uses sequence checksums when available for exact identification
- **Fuzzy matching**: Falls back to name+length matching when MD5s are missing
- **Rename detection**: Identifies when files differ only in contig naming (chr1 vs 1)
- **Order detection**: Detects when contigs are reordered vs. reference
- **Conflict detection**: Identifies problematic differences (e.g., wrong mitochondrial sequence)
- **Actionable suggestions**: Provides commands to fix issues using fgbio/Picard tools
- **Web interface**: Interactive browser-based UI for pasting headers
- **Embedded catalog**: 15+ common human reference genomes built-in

## Supported References

The built-in catalog includes 15 human reference genomes:

### GRCh38/hg38 Family

| ID | Name | Source | Notes |
|----|------|--------|-------|
| `hg38_ucsc` | hg38 (UCSC) | UCSC | Standard UCSC reference, chr-prefixed naming |
| `grch38_ncbi` | GRCh38 (NCBI) | NCBI | NCBI numeric naming (1, 2, ..., X, Y, MT) |
| `grch38_broad_analysis_set` | GRCh38 Broad Analysis Set | Broad | GATK Best Practices reference |
| `hs38` | hs38 (no-ALT) | lh3/ref-gen | **Recommended.** Primary + unplaced/unlocalized (195 seq) |
| `hs38DH` | hs38DH | lh3/bwakit | Full set: ALT (261) + decoy (2386) + HLA (525) = 3155 seq |
| `grch38_1kg_analysis` | GRCh38 1KG Analysis Set | 1000 Genomes | 1000 Genomes Project analysis set with decoy/HLA |
| `grch38_gdc` | GRCh38.d1.vd1 | NCI GDC | TCGA/TARGET reference with 10 viral genomes (2779 seq) |
| `grch38_dragen` | GRCh38 DRAGEN | Illumina | Standard Illumina DRAGEN reference |
| `grch38_dragen_altmasked` | GRCh38 DRAGEN ALT-masked | Illumina | DRAGEN v3.9+ with ALT regions N-masked |

### GRCh37/hg19 Family

| ID | Name | Source | Notes |
|----|------|--------|-------|
| `hg19_ucsc` | hg19 (UCSC) | UCSC | ⚠️ Old Cambridge chrM (16571bp), not rCRS |
| `grch37_ncbi` | GRCh37 (NCBI) | NCBI | NCBI naming with rCRS mitochondrial |
| `hs37` | hs37 (minimal) | lh3/ref-gen | **Recommended.** 25 primary sequences only |
| `hs37d5` | hs37d5 | 1000 Genomes | With hs37d5 decoy sequence |
| `b37_broad` | b37 (Broad) | Broad | Legacy GATK Best Practices |

### T2T-CHM13

| ID | Name | Source | Notes |
|----|------|--------|-------|
| `chm13v2` | T2T-CHM13v2.0 | T2T Consortium | Complete gapless assembly, all centromeres resolved |

### Quick Reference Selection Guide

| Use Case | GRCh38 | GRCh37 |
|----------|--------|--------|
| Standard analysis (BWA-MEM, GATK) | `hs38` | `hs37` |
| ALT-aware alignment (BWA-MEM2) | `hs38DH` | `hs37d5` |
| Illumina DRAGEN | `grch38_dragen_altmasked` | — |
| TCGA/GDC compatibility | `grch38_gdc` | — |
| Legacy pipelines | `hg38_ucsc` | `hg19_ucsc` |

## Commands

### `identify`
Identify the reference genome from a BAM/SAM file.

```bash
ref-solver identify [OPTIONS] <INPUT>

Arguments:
  <INPUT>  Input file (BAM, SAM, CRAM, .dict, or TSV). Use '-' for stdin.

Options:
  -n, --max-matches <N>  Number of matches to show [default: 5]
      --exact-only       Only show exact matches
      --catalog <PATH>   Path to custom catalog file
      --input-format <FORMAT>  Override auto-detection [sam, bam, cram, dict, tsv, csv]
```

### `compare`
Compare two headers or a header against a known reference.

```bash
ref-solver compare [OPTIONS] <INPUT_A> <INPUT_B>

Arguments:
  <INPUT_A>  First input file
  <INPUT_B>  Second input file, or reference ID from catalog

Options:
      --reference        Treat INPUT_B as a reference ID from the catalog
      --catalog <PATH>   Path to custom catalog file
```

### `catalog`
Manage the reference catalog.

```bash
ref-solver catalog <COMMAND>

Commands:
  list    List all references in the catalog
  show    Show details of a specific reference
  export  Export the catalog to a JSON file
```

### `score`
Compare two files directly without using the catalog. Useful for comparing arbitrary files. By default, scoring is asymmetric: it measures how well the query matches the reference.

```bash
ref-solver score [OPTIONS] <QUERY> <REFERENCE>

Arguments:
  <QUERY>      Query file (BAM, SAM, CRAM, FASTA, FAI, VCF, .dict, TSV, CSV)
  <REFERENCE>  Reference file to compare against

Options:
      --symmetric            Compute both directions (query→reference and reference→query)
      --weight-match <N>     Weight for contig match score (0-100) [default: 70]
      --weight-coverage <N>  Weight for coverage score (0-100) [default: 20]
      --weight-order <N>     Weight for order score (0-100) [default: 10]
```

Example:
```bash
# Compare a BAM against a reference FASTA index
ref-solver score sample.bam reference.fa.fai

# Compare in both directions
ref-solver score --symmetric file_a.dict file_b.dict

# Custom scoring weights (emphasize coverage)
ref-solver score --weight-match 50 --weight-coverage 40 --weight-order 10 query.bam ref.dict
```

### `serve`
Start the web interface.

```bash
ref-solver serve [OPTIONS]

Options:
  -p, --port <PORT>      Port to listen on [default: 8080]
  -a, --address <ADDR>   Address to bind to [default: 127.0.0.1]
      --open             Open browser automatically
```

## Output Formats

Use `--format` to control output:

- `text` (default): Human-readable tabular output
- `json`: Structured JSON for scripting
- `tsv`: Tab-separated values

## Understanding Results

### Match Types

| Type | Meaning |
|------|---------|
| `Exact` | All contigs match exactly (name, length, MD5) |
| `Renamed` | Same sequences, different naming convention |
| `Reordered` | Same contigs, different order |
| `Partial` | Most contigs match, some differences |
| `Mixed` | Contigs appear to come from multiple references |
| `NoMatch` | No good match found |

### Confidence Levels

| Level | Score Range | Meaning |
|-------|-------------|---------|
| `Exact` | 100% | Perfect match |
| `High` | ≥95% | Very confident |
| `Medium` | ≥80% | Likely match |
| `Low` | <80% | Uncertain |

## Common Issues

### Wrong Mitochondrial Sequence
The UCSC hg19 `chrM` is 16,571 bp (old Cambridge sequence), while most modern references use rCRS (16,569 bp). This is a real sequence difference, not just naming.

### Chr Prefix Mismatch
UCSC uses `chr1`, NCBI/Ensembl use `1`. Use fgbio to rename:
```bash
fgbio UpdateSequenceDictionary -i input.bam -o output.bam -s reference.dict
```

### Contig Order Differences
Some tools are sensitive to contig order. Use Picard to reorder:
```bash
picard ReorderSam I=input.bam O=output.bam R=reference.fa
```

## Custom Catalogs

Export and modify the built-in catalog:

```bash
# Export current catalog
ref-solver catalog export my_catalog.json

# Edit to add custom references...

# Use custom catalog
ref-solver identify --catalog my_catalog.json sample.bam
```

## UCSC-Style Naming for Patches

ref-solver automatically generates UCSC-style names for fix-patches and novel-patches in GRCh38 assembly reports. This is particularly important for assembly reports prior to p13, where the UCSC-style-name column shows "na" for patches.

### Naming Convention

UCSC uses the following format for patch contigs:

| Patch Type | Format | Example |
|------------|--------|---------|
| Fix patches | `chr{chr}_{accession}v{version}_fix` | `chr1_KN196472v1_fix` |
| Novel patches | `chr{chr}_{accession}v{version}_alt` | `chr1_KQ458382v1_alt` |

The transformation converts NCBI GenBank accessions (e.g., `KN196472.1`) to UCSC-style names by:
1. Replacing `.` with `v` in the accession
2. Prepending `chr{chromosome}_`
3. Appending `_fix` or `_alt` based on patch type

### Examples

| NCBI Accession | Patch Type | Chromosome | UCSC Name |
|----------------|------------|------------|-----------|
| `KN196472.1` | fix-patch | 1 | `chr1_KN196472v1_fix` |
| `KQ458382.1` | novel-patch | 1 | `chr1_KQ458382v1_alt` |
| `KN196487.1` | fix-patch | Y | `chrY_KN196487v1_fix` |
| `KV766199.1` | novel-patch | X | `chrX_KV766199v1_alt` |

### Disabling UCSC Name Generation

When building custom catalogs, you can disable automatic UCSC name generation:

```bash
ref-solver catalog build \
  --id my_ref \
  --name "My Reference" \
  --input assembly_report.txt \
  --no-generate-ucsc-names
```

This is useful when you want strict adherence to names in the assembly report.

### Official Documentation

- [UCSC FAQ - Downloads](https://genome.ucsc.edu/FAQ/FAQdownloads.html)
- [UCSC Blog - Patches Explained](https://genome-blog.soe.ucsc.edu/blog/2019/02/22/patches/)
- [UCSC hg38 Downloads](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/)
- [GRC Patches Documentation](https://www.ncbi.nlm.nih.gov/grc/help/patches/)

## Library Usage

```rust
use ref_solver::{ReferenceCatalog, MatchingEngine, MatchingConfig, QueryHeader};
use ref_solver::parsing::sam::parse_header_text;

// Load catalog
let catalog = ReferenceCatalog::load_embedded()?;

// Parse a header
let header_text = "@SQ\tSN:chr1\tLN:248956422\tM5:6aef897c3d6ff0c78aff06ac189178dd\n";
let query = parse_header_text(header_text)?;

// Find matches
let engine = MatchingEngine::new(&catalog, MatchingConfig::default());
let matches = engine.find_matches(&query, 5);

for m in matches {
    println!("{}: {:.1}%", m.reference.display_name, m.score.composite * 100.0);
}
```

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Adding New References

To add a new reference to the catalog:

1. Obtain the sequence dictionary (run `samtools dict reference.fa`)
2. Add the reference to `catalogs/human_references.json`
3. Include MD5 checksums for all contigs
4. Run tests to verify matching works correctly

## License

MIT License. See [LICENSE](LICENSE) for details.

## Acknowledgments

- [Heng Li's ref-gen](https://github.com/lh3/ref-gen) for reference genome recommendations (hs37, hs38, hs38DH)
- [NCI GDC](https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files) for GRCh38.d1.vd1 reference documentation
- [Illumina DRAGEN](https://support.illumina.com/) for DRAGEN reference specifications
- [noodles](https://github.com/zaeleus/noodles) for SAM/BAM parsing
- [GATK](https://gatk.broadinstitute.org/) for reference genome documentation
- [T2T Consortium](https://github.com/marbl/CHM13) for CHM13 resources
- [1000 Genomes Project](https://www.internationalgenome.org/) for hs37d5 and analysis sets

## Citation

If you use ref-solver in your research, please cite:

```
ref-solver: A tool for identifying human reference genomes from BAM files
https://github.com/fulcrumgenomics/ref-solver
```
