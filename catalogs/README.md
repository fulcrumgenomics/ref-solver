# Reference Catalog Build Guide

This document describes how to build and maintain the human reference genome catalog
used by ref-solver.

## Overview

Each reference in the catalog is built from three inputs:

1. **Sequence dictionary (`.dict`)** — provides contig names, lengths, MD5 checksums,
   and SAM metadata (AS, SP, UR, AN tags)
2. **FASTA file** — provides sha512t24u digests (and MD5 for validation)
3. **NCBI assembly report** — provides sequence roles, aliases, and cross-reference
   between naming conventions

The build script (`build_catalog.sh`) runs `ref-solver catalog build` for each
reference, then merges the per-reference JSON files into `human_references.json`.

## Reference Data Location

Reference FASTA and dict files live on scratch storage:

```text
/Volumes/scratch-00001/data/references/
  Homo_sapiens_assembly38/    # Broad GRCh38 analysis set
  GRCh38_reference_genome/    # 1000 Genomes GRCh38
  hs38DH/                     # bwakit GRCh38 + ALT + decoy + HLA
  hg38/                       # UCSC hg38
  hg38_p12_ucsc/              # UCSC hg38 with p12 patches
  GRCh38.p13/                 # NCBI GRCh38.p13
  grch38_dragen/              # Illumina DRAGEN GRCh38
  grch38_dragen_altmasked/    # Illumina DRAGEN GRCh38 ALT-masked
  grch38_gdc/                 # NCI GDC GRCh38.d1.vd1
  Homo_sapiens_assembly19/    # Broad b37
  hs37d5/                     # 1000 Genomes GRCh37 + decoy
  hg19/                       # UCSC hg19
  grch37_ncbi/                # NCBI GRCh37.p13
  chm13v2/                    # T2T-CHM13v2.0
```

## Assembly Reports

Assembly reports are stored in `catalogs/sources/` and map between naming conventions
(UCSC, RefSeq, GenBank, sequence name) for each assembly version.

### Which report for which reference

**Critical**: Each reference must use the assembly report matching its patch level.
Using the wrong patch report causes missing AN tags for patch-specific contigs.

| Report | Patch | References |
|--------|-------|------------|
| `GCF_000001405.40_GRCh38.p14_assembly_report.txt` | GRCh38.p14 | Broad, 1kg, hs38DH, hg38, DRAGEN, DRAGEN-altmasked, GDC |
| `GCF_000001405.39_GRCh38.p13_assembly_report.txt` | GRCh38.p13 | grch38_ncbi (FASTA is p13) |
| `GCF_000001405.38_GRCh38.p12_assembly_report.txt` | GRCh38.p12 | hg38_p12_ucsc |
| `GCF_000001405.25_GRCh37.p13_assembly_report.txt` | GRCh37.p13 | b37, hs37d5, hg19, grch37_ncbi |
| `GCF_009914755.1_T2T-CHM13v2.0_assembly_report.txt` | T2T-CHM13v2.0 | chm13v2 |

References without patches (Broad, 1kg, hg38, DRAGEN, etc.) can use the latest p14
report since primary assembly contigs don't change between patches.

### Additional assembly reports for non-primary contigs

| Report | Assembly | Contigs |
|--------|----------|---------|
| `GCA_000786075.2_hs38d1_assembly_report.txt` | hs38d1 decoys | KN707\* (387) and JTFH\* (1998) decoy contigs |
| `GCF_002402265.1_ASM240226v1_assembly_report.txt` | EBV | NC_007605.1 / chrEBV (1 contig) |

These are used when adding AN tags to dict files (see below) but are not directly
consumed by `ref-solver catalog build`.

## Building Dict Files

Each dict file should have complete SAM metadata on every `@SQ` line:

- **SN** — sequence name (from `samtools dict`)
- **LN** — sequence length
- **M5** — MD5 checksum
- **AS** — assembly identifier (GRCh38, GRCh37, or T2T-CHM13v2.0)
- **SP** — species (`Homo sapiens`)
- **UR** — canonical remote download URL for the FASTA
- **AN** — alternate contig names from assembly reports

### Step 1: Create base dict

```bash
samtools dict <reference.fasta> > <reference.dict>
```

This provides SN, LN, M5. The UR tag will be the local FASTA path (overridden later).

### Step 2: Set AS, SP, UR tags

Update all `@SQ` lines with the correct assembly, species, and remote URL. A simple
Python script can iterate over `@SQ` lines and replace/add the AS, SP, and UR tags.
The UR value should match the `--download-url` in `build_catalog.sh`.

| Reference | AS value |
|-----------|----------|
| All GRCh38-based | GRCh38 |
| All GRCh37-based | GRCh37 |
| T2T-CHM13v2.0 | T2T-CHM13v2.0 |

### Step 3: Add AN tags with fgbio

Use `fgbio CollectAlternateContigNames` to add alternate contig names from assembly
reports. Install fgbio via the pixi.toml in the repo root (`pixi run fgbio ...`).

```bash
pixi run fgbio CollectAlternateContigNames \
  --existing <dict> \
  --assembly-report <assembly_report> \
  --output <dict> \
  --primary <primary_column> \
  --alternates <alt1> <alt2> <alt3> \
  --sequence-roles AssembledMolecule UnlocalizedScaffold UnplacedScaffold \
                   FixPatch NovelPatch AltScaffold
```

The `--primary` flag must match the naming convention used by the dict's SN values:

| Dict naming convention | `--primary` | `--alternates` | Examples |
|----------------------|-------------|----------------|----------|
| UCSC (`chr1`) | `UcscName` | `AssignedMolecule GenBankAccession RefSeqAccession` | hg38, hg19, Broad, DRAGEN, GDC |
| Numeric (`1`) | `AssignedMolecule` | `UcscName GenBankAccession RefSeqAccession` | b37, hs37d5 |
| RefSeq (`NC_000001.11`) | `RefSeqAccession` | `AssignedMolecule GenBankAccession UcscName` | GRCh38.p13, grch37_ncbi, chm13v2 |

### Mixed naming conventions (b37, hs37d5)

The b37 and hs37d5 dicts use `AssignedMolecule` for chromosomes (`1`, `2`, ..., `MT`)
but `GenBankAccession` for scaffolds (`GL000207.1`, etc.). No single `--primary` column
matches all contigs. Use a two-pass approach:

1. Pass 1: `--primary AssignedMolecule` (picks up chromosomes: 1-22, X, Y, MT)
2. Pass 2: `--primary GenBankAccession` (picks up GL* scaffolds)

### Decoy contigs (hs38d1)

Decoy contigs (KN707\*, JTFH\*) come from the hs38d1 assembly (`GCA_000786075.2`), not
the main GRCh38 assembly. The dict SN names use UCSC-style format
(`chrUn_KN707606v1_decoy`) but the assembly report uses GenBank accessions
(`KN707606.1`). There is no direct column match, so fgbio cannot be used directly.

Instead, a Python script maps between the naming conventions:
- Strip `chrUn_` prefix and `_decoy` suffix from the SN
- Convert version format: `KN707606v1` -> `KN707606.1`
- Look up in the hs38d1 assembly report's GenBank column
- Build AN from SequenceName (`decoy00416`) and GenBank accession (`KN707606.1`)

### EBV

EBV (`chrEBV` / `NC_007605`) has its own assembly report
(`GCF_002402265.1_ASM240226v1`). For UCSC-named dicts: `AN:NC_007605.1,AJ507799.2`.
For RefSeq/numeric-named dicts: `AN:AJ507799.2,chrEBV`.

### HLA alleles

HLA contigs (`HLA-A*01:01:01:01`, etc.) come from the IPD-IMGT/HLA database. There
is no NCBI assembly report for these — they cannot get AN tags. This affects Broad,
1kg, hs38DH, and DRAGEN-altmasked (525 contigs each).

### Viral contigs (GDC only)

The GDC reference includes 199 viral sequences (HPV, CMV, HBV, HCV, HIV). These have
no NCBI assembly reports and cannot get AN tags.

## Building the Catalog

### Prerequisites

- Release binary: `cargo build --release`
- Reference FASTA/dict files on scratch storage
- Assembly reports in `catalogs/sources/`

### Run the build

```bash
./catalogs/build_catalog.sh /Volumes/scratch-00001/data/references
```

This runs `ref-solver catalog build` for each reference with:
- `--download-url` — overrides local dict UR with the canonical remote URL
- `--species "Homo sapiens"` — sets the species on all contigs
- `--assembly` — sets the assembly enum and propagates to contig assembly fields
- `--assembly-report-url` — stored as metadata on the reference

The script merges all per-reference JSON files into `catalogs/human_references.json`,
which is embedded into the binary at compile time.

### After rebuilding

```bash
cargo build --release    # Embed new catalog
cargo ci-test            # Verify tests pass
cargo ci-fmt             # Check formatting
cargo ci-lint            # Check clippy
```

## Adding a New Reference

1. Obtain the FASTA and create a dict with `samtools dict`
2. Update the dict with correct AS/SP/UR tags
3. Add AN tags using fgbio with the appropriate assembly report
4. Add a new entry to `build_catalog.sh` following the existing pattern
5. Add the output JSON filename to the `jq -s` merge command at the bottom
6. Run the build and verify with `cargo ci-test`
