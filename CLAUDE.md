# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

ref-solver is a Rust bioinformatics tool that identifies which human reference genome was used to align a BAM/SAM/CRAM file. It matches the sequence dictionary (contig names, lengths, MD5 checksums) from an alignment file against a catalog of 15+ known human reference genomes (GRCh38/hg38, GRCh37/hg19, T2T-CHM13, and variants from UCSC, NCBI, Broad, DRAGEN, Illumina, 1000 Genomes).

## Build & Development Commands

```bash
cargo build                    # Debug build
cargo build --release          # Release build (LTO enabled, stripped)
cargo ci-test                  # Run tests with nextest (CI alias, --all-features)
cargo t                        # Run tests with cargo test (--all-features)
cargo ci-fmt                   # Check formatting (CI alias)
cargo ci-lint                  # Run clippy with pedantic warnings as errors (CI alias)
cargo nextest run <test_name>  # Run a single test by name
```

CI runs: `ci-test`, `ci-fmt`, and `ci-lint`. All three must pass. Clippy is configured with `-D warnings -W clippy::pedantic`.

## Architecture

The data flow is: **parsing** → **core types** → **catalog lookup** → **matching/scoring** → **CLI output or web response**.

### Core types (`src/core/`)
- `Contig`: A sequence with name, length, optional MD5 checksum, and aliases
- `QueryHeader`: The sequence dictionary extracted from an input file, with pre-computed indexes for efficient matching
- `KnownReference`: A catalog entry representing a known reference genome
- Key enums: `Assembly` (GRCh37/GRCh38/T2T_CHM13/Other), `ReferenceSource`, `NamingConvention` (UCSC/NCBI/Mixed), `MatchType`, `Confidence`

### Parsing (`src/parsing/`)
Parsers for SAM, BAM, CRAM, FASTA/FAI, dict, VCF, TSV/CSV, and NCBI assembly reports. Format is auto-detected from file extension or content. Uses the `noodles` crate for bioinformatics formats.

### Catalog (`src/catalog/`)
- `ReferenceCatalog` / `CatalogStore`: Loads the embedded catalog (`catalogs/human_references.json`, ~7.6MB) compiled into the binary at build time. Also supports custom JSON catalogs.
- `CatalogIndex`: Fast lookup via MD5, name+length, and signature matching
- `build.rs` validates catalog JSON structure at compile time; rebuild triggers on catalog changes

### Matching (`src/matching/`)
- `MatchingEngine`: Matches a `QueryHeader` against the catalog. `HierarchicalEngine` optimizes via catalog structure.
- Scoring is multi-factor: contig match quality (70%), reference coverage (20%), contig ordering (10%). Weights are configurable via `MatchingConfig`.
- Matching strategies in priority order: exact signature match (sorted MD5 hash), MD5-based Jaccard similarity, name+length fallback, alias matching
- `Diagnosis`: Generates actionable suggestions for mismatches (renaming, reordering, missing contigs)

### CLI (`src/cli/`)
Five subcommands: `identify`, `compare`, `score`, `catalog`, `serve`. Output formats: text, JSON, TSV.

### Web server (`src/web/`)
Axum-based async server with rate limiting (tower_governor), file size limits, request timeouts, and concurrency limits. Endpoints: `POST /api/identify`, `GET /api/catalog`. Includes an interactive HTML UI.

## Key Domain Concepts

- Contig naming varies across references: `chr1` (UCSC) vs `1` (NCBI) vs `NC_000001.11` (GenBank accession)
- MD5 checksums are the gold standard for sequence identity; name+length matching is a fallback
- Patch contigs (fix-patches, novel-patches) use UCSC-style naming conventions
- Mitochondrial sequence has two versions: rCRS vs old Cambridge reference sequence
