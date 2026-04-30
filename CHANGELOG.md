# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.3.1](https://github.com/fulcrumgenomics/ref-solver/compare/v0.3.0...v0.3.1) - 2026-04-30

### Fixed

- unblock CI (rust-1.95 clippy + SHA-pinned actions) ([#19](https://github.com/fulcrumgenomics/ref-solver/pull/19))

### Other

- *(readme)* use absolute URLs for logo images ([#21](https://github.com/fulcrumgenomics/ref-solver/pull/21))

## [0.3.0](https://github.com/fulcrumgenomics/ref-solver/compare/v0.2.0...v0.3.0) - 2026-03-28

### Added

- stream binary uploads to read only headers from BAM/CRAM files ([#18](https://github.com/fulcrumgenomics/ref-solver/pull/18))
- show upload size limit and header extraction info in web UI
- extract BAM headers client-side before upload
- cap binary file uploads at 4MB for header-only parsing
- add reader-based BAM/CRAM header parsing functions
- normalize space-separated SAM headers in web UI ([#13](https://github.com/fulcrumgenomics/ref-solver/pull/13))
- add refget integration for unknown contig lookup ([#9](https://github.com/fulcrumgenomics/ref-solver/pull/9))

### Other

- replace stringly-typed error_type with enum ([#16](https://github.com/fulcrumgenomics/ref-solver/pull/16))
- update fulcrum genomics logo with light/dark theme support ([#17](https://github.com/fulcrumgenomics/ref-solver/pull/17))
- simplify headerExtractor, fix review findings
- add integration tests for reader-based BAM header parsing
- use reader-based parsing for BAM/CRAM, eliminate temp files

## [0.2.0](https://github.com/fulcrumgenomics/ref-solver/compare/v0.1.0...v0.2.0) - 2026-03-17

### Added

- add pre-commit hook for format and lint checks ([#7](https://github.com/fulcrumgenomics/ref-solver/pull/7))
- display sha512t24u digests in web UI contig details
- rebuild catalog with sha512t24u digests
- add sha512t24u (GA4GH refget) digest computation

### Fixed

- clarify web UI title to specify human genomes ([#5](https://github.com/fulcrumgenomics/ref-solver/pull/5))

### Other

- add #[non_exhaustive] to public structs and enums ([#8](https://github.com/fulcrumgenomics/ref-solver/pull/8))
- Add Zenodo DOI badge to README ([#4](https://github.com/fulcrumgenomics/ref-solver/pull/4))
- Add bioconda badge to README ([#3](https://github.com/fulcrumgenomics/ref-solver/pull/3))
- release v0.1.0 ([#1](https://github.com/fulcrumgenomics/ref-solver/pull/1))

## [0.1.0](https://github.com/fulcrumgenomics/ref-solver/releases/tag/v0.1.0) - 2026-02-13

### Other

- Prepare for first crates.io release
- Add Fulcrum Genomics branding to README and add CLAUDE.md
- Add local server deployment workflow
- Update catalog with UCSC-style aliases for patches
- Add UCSC name generation for patches in catalog builder
- Fix clippy warnings for CI
- Add documentation and unit tests for score command
- Add score command for direct file-to-file comparison
- Implement new per-contig scoring algorithm with configurable weights
- Handle all 4 name/alias matching combinations equally
- Treat alias matches equally as direct name matches
- Fix alias-based contig matching for NCBI references
- Fix composite score exceeding 100% when MD5s unavailable
- Add FAI and FASTA format support to web UI
- Set threshold to 0 when loading examples
- Add Caddy setup script for HTTPS with custom domains
- Improve web UI example data for better matching results
- Fix deploy workflow: use correct rust-toolchain action name
- Fix test data: remove underscores from integer strings
- Fix clippy pedantic warnings and add error documentation
- Add GitHub Actions workflow for Lightsail deployment
- Initial commit
