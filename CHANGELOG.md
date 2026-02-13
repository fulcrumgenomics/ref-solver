# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
