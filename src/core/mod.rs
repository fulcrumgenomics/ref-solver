//! Core data types for reference genome identification.
//!
//! This module provides the fundamental types used throughout the library:
//!
//! - [`Contig`]: Represents a single sequence/chromosome with name, length, and optional MD5
//! - [`QueryHeader`]: A sequence dictionary extracted from a BAM/SAM/CRAM file
//! - [`KnownReference`]: A reference genome definition from the catalog
//! - [`ReferenceId`], [`Assembly`], [`ReferenceSource`]: Reference metadata types
//! - [`MatchType`], [`Confidence`]: Result classification types
//!
//! ## Contig Naming
//!
//! Different reference sources use different naming conventions:
//!
//! | Source | Chromosome 1 | Mitochondrial |
//! |--------|--------------|---------------|
//! | UCSC   | chr1         | chrM          |
//! | NCBI   | 1            | MT            |
//! | Ensembl| 1            | MT            |
//!
//! Matching uses **exact names** - name equivalence is defined only through explicit
//! aliases (from SAM AN tag or NCBI assembly reports).

pub mod assembly;
pub mod contig;
pub mod header;
pub mod reference;
pub mod types;
