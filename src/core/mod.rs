//! Core data types for reference genome identification.
//!
//! This module provides the fundamental types used throughout the library:
//!
//! - [`contig::Contig`]: Represents a single sequence/chromosome with name, length, and optional MD5
//! - [`header::QueryHeader`]: A sequence dictionary extracted from a BAM/SAM/CRAM file
//! - [`reference::KnownReference`]: A reference genome definition from the catalog
//! - [`types::ReferenceId`], [`types::Assembly`], [`types::ReferenceSource`]: Reference metadata types
//! - [`types::MatchType`], [`types::Confidence`]: Result classification types
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
