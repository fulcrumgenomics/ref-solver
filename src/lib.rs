//! # ref-solver
//!
//! A library for identifying human reference genomes from BAM/SAM/CRAM headers.
//!
//! When working with alignment files from external sources, it's often unclear exactly
//! which reference genome was used. While a reference might be labeled "`GRCh38`" or "hg19",
//! there are dozens of variations with different naming conventions, contig sets, and
//! sequence versions.
//!
//! `ref-solver` solves this by matching the sequence dictionary from your alignment file
//! against a catalog of known human reference genomes.
//!
//! ## Features
//!
//! - **MD5-based matching**: Uses sequence checksums for exact identification
//! - **Fuzzy matching**: Falls back to name+length matching when MD5s are missing
//! - **Rename detection**: Identifies when files differ only in contig naming
//! - **Order detection**: Detects when contigs are reordered vs. reference
//! - **Conflict detection**: Identifies problematic differences
//! - **Actionable suggestions**: Provides commands to fix issues
//!
//! ## Example
//!
//! ```rust,no_run
//! use ref_solver::{ReferenceCatalog, MatchingEngine, QueryHeader};
//! use ref_solver::parsing::sam::parse_header_text;
//!
//! // Load the embedded catalog of known references
//! let catalog = ReferenceCatalog::load_embedded().unwrap();
//!
//! // Parse a SAM header
//! let header_text = "@SQ\tSN:chr1\tLN:248_956_422\tM5:6aef897c3d6ff0c78aff06ac189178dd\n";
//! let query = parse_header_text(header_text).unwrap();
//!
//! // Find matching references
//! let engine = MatchingEngine::new(&catalog);
//! let matches = engine.find_matches(&query, 5);
//!
//! for m in matches {
//!     println!("{}: {:.1}%", m.reference.display_name, m.score.composite * 100.0);
//! }
//! ```
//!
//! ## Modules
//!
//! - [`catalog`]: Reference catalog storage and indexing
//! - [`core`]: Core data types for contigs, references, and headers
//! - [`matching`]: Matching engine and scoring algorithms
//! - [`parsing`]: Parsers for SAM/BAM/CRAM, dict, and TSV files
//! - [`cli`]: Command-line interface implementation
//! - [`web`]: Web server for browser-based identification

pub mod catalog;
pub mod cli;
pub mod core;
pub mod matching;
pub mod parsing;
pub mod utils;
pub mod web;

// Re-export commonly used types for convenience
pub use catalog::store::ReferenceCatalog;
pub use core::contig::Contig;
pub use core::header::QueryHeader;
pub use core::reference::KnownReference;
pub use core::types::*;
pub use matching::engine::{MatchResult, MatchingEngine};
