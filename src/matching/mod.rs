//! Reference genome matching engine and scoring algorithms.
//!
//! This module provides the core matching functionality:
//!
//! - [`MatchingEngine`]: Main entry point for finding reference matches
//! - [`MatchScore`]: Detailed similarity scores between a query and reference
//! - [`MatchDiagnosis`]: Detailed analysis of differences and suggestions
//!
//! ## Matching Algorithm
//!
//! The matching process uses multiple strategies:
//!
//! 1. **Signature matching**: Exact match via sorted MD5 hash signature
//! 2. **MD5-based scoring**: Jaccard similarity of MD5 checksum sets
//! 3. **Name+length fallback**: When MD5s are missing, uses contig names and lengths
//! 4. **Order analysis**: Detects if contigs are reordered vs. reference
//!
//! ## Scoring
//!
//! The composite score combines multiple factors:
//!
//! - **MD5 Jaccard**: Set similarity of sequence checksums
//! - **Name+Length Jaccard**: Set similarity of (name, length) pairs
//! - **Query coverage**: Fraction of query contigs matched
//! - **Order score**: Fraction of contigs in correct relative order
//!
//! ## Example
//!
//! ```rust,no_run
//! use ref_solver::{ReferenceCatalog, MatchingEngine, MatchingConfig, QueryHeader};
//! use ref_solver::parsing::sam::parse_header_text;
//!
//! let catalog = ReferenceCatalog::load_embedded().unwrap();
//! let query = parse_header_text("@SQ\tSN:chr1\tLN:248_956_422\n").unwrap();
//!
//! let engine = MatchingEngine::new(&catalog, MatchingConfig::default());
//! let matches = engine.find_matches(&query, 5);
//!
//! for m in &matches {
//!     println!("{}: {:?} ({:.1}%)",
//!         m.reference.display_name,
//!         m.diagnosis.match_type,
//!         m.score.composite * 100.0
//!     );
//! }
//! ```

pub mod diagnosis;
pub mod engine;
pub mod hierarchical_engine;
pub mod scoring;

pub use diagnosis::Suggestion;
