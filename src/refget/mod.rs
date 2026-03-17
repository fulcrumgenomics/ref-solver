//! Refget integration for enriching unknown contigs with metadata from a refget server.
//!
//! After matching, contigs that don't match any known reference can be queried against
//! a refget server (e.g., EBI's ENA CRAM server) to retrieve aliases and other metadata.
//! This helps users identify what unknown contigs actually are.

pub mod enrichment;

use std::time::Duration;

/// Default refget server URL (EBI's ENA CRAM refget endpoint).
pub const DEFAULT_REFGET_SERVER: &str = "https://www.ebi.ac.uk/ena/cram";

/// Default per-request timeout for refget lookups.
const DEFAULT_TIMEOUT_SECS: u64 = 5;

/// Default maximum number of concurrent refget requests.
const DEFAULT_MAX_CONCURRENT: usize = 5;

/// Configuration for refget server lookups.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct RefgetConfig {
    /// Base URL of the refget server.
    pub server_url: String,
    /// Timeout for each individual HTTP request.
    pub timeout: Duration,
    /// Maximum number of concurrent requests to the refget server.
    pub max_concurrent: usize,
}

impl RefgetConfig {
    /// Create a new `RefgetConfig` with the given server URL and default settings.
    #[must_use]
    pub fn new(server_url: &str) -> Self {
        Self {
            server_url: server_url.to_string(),
            timeout: Duration::from_secs(DEFAULT_TIMEOUT_SECS),
            max_concurrent: DEFAULT_MAX_CONCURRENT,
        }
    }
}

impl Default for RefgetConfig {
    fn default() -> Self {
        Self::new(DEFAULT_REFGET_SERVER)
    }
}

/// Result of looking up a single contig in refget.
#[derive(Debug, Clone, serde::Serialize)]
#[serde(tag = "status", rename_all = "snake_case")]
#[non_exhaustive]
pub enum RefgetLookupResult {
    /// Metadata was found for this contig.
    Found {
        /// Known aliases for this sequence.
        aliases: Vec<RefgetAlias>,
        /// GA4GH sha512t24u digest.
        sha512t24u: String,
        /// Whether the sequence is circular (e.g., mitochondrial).
        circular: bool,
    },
    /// No metadata found for this digest.
    NotFound,
    /// An error occurred during lookup.
    Error {
        /// Description of the error.
        message: String,
    },
}

/// A naming-authority alias for a sequence, as returned by refget.
#[derive(Debug, Clone, serde::Serialize)]
#[non_exhaustive]
pub struct RefgetAlias {
    /// The naming authority (e.g., "insdc", "ensembl").
    pub naming_authority: String,
    /// The identifier value within that authority.
    pub value: String,
}

/// A contig enriched with optional refget metadata.
#[derive(Debug, Clone, serde::Serialize)]
#[non_exhaustive]
pub struct EnrichedContig {
    /// Name of the contig.
    pub name: String,
    /// MD5 digest used for the lookup, if available.
    pub md5: Option<String>,
    /// Result of the refget lookup.
    pub refget_metadata: RefgetLookupResult,
}
