use serde::{Deserialize, Serialize};

/// Unique identifier for a reference in the catalog
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct ReferenceId(pub String);

impl ReferenceId {
    pub fn new(s: impl Into<String>) -> Self {
        Self(s.into())
    }
}

impl std::fmt::Display for ReferenceId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// Source organization for a reference genome
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ReferenceSource {
    Ucsc,
    Ncbi,
    Broad,
    Ensembl,
    OneThousandGenomes,
    Illumina,
    Custom(String),
}

impl std::fmt::Display for ReferenceSource {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Ucsc => write!(f, "UCSC"),
            Self::Ncbi => write!(f, "NCBI/GRC"),
            Self::Broad => write!(f, "Broad Institute"),
            Self::Ensembl => write!(f, "Ensembl"),
            Self::OneThousandGenomes => write!(f, "1000 Genomes"),
            Self::Illumina => write!(f, "Illumina"),
            Self::Custom(name) => write!(f, "{name}"),
        }
    }
}

/// Assembly version (e.g., `GRCh37`, `GRCh38`)
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum Assembly {
    Grch37,
    Grch38,
    Other(String),
}

impl std::fmt::Display for Assembly {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Grch37 => write!(f, "GRCh37"),
            Self::Grch38 => write!(f, "GRCh38"),
            Self::Other(name) => write!(f, "{name}"),
        }
    }
}

/// Naming convention used for contigs
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum NamingConvention {
    /// UCSC style: chr1, chr2, ..., chrX, chrY, chrM
    Ucsc,
    /// NCBI/Ensembl style: 1, 2, ..., X, Y, MT
    Ncbi,
    /// Mixed or unknown
    Mixed,
}

/// Type of match found
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MatchType {
    /// All contigs match exactly (same name, length, MD5)
    Exact,
    /// All MD5s match but order differs
    Reordered,
    /// All lengths match, names follow known transformation
    Renamed,
    /// Reordered and renamed
    ReorderedAndRenamed,
    /// Most contigs match, some differences
    Partial,
    /// Contigs appear to come from multiple references
    Mixed,
    /// No good match found
    NoMatch,
}

/// Confidence level for a match
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Confidence {
    Low,
    Medium,
    High,
    Exact,
}

impl Confidence {
    #[must_use]
    pub fn from_score(score: f64) -> Self {
        if score >= 1.0 {
            Self::Exact
        } else if score >= 0.95 {
            Self::High
        } else if score >= 0.80 {
            Self::Medium
        } else {
            Self::Low
        }
    }
}
