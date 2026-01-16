//! Hierarchical assembly types for the catalog.
//!
//! This module provides types for representing the hierarchy:
//! Assembly → `AssemblyVersion` → `FastaDistribution` → `FastaContig`
//!
//! This allows tracking:
//! - Multiple FASTA distributions per assembly version (e.g., UCSC hg38 vs 1KG hs38DH)
//! - Per-contig presence tracking (in report, in FASTA, or both)
//! - Report provenance (official NCBI vs derived from FASTA)

use serde::{Deserialize, Serialize};
use std::collections::HashSet;

use crate::core::contig::SequenceRole;
use crate::core::types::ReferenceSource;

/// Top-level assembly (e.g., `GRCh38`, `GRCh37`, CHM13)
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct HierarchicalAssembly {
    /// Unique identifier (e.g., "grch38")
    pub id: String,
    /// Display name (e.g., "`GRCh38`")
    pub name: String,
    /// Organism (e.g., "Homo sapiens")
    pub organism: String,
    /// All versions/patches of this assembly
    pub versions: Vec<AssemblyVersion>,
}

/// A specific version/patch of an assembly (has one assembly report)
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct AssemblyVersion {
    /// Unique identifier (e.g., "`grch38_p14`")
    pub id: String,
    /// Version string (e.g., "p14")
    pub version: String,
    /// Source/provenance of the assembly report
    pub source: ReportSource,
    /// Canonical contigs from assembly report
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub report_contigs: Vec<ReportContig>,
    /// FASTA distributions for this version
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub fasta_distributions: Vec<FastaDistribution>,
}

/// Provenance tracking for assembly reports
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum ReportSource {
    /// Official NCBI assembly report
    Ncbi {
        /// NCBI accession (e.g., "`GCF_000001405.40`")
        accession: String,
        /// Download URL
        #[serde(default, skip_serializing_if = "Option::is_none")]
        url: Option<String>,
        /// Report date
        #[serde(default, skip_serializing_if = "Option::is_none")]
        date: Option<String>,
    },
    /// Derived from FASTA/dict file (no official report)
    DerivedFromFasta {
        /// Source files used to build this
        source_files: Vec<String>,
        /// Inferred base assembly (e.g., "`grch38_p14`")
        #[serde(default, skip_serializing_if = "Option::is_none")]
        base_assembly: Option<String>,
    },
    /// Manually curated
    Manual {
        #[serde(default, skip_serializing_if = "Option::is_none")]
        curator: Option<String>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        notes: Option<String>,
    },
}

/// A contig from an assembly report (canonical definition)
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct ReportContig {
    /// Internal ID for linking to `FastaContig`
    pub id: u32,
    /// Primary sequence name
    pub sequence_name: String,
    /// Sequence length
    pub length: u64,
    /// MD5 checksum if known
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub md5: Option<String>,
    /// `RefSeq` accession (e.g., "`NC_000001.11`")
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub refseq_accn: Option<String>,
    /// `GenBank` accession (e.g., "CM000663.2")
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genbank_accn: Option<String>,
    /// UCSC-style name (e.g., "chr1")
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ucsc_name: Option<String>,
    /// Sequence role
    #[serde(default)]
    pub sequence_role: SequenceRole,
    /// Assigned molecule (e.g., "1", "X", "MT")
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub assigned_molecule: Option<String>,
}

/// A FASTA distribution (specific file with contigs)
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct FastaDistribution {
    /// Unique identifier (e.g., "hs38DH")
    pub id: String,
    /// Display name (e.g., "hs38DH (1KG)")
    pub display_name: String,
    /// Source organization
    pub source: ReferenceSource,
    /// Download URL
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub download_url: Option<String>,
    /// Tags (e.g., "`analysis_set`", "`with_decoy`")
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub tags: Vec<String>,
    /// Contigs in this FASTA
    pub contigs: Vec<FastaContig>,
}

impl FastaDistribution {
    /// Count contigs by presence type
    #[must_use]
    pub fn presence_counts(&self) -> PresenceCounts {
        let mut counts = PresenceCounts::default();
        for contig in &self.contigs {
            if contig.report_contig_id.is_some() {
                counts.in_both += 1;
            } else {
                counts.fasta_only += 1;
            }
        }
        counts
    }
}

/// A contig in a FASTA distribution
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct FastaContig {
    /// Contig name as it appears in the FASTA
    pub name: String,
    /// Sequence length
    pub length: u64,
    /// MD5 checksum (required for matching)
    pub md5: String,
    /// Sort order (position in original file)
    pub sort_order: u32,
    /// Link to `ReportContig` (None for decoy/HLA not in assembly report)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub report_contig_id: Option<u32>,
    /// Alternative names
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub aliases: Vec<String>,
}

impl FastaContig {
    #[cfg(test)]
    pub fn new(name: impl Into<String>, length: u64, md5: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            length,
            md5: md5.into(),
            sort_order: 0,
            report_contig_id: None,
            aliases: Vec::new(),
        }
    }

    /// Merge another contig's metadata into this one
    ///
    /// Used when building from multiple input files (e.g., dict + NCBI report)
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - Names don't match (`ContigMergeError::NameMismatch`)
    /// - Lengths don't match (`ContigMergeError::LengthMismatch`)
    /// - Both contigs have different non-empty MD5 checksums (`ContigMergeError::Md5Conflict`)
    pub fn merge(&mut self, other: &FastaContig) -> Result<(), ContigMergeError> {
        // Validate name and length match
        if self.name != other.name {
            return Err(ContigMergeError::NameMismatch {
                expected: self.name.clone(),
                found: other.name.clone(),
            });
        }
        if self.length != other.length {
            return Err(ContigMergeError::LengthMismatch {
                name: self.name.clone(),
                expected: self.length,
                found: other.length,
            });
        }

        // MD5: take first non-empty, error if conflict
        if !other.md5.is_empty() {
            if self.md5.is_empty() {
                self.md5.clone_from(&other.md5);
            } else if self.md5 != other.md5 {
                return Err(ContigMergeError::Md5Conflict {
                    name: self.name.clone(),
                    existing: self.md5.clone(),
                    incoming: other.md5.clone(),
                });
            }
        }

        // Aliases: union (excluding the contig name itself)
        let existing: HashSet<_> = self.aliases.iter().cloned().collect();
        for alias in &other.aliases {
            if !existing.contains(alias) && alias != &self.name {
                self.aliases.push(alias.clone());
            }
        }

        // report_contig_id: take first non-null
        if self.report_contig_id.is_none() {
            self.report_contig_id = other.report_contig_id;
        }

        Ok(())
    }
}

/// Counts of contigs by presence type
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct PresenceCounts {
    /// Contigs in both report and FASTA
    pub in_both: usize,
    /// Contigs only in assembly report
    pub report_only: usize,
    /// Contigs only in FASTA (decoy, HLA, etc.)
    pub fasta_only: usize,
}

/// Error when merging contig metadata
#[derive(Debug, Clone, PartialEq)]
pub enum ContigMergeError {
    NameMismatch {
        expected: String,
        found: String,
    },
    LengthMismatch {
        name: String,
        expected: u64,
        found: u64,
    },
    Md5Conflict {
        name: String,
        existing: String,
        incoming: String,
    },
}

impl std::fmt::Display for ContigMergeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NameMismatch { expected, found } => {
                write!(
                    f,
                    "Contig name mismatch: expected '{expected}', found '{found}'"
                )
            }
            Self::LengthMismatch {
                name,
                expected,
                found,
            } => {
                write!(
                    f,
                    "Length mismatch for '{name}': expected {expected}, found {found}"
                )
            }
            Self::Md5Conflict {
                name,
                existing,
                incoming,
            } => {
                write!(
                    f,
                    "MD5 conflict for '{name}': existing={existing}, incoming={incoming}"
                )
            }
        }
    }
}

impl std::error::Error for ContigMergeError {}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_contig(name: &str, length: u64, md5: &str, aliases: Vec<&str>) -> FastaContig {
        FastaContig {
            name: name.to_string(),
            length,
            md5: md5.to_string(),
            sort_order: 0,
            report_contig_id: None,
            aliases: aliases.into_iter().map(String::from).collect(),
        }
    }

    #[test]
    fn test_merge_adds_md5_when_missing() {
        let mut base = make_contig("chr1", 1000, "", vec![]);
        let other = make_contig("chr1", 1000, "abc123", vec![]);

        base.merge(&other).unwrap();

        assert_eq!(base.md5, "abc123");
    }

    #[test]
    fn test_merge_keeps_existing_md5() {
        let mut base = make_contig("chr1", 1000, "abc123", vec![]);
        let other = make_contig("chr1", 1000, "", vec![]);

        base.merge(&other).unwrap();

        assert_eq!(base.md5, "abc123");
    }

    #[test]
    fn test_merge_md5_conflict_errors() {
        let mut base = make_contig("chr1", 1000, "abc123", vec![]);
        let other = make_contig("chr1", 1000, "def456", vec![]);

        let result = base.merge(&other);

        assert!(matches!(result, Err(ContigMergeError::Md5Conflict { .. })));
    }

    #[test]
    fn test_merge_md5_same_value_ok() {
        let mut base = make_contig("chr1", 1000, "abc123", vec![]);
        let other = make_contig("chr1", 1000, "abc123", vec![]);

        base.merge(&other).unwrap();

        assert_eq!(base.md5, "abc123");
    }

    #[test]
    fn test_merge_unions_aliases() {
        let mut base = make_contig("chr1", 1000, "", vec!["1", "NC_000001"]);
        let other = make_contig("chr1", 1000, "", vec!["NC_000001", "CM000663"]);

        base.merge(&other).unwrap();

        assert_eq!(base.aliases.len(), 3);
        assert!(base.aliases.contains(&"1".to_string()));
        assert!(base.aliases.contains(&"NC_000001".to_string()));
        assert!(base.aliases.contains(&"CM000663".to_string()));
    }

    #[test]
    fn test_merge_excludes_name_from_aliases() {
        let mut base = make_contig("chr1", 1000, "", vec![]);
        let other = make_contig("chr1", 1000, "", vec!["chr1", "1"]);

        base.merge(&other).unwrap();

        assert_eq!(base.aliases, vec!["1"]);
        assert!(!base.aliases.contains(&"chr1".to_string()));
    }

    #[test]
    fn test_merge_length_mismatch_errors() {
        let mut base = make_contig("chr1", 1000, "", vec![]);
        let other = make_contig("chr1", 2000, "", vec![]);

        let result = base.merge(&other);

        assert!(matches!(
            result,
            Err(ContigMergeError::LengthMismatch { .. })
        ));
    }

    #[test]
    fn test_merge_name_mismatch_errors() {
        let mut base = make_contig("chr1", 1000, "", vec![]);
        let other = make_contig("chr2", 1000, "", vec![]);

        let result = base.merge(&other);

        assert!(matches!(result, Err(ContigMergeError::NameMismatch { .. })));
    }

    #[test]
    fn test_merge_takes_first_report_contig_id() {
        let mut base = make_contig("chr1", 1000, "", vec![]);
        base.report_contig_id = Some(42);

        let mut other = make_contig("chr1", 1000, "", vec![]);
        other.report_contig_id = Some(99);

        base.merge(&other).unwrap();

        assert_eq!(base.report_contig_id, Some(42));
    }

    #[test]
    fn test_merge_fills_missing_report_contig_id() {
        let mut base = make_contig("chr1", 1000, "", vec![]);

        let mut other = make_contig("chr1", 1000, "", vec![]);
        other.report_contig_id = Some(42);

        base.merge(&other).unwrap();

        assert_eq!(base.report_contig_id, Some(42));
    }
}
