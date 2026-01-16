use serde::{Deserialize, Serialize};
use std::collections::HashSet;

use crate::core::contig::{detect_naming_convention, Contig};
use crate::core::types::NamingConvention;
use crate::utils::validation::compute_signature;

/// Helper function to convert usize count to f64 with explicit precision loss allowance
#[inline]
fn count_to_f64(count: usize) -> f64 {
    #[allow(clippy::cast_precision_loss)]
    {
        count as f64
    }
}

/// A query header extracted from a BAM/SAM/CRAM file
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QueryHeader {
    /// Source file path (if known)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source: Option<String>,

    /// All contigs from @SQ lines
    pub contigs: Vec<Contig>,

    /// Detected naming convention
    pub naming_convention: NamingConvention,

    // === Pre-computed for matching ===
    /// Set of MD5s present in header
    #[serde(skip)]
    pub md5_set: HashSet<String>,

    /// Set of (`exact_name`, length) pairs for matching
    #[serde(skip)]
    pub name_length_set: HashSet<(String, u64)>,

    /// Set of (alias, length) pairs for alias-based matching
    #[serde(skip)]
    pub alias_length_set: HashSet<(String, u64)>,

    /// Signature for exact matching
    #[serde(skip)]
    pub signature: Option<String>,
}

impl QueryHeader {
    #[must_use]
    pub fn new(contigs: Vec<Contig>) -> Self {
        let naming_convention = detect_naming_convention(&contigs);

        let mut header = Self {
            source: None,
            contigs,
            naming_convention,
            md5_set: HashSet::new(),
            name_length_set: HashSet::new(),
            alias_length_set: HashSet::new(),
            signature: None,
        };

        header.rebuild_indexes();
        header
    }

    #[must_use]
    pub fn with_source(mut self, source: impl Into<String>) -> Self {
        self.source = Some(source.into());
        self
    }

    pub fn rebuild_indexes(&mut self) {
        self.md5_set.clear();
        self.name_length_set.clear();
        self.alias_length_set.clear();

        for contig in &self.contigs {
            if let Some(md5) = &contig.md5 {
                self.md5_set.insert(md5.clone());
            }
            // Use exact name for matching (no normalization)
            self.name_length_set
                .insert((contig.name.clone(), contig.length));

            // Also index aliases for reverse matching (query alias -> catalog name)
            for alias in &contig.aliases {
                self.alias_length_set.insert((alias.clone(), contig.length));
            }
        }

        // Compute signature using centralized helper
        let sig = compute_signature(&self.md5_set);
        if !sig.is_empty() {
            self.signature = Some(sig);
        }
    }

    /// Check if header has MD5 information
    #[must_use]
    pub fn has_md5s(&self) -> bool {
        !self.md5_set.is_empty()
    }

    /// Fraction of contigs with MD5 checksums
    #[must_use]
    pub fn md5_coverage(&self) -> f64 {
        if self.contigs.is_empty() {
            return 0.0;
        }
        count_to_f64(self.contigs.iter().filter(|c| c.md5.is_some()).count())
            / count_to_f64(self.contigs.len())
    }

    /// Get only primary chromosomes (1-22, X, Y, MT)
    #[cfg(test)]
    #[must_use]
    pub fn primary_contigs(&self) -> Vec<&Contig> {
        self.contigs
            .iter()
            .filter(|c| c.is_primary_chromosome() || c.is_mitochondrial())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_query_header_new() {
        let contigs = vec![Contig::new("chr1", 100), Contig::new("chr2", 200)];
        let header = QueryHeader::new(contigs);

        assert_eq!(header.contigs.len(), 2);
        assert_eq!(header.naming_convention, NamingConvention::Ucsc);
    }

    #[test]
    fn test_query_header_md5_set() {
        let contigs = vec![
            Contig::new("chr1", 100).with_md5("abc123"),
            Contig::new("chr2", 200).with_md5("def456"),
        ];
        let header = QueryHeader::new(contigs);

        assert!(header.has_md5s());
        assert_eq!(header.md5_set.len(), 2);
        assert!(header.md5_set.contains("abc123"));
        assert!(header.md5_set.contains("def456"));
    }

    #[test]
    fn test_query_header_no_md5() {
        let contigs = vec![Contig::new("chr1", 100), Contig::new("chr2", 200)];
        let header = QueryHeader::new(contigs);

        assert!(!header.has_md5s());
        assert!(header.md5_set.is_empty());
        assert!(header.signature.is_none());
    }

    #[test]
    fn test_md5_coverage() {
        let contigs = vec![
            Contig::new("chr1", 100).with_md5("abc123"),
            Contig::new("chr2", 200),
            Contig::new("chr3", 300).with_md5("ghi789"),
            Contig::new("chr4", 400),
        ];
        let header = QueryHeader::new(contigs);

        assert!((header.md5_coverage() - 0.5).abs() < 0.01);
    }

    #[test]
    fn test_md5_coverage_empty() {
        let header = QueryHeader::new(vec![]);
        assert!((header.md5_coverage() - 0.0).abs() < 0.01);
    }

    #[test]
    fn test_md5_coverage_full() {
        let contigs = vec![
            Contig::new("chr1", 100).with_md5("abc123"),
            Contig::new("chr2", 200).with_md5("def456"),
        ];
        let header = QueryHeader::new(contigs);

        assert!((header.md5_coverage() - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_name_length_set() {
        let contigs = vec![Contig::new("chr1", 100), Contig::new("chr2", 200)];
        let header = QueryHeader::new(contigs);

        assert_eq!(header.name_length_set.len(), 2);
        // Names are exact (no normalization)
        assert!(header.name_length_set.contains(&("chr1".to_string(), 100)));
        assert!(header.name_length_set.contains(&("chr2".to_string(), 200)));
    }

    #[test]
    fn test_signature_computed() {
        let contigs = vec![
            Contig::new("chr1", 100).with_md5("abc123"),
            Contig::new("chr2", 200).with_md5("def456"),
        ];
        let header = QueryHeader::new(contigs);

        assert!(header.signature.is_some());
        assert_eq!(header.signature.as_ref().unwrap().len(), 32);
    }

    #[test]
    fn test_signature_deterministic() {
        let contigs1 = vec![
            Contig::new("chr1", 100).with_md5("abc"),
            Contig::new("chr2", 200).with_md5("def"),
        ];
        let contigs2 = vec![
            Contig::new("chr2", 200).with_md5("def"),
            Contig::new("chr1", 100).with_md5("abc"),
        ];
        let header1 = QueryHeader::new(contigs1);
        let header2 = QueryHeader::new(contigs2);

        // Same MD5s in different order should produce same signature
        assert_eq!(header1.signature, header2.signature);
    }

    #[test]
    fn test_with_source() {
        let contigs = vec![Contig::new("chr1", 100)];
        let header = QueryHeader::new(contigs).with_source("/path/to/file.bam");

        assert_eq!(header.source, Some("/path/to/file.bam".to_string()));
    }

    #[test]
    fn test_primary_contigs() {
        let contigs = vec![
            Contig::new("chr1", 100),
            Contig::new("chr22", 200),
            Contig::new("chrX", 300),
            Contig::new("chrY", 400),
            Contig::new("chrM", 500),
            Contig::new("chr1_random", 600),
            Contig::new("chrUn_gl000220", 700),
        ];
        let header = QueryHeader::new(contigs);

        let primary = header.primary_contigs();
        assert_eq!(primary.len(), 5); // chr1, chr22, chrX, chrY, chrM
    }

    #[test]
    fn test_detect_ucsc_naming() {
        let contigs = vec![Contig::new("chr1", 100), Contig::new("chr2", 200)];
        let header = QueryHeader::new(contigs);
        assert_eq!(header.naming_convention, NamingConvention::Ucsc);
    }

    #[test]
    fn test_detect_ncbi_naming() {
        let contigs = vec![Contig::new("1", 100), Contig::new("2", 200)];
        let header = QueryHeader::new(contigs);
        assert_eq!(header.naming_convention, NamingConvention::Ncbi);
    }

    #[test]
    fn test_detect_mixed_naming() {
        let contigs = vec![Contig::new("chr1", 100), Contig::new("2", 200)];
        let header = QueryHeader::new(contigs);
        assert_eq!(header.naming_convention, NamingConvention::Mixed);
    }
}
