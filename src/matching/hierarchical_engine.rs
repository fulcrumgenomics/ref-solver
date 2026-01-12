//! Hierarchical matching engine for the new catalog structure.
//!
//! This engine matches queries against the hierarchical catalog which
//! contains Assembly → AssemblyVersion → FastaDistribution → FastaContig.

use std::collections::HashMap;

use crate::catalog::hierarchical::{CatalogIndex, HierarchicalCatalog};
use crate::core::assembly::PresenceCounts;
use crate::core::header::QueryHeader;
use crate::core::types::MatchType;

/// Helper function to convert usize count to f64 with explicit precision loss allowance
#[inline]
fn count_to_f64(count: usize) -> f64 {
    #[allow(clippy::cast_precision_loss)]
    {
        count as f64
    }
}

/// Result of matching against the hierarchical catalog
#[derive(Debug, Clone)]
pub struct HierarchicalMatchResult {
    /// The matched distribution
    pub distribution_id: String,
    /// Display name of the distribution
    pub display_name: String,
    /// Assembly ID (empty for standalone distributions)
    pub assembly_id: String,
    /// Assembly name
    pub assembly_name: String,
    /// Version ID
    pub version_id: String,
    /// Version string
    pub version_string: String,
    /// Match type
    pub match_type: MatchType,
    /// Match score (0.0 - 1.0)
    pub score: f64,
    /// Number of query contigs matched
    pub matched_contigs: usize,
    /// Total query contigs
    pub total_query_contigs: usize,
    /// Number of distribution contigs
    pub total_distribution_contigs: usize,
    /// Presence breakdown
    pub presence_counts: PresenceCounts,
    /// Extra contigs in query not in distribution
    pub extra_in_query: usize,
    /// Missing contigs from query that are in distribution
    pub missing_from_query: usize,
}

impl HierarchicalMatchResult {
    /// Get match percentage
    pub fn match_percentage(&self) -> f64 {
        self.score * 100.0
    }

    /// Check if this is an exact match
    #[cfg(test)]
    pub fn is_exact(&self) -> bool {
        matches!(self.match_type, MatchType::Exact)
    }
}

/// Matching engine for hierarchical catalogs
pub struct HierarchicalMatchingEngine<'a> {
    catalog: &'a HierarchicalCatalog,
    index: CatalogIndex,
    min_score: f64,
}

impl<'a> HierarchicalMatchingEngine<'a> {
    /// Create a new engine from a hierarchical catalog
    pub fn new(catalog: &'a HierarchicalCatalog) -> Self {
        let index = catalog.build_index();
        Self {
            catalog,
            index,
            min_score: 0.1,
        }
    }

    /// Find the best matching distributions for a query
    pub fn find_matches(&self, query: &QueryHeader, limit: usize) -> Vec<HierarchicalMatchResult> {
        // Score each distribution based on MD5 and name+length matches
        let mut scores: HashMap<String, (usize, usize)> = HashMap::new(); // dist_id -> (md5_matches, name_len_matches)

        // Count matches for each distribution
        for contig in &query.contigs {
            // Try MD5 matching first
            if let Some(md5) = &contig.md5 {
                for location in self.index.find_by_md5(md5) {
                    let entry = scores
                        .entry(location.distribution_id.clone())
                        .or_insert((0, 0));
                    entry.0 += 1;
                }
            }

            // Also try name+length matching
            for location in self.index.find_by_name_length(&contig.name, contig.length) {
                let entry = scores
                    .entry(location.distribution_id.clone())
                    .or_insert((0, 0));
                entry.1 += 1;
            }
        }

        // Build results
        let mut results: Vec<HierarchicalMatchResult> = scores
            .into_iter()
            .filter_map(|(dist_id, (md5_matches, name_len_matches))| {
                self.build_result(&dist_id, query, md5_matches, name_len_matches)
            })
            .filter(|r| r.score >= self.min_score)
            .collect();

        // Sort by score descending
        results.sort_by(|a, b| {
            b.score
                .partial_cmp(&a.score)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        results.truncate(limit);
        results
    }

    /// Build a match result for a distribution
    fn build_result(
        &self,
        dist_id: &str,
        query: &QueryHeader,
        md5_matches: usize,
        name_len_matches: usize,
    ) -> Option<HierarchicalMatchResult> {
        let dist_ref = self.catalog.get_distribution(dist_id)?;
        let dist = dist_ref.distribution;

        // Use the better of MD5 or name+length matches
        let matched = md5_matches.max(name_len_matches);

        // Calculate Jaccard similarity
        let query_size = query.contigs.len();
        let dist_size = dist.contigs.len();
        let union = query_size + dist_size - matched;
        let score = if union > 0 {
            count_to_f64(matched) / count_to_f64(union)
        } else {
            0.0
        };

        // Determine match type
        let match_type = if matched == query_size && matched == dist_size {
            MatchType::Exact
        } else if matched > 0 {
            MatchType::Partial
        } else {
            MatchType::NoMatch
        };

        let presence_counts = dist.presence_counts();

        Some(HierarchicalMatchResult {
            distribution_id: dist_id.to_string(),
            display_name: dist.display_name.clone(),
            assembly_id: dist_ref.assembly.map(|a| a.id.clone()).unwrap_or_default(),
            assembly_name: dist_ref
                .assembly
                .map(|a| a.name.clone())
                .unwrap_or_default(),
            version_id: dist_ref.version.map(|v| v.id.clone()).unwrap_or_default(),
            version_string: dist_ref
                .version
                .map(|v| v.version.clone())
                .unwrap_or_default(),
            match_type,
            score,
            matched_contigs: matched,
            total_query_contigs: query_size,
            total_distribution_contigs: dist_size,
            presence_counts,
            extra_in_query: query_size.saturating_sub(matched),
            missing_from_query: dist_size.saturating_sub(matched),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::assembly::{
        AssemblyVersion, FastaContig, FastaDistribution, HierarchicalAssembly, ReportSource,
    };
    use crate::core::contig::Contig;
    use crate::core::types::ReferenceSource;

    fn make_test_catalog() -> HierarchicalCatalog {
        HierarchicalCatalog {
            version: "1.0.0".to_string(),
            created_at: "2024-01-01T00:00:00Z".to_string(),
            assemblies: vec![HierarchicalAssembly {
                id: "grch38".to_string(),
                name: "GRCh38".to_string(),
                organism: "Homo sapiens".to_string(),
                versions: vec![AssemblyVersion {
                    id: "grch38_p14".to_string(),
                    version: "p14".to_string(),
                    source: ReportSource::Ncbi {
                        accession: "GCF_000001405.40".to_string(),
                        url: None,
                        date: None,
                    },
                    report_contigs: vec![],
                    fasta_distributions: vec![
                        FastaDistribution {
                            id: "hg38_ucsc".to_string(),
                            display_name: "hg38 (UCSC)".to_string(),
                            source: ReferenceSource::Ucsc,
                            download_url: None,
                            tags: vec![],
                            contigs: vec![
                                FastaContig::new(
                                    "chr1",
                                    248_956_422,
                                    "6aef897c3d6ff0c78aff06ac189178dd",
                                ),
                                FastaContig::new(
                                    "chr2",
                                    242_193_529,
                                    "f98db672eb0993dcfdabafe2a882905c",
                                ),
                            ],
                        },
                        FastaDistribution {
                            id: "hs38DH".to_string(),
                            display_name: "hs38DH (1KG)".to_string(),
                            source: ReferenceSource::OneThousandGenomes,
                            download_url: None,
                            tags: vec![],
                            contigs: vec![
                                FastaContig::new(
                                    "chr1",
                                    248_956_422,
                                    "6aef897c3d6ff0c78aff06ac189178dd",
                                ),
                                FastaContig::new(
                                    "chr2",
                                    242_193_529,
                                    "f98db672eb0993dcfdabafe2a882905c",
                                ),
                                FastaContig::new(
                                    "chr1_decoy",
                                    5000,
                                    "decoy_md5_hash_here_12345678",
                                ),
                            ],
                        },
                    ],
                }],
            }],
            standalone_distributions: vec![],
        }
    }

    #[test]
    fn test_exact_match() {
        let catalog = make_test_catalog();
        let engine = HierarchicalMatchingEngine::new(&catalog);

        // Query exactly matching hg38_ucsc
        let query = QueryHeader::new(vec![
            Contig::new("chr1", 248_956_422).with_md5("6aef897c3d6ff0c78aff06ac189178dd"),
            Contig::new("chr2", 242_193_529).with_md5("f98db672eb0993dcfdabafe2a882905c"),
        ]);

        let matches = engine.find_matches(&query, 5);
        assert!(!matches.is_empty());

        // Should find hg38_ucsc as exact match
        let hg38 = matches.iter().find(|m| m.distribution_id == "hg38_ucsc");
        assert!(hg38.is_some());
        let hg38 = hg38.unwrap();
        assert!(hg38.is_exact());
        assert_eq!(hg38.assembly_name, "GRCh38");
    }

    #[test]
    fn test_subset_match() {
        let catalog = make_test_catalog();
        let engine = HierarchicalMatchingEngine::new(&catalog);

        // Query with only chr1
        let query = QueryHeader::new(vec![
            Contig::new("chr1", 248_956_422).with_md5("6aef897c3d6ff0c78aff06ac189178dd")
        ]);

        let matches = engine.find_matches(&query, 5);
        assert!(!matches.is_empty());

        // Query is a partial match (subset of both distributions)
        let best = &matches[0];
        assert!(matches!(best.match_type, MatchType::Partial));
        assert_eq!(best.matched_contigs, 1);
        assert!(best.missing_from_query > 0);
    }

    #[test]
    fn test_superset_match() {
        let catalog = make_test_catalog();
        let engine = HierarchicalMatchingEngine::new(&catalog);

        // Query with extra contig not in hg38_ucsc
        let query = QueryHeader::new(vec![
            Contig::new("chr1", 248_956_422).with_md5("6aef897c3d6ff0c78aff06ac189178dd"),
            Contig::new("chr2", 242_193_529).with_md5("f98db672eb0993dcfdabafe2a882905c"),
            Contig::new("chr1_decoy", 5000).with_md5("decoy_md5_hash_here_12345678"),
        ]);

        let matches = engine.find_matches(&query, 5);
        assert!(!matches.is_empty());

        // Should be exact match for hs38DH
        let hs38dh = matches.iter().find(|m| m.distribution_id == "hs38DH");
        assert!(hs38dh.is_some());
        assert!(hs38dh.unwrap().is_exact());

        // hg38_ucsc is a partial match (query has extra contigs)
        let hg38 = matches.iter().find(|m| m.distribution_id == "hg38_ucsc");
        assert!(hg38.is_some());
        assert!(matches!(hg38.unwrap().match_type, MatchType::Partial));
        assert!(hg38.unwrap().extra_in_query > 0);
    }

    #[test]
    fn test_name_length_fallback() {
        let catalog = make_test_catalog();
        let engine = HierarchicalMatchingEngine::new(&catalog);

        // Query without MD5
        let query = QueryHeader::new(vec![
            Contig::new("chr1", 248_956_422),
            Contig::new("chr2", 242_193_529),
        ]);

        let matches = engine.find_matches(&query, 5);
        assert!(!matches.is_empty());

        // Should still find matches via name+length
        let best = &matches[0];
        assert!(best.score > 0.5);
    }

    #[test]
    fn test_assembly_info_populated() {
        let catalog = make_test_catalog();
        let engine = HierarchicalMatchingEngine::new(&catalog);

        let query = QueryHeader::new(vec![
            Contig::new("chr1", 248_956_422).with_md5("6aef897c3d6ff0c78aff06ac189178dd")
        ]);

        let matches = engine.find_matches(&query, 5);
        let best = &matches[0];

        assert_eq!(best.assembly_id, "grch38");
        assert_eq!(best.assembly_name, "GRCh38");
        assert_eq!(best.version_id, "grch38_p14");
        assert_eq!(best.version_string, "p14");
    }
}
