//! Hierarchical catalog structure for assemblies and distributions.
//!
//! This module provides a hierarchical catalog structure:
//! `HierarchicalCatalog` → `HierarchicalAssembly` → `AssemblyVersion` → `FastaDistribution`
//!
//! This structure allows:
//! - Clear relationships between assemblies and their versions
//! - Multiple FASTA distributions per assembly version
//! - Per-contig presence tracking (in report, in FASTA, or both)
//! - Fast lookup by MD5 or name+length

use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::path::Path;

use crate::core::assembly::{
    AssemblyVersion, FastaContig, FastaDistribution, HierarchicalAssembly,
};

/// Helper function to convert usize count to f64 with explicit precision loss allowance
#[inline]
fn count_to_f64(count: usize) -> f64 {
    #[allow(clippy::cast_precision_loss)]
    {
        count as f64
    }
}

/// Hierarchical catalog with assemblies, versions, and distributions
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HierarchicalCatalog {
    /// Catalog schema version
    pub version: String,
    /// Creation timestamp
    pub created_at: String,
    /// All assemblies in the catalog
    pub assemblies: Vec<HierarchicalAssembly>,
    /// Standalone distributions without assembly reports
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub standalone_distributions: Vec<FastaDistribution>,
}

impl Default for HierarchicalCatalog {
    fn default() -> Self {
        Self::new()
    }
}

impl HierarchicalCatalog {
    #[must_use]
    pub fn new() -> Self {
        Self {
            version: "1.0.0".to_string(),
            created_at: chrono::Utc::now().to_rfc3339(),
            assemblies: Vec::new(),
            standalone_distributions: Vec::new(),
        }
    }

    /// Add a standalone distribution (no assembly report)
    #[must_use]
    pub fn with_standalone_distribution(mut self, dist: FastaDistribution) -> Self {
        self.standalone_distributions.push(dist);
        self
    }

    /// Build indexes for fast matching
    #[must_use]
    pub fn build_index(&self) -> CatalogIndex {
        let mut index = CatalogIndex::default();

        for assembly in &self.assemblies {
            for version in &assembly.versions {
                for dist in &version.fasta_distributions {
                    for (i, contig) in dist.contigs.iter().enumerate() {
                        let location = ContigLocation {
                            assembly_id: assembly.id.clone(),
                            version_id: version.id.clone(),
                            distribution_id: dist.id.clone(),
                            contig_index: i,
                        };

                        // Index by MD5
                        if !contig.md5.is_empty() {
                            index
                                .by_md5
                                .entry(contig.md5.clone())
                                .or_default()
                                .push(location.clone());
                        }

                        // Index by name+length
                        index
                            .by_name_length
                            .entry((contig.name.clone(), contig.length))
                            .or_default()
                            .push(location);
                    }
                }
            }
        }

        // Also index standalone distributions
        for dist in &self.standalone_distributions {
            for (i, contig) in dist.contigs.iter().enumerate() {
                let location = ContigLocation {
                    assembly_id: String::new(),
                    version_id: String::new(),
                    distribution_id: dist.id.clone(),
                    contig_index: i,
                };

                if !contig.md5.is_empty() {
                    index
                        .by_md5
                        .entry(contig.md5.clone())
                        .or_default()
                        .push(location.clone());
                }

                index
                    .by_name_length
                    .entry((contig.name.clone(), contig.length))
                    .or_default()
                    .push(location);
            }
        }

        index
    }

    /// Get a distribution by ID
    #[must_use]
    pub fn get_distribution(&self, id: &str) -> Option<DistributionRef<'_>> {
        for assembly in &self.assemblies {
            for version in &assembly.versions {
                for dist in &version.fasta_distributions {
                    if dist.id == id {
                        return Some(DistributionRef {
                            assembly: Some(assembly),
                            version: Some(version),
                            distribution: dist,
                        });
                    }
                }
            }
        }

        for dist in &self.standalone_distributions {
            if dist.id == id {
                return Some(DistributionRef {
                    assembly: None,
                    version: None,
                    distribution: dist,
                });
            }
        }

        None
    }

    /// Load from JSON file
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be read or contains invalid JSON.
    pub fn load(path: &Path) -> Result<Self, std::io::Error> {
        let file = std::fs::File::open(path)?;
        let reader = std::io::BufReader::new(file);
        serde_json::from_reader(reader)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))
    }

    /// Save to JSON file
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be created or written.
    pub fn save(&self, path: &Path) -> Result<(), std::io::Error> {
        let file = std::fs::File::create(path)?;
        let writer = std::io::BufWriter::new(file);
        serde_json::to_writer_pretty(writer, self)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))
    }

    /// Infer the base assembly version for a set of contigs by MD5 matching.
    ///
    /// This method compares MD5 checksums from the input contigs against all known
    /// distributions in the catalog. Returns the best matching assembly version
    /// if the match rate exceeds the threshold.
    ///
    /// # Arguments
    /// * `contigs` - The contigs to match (must have MD5 checksums)
    /// * `min_match_rate` - Minimum fraction of contigs that must match (0.0 - 1.0)
    ///
    /// # Returns
    /// * `Some(InferredAssembly)` if a match is found above the threshold
    /// * `None` if no assembly matches well enough
    #[must_use]
    pub fn infer_base_assembly(
        &self,
        contigs: &[FastaContig],
        min_match_rate: f64,
    ) -> Option<InferredAssembly> {
        // Collect MD5s from input contigs
        let input_md5s: HashSet<&str> = contigs
            .iter()
            .filter(|c| !c.md5.is_empty())
            .map(|c| c.md5.as_str())
            .collect();

        if input_md5s.is_empty() {
            return None;
        }

        let mut best_match: Option<InferredAssembly> = None;

        // Check each assembly version's distributions
        for assembly in &self.assemblies {
            for version in &assembly.versions {
                // Collect all unique MD5s from this version's distributions
                let mut version_md5s: HashSet<&str> = HashSet::new();
                for dist in &version.fasta_distributions {
                    for contig in &dist.contigs {
                        if !contig.md5.is_empty() {
                            version_md5s.insert(&contig.md5);
                        }
                    }
                }

                if version_md5s.is_empty() {
                    continue;
                }

                // Count matches
                let matches = input_md5s.intersection(&version_md5s).count();
                let match_rate = count_to_f64(matches) / count_to_f64(input_md5s.len());

                if match_rate >= min_match_rate {
                    let is_better = best_match
                        .as_ref()
                        .map_or(true, |b| match_rate > b.match_rate);

                    if is_better {
                        best_match = Some(InferredAssembly {
                            assembly_id: assembly.id.clone(),
                            assembly_name: assembly.name.clone(),
                            version_id: version.id.clone(),
                            version_string: version.version.clone(),
                            match_rate,
                            matched_contigs: matches,
                            total_input_contigs: input_md5s.len(),
                        });
                    }
                }
            }
        }

        best_match
    }

    /// Infer base assembly with default threshold (90%)
    #[must_use]
    pub fn infer_base_assembly_default(&self, contigs: &[FastaContig]) -> Option<InferredAssembly> {
        self.infer_base_assembly(contigs, 0.9)
    }
}

/// Result of inferring the base assembly for a set of contigs
#[derive(Debug, Clone)]
pub struct InferredAssembly {
    /// Assembly ID (e.g., "grch38")
    pub assembly_id: String,
    /// Assembly name (e.g., "`GRCh38`")
    pub assembly_name: String,
    /// Version ID (e.g., "`grch38_p14`")
    pub version_id: String,
    /// Version string (e.g., "p14")
    pub version_string: String,
    /// Fraction of input contigs that matched
    pub match_rate: f64,
    /// Number of contigs that matched
    pub matched_contigs: usize,
    /// Total input contigs with MD5
    pub total_input_contigs: usize,
}

/// Reference to a distribution with its context
#[derive(Debug, Clone)]
pub struct DistributionRef<'a> {
    pub assembly: Option<&'a HierarchicalAssembly>,
    pub version: Option<&'a AssemblyVersion>,
    pub distribution: &'a FastaDistribution,
}

/// Indexes for fast lookup
#[derive(Debug, Default)]
pub struct CatalogIndex {
    /// MD5 -> locations
    pub by_md5: HashMap<String, Vec<ContigLocation>>,
    /// (name, length) -> locations
    pub by_name_length: HashMap<(String, u64), Vec<ContigLocation>>,
}

impl CatalogIndex {
    /// Find all locations for a given MD5
    #[must_use]
    pub fn find_by_md5(&self, md5: &str) -> &[ContigLocation] {
        self.by_md5.get(md5).map_or(&[], std::vec::Vec::as_slice)
    }

    /// Find all locations for a given name+length
    #[must_use]
    pub fn find_by_name_length(&self, name: &str, length: u64) -> &[ContigLocation] {
        self.by_name_length
            .get(&(name.to_string(), length))
            .map_or(&[], std::vec::Vec::as_slice)
    }
}

/// Location of a contig in the catalog
#[derive(Debug, Clone, PartialEq)]
pub struct ContigLocation {
    /// Assembly ID (empty for standalone distributions)
    pub assembly_id: String,
    /// Version ID (empty for standalone distributions)
    pub version_id: String,
    /// Distribution ID
    pub distribution_id: String,
    /// Index in distribution's contig list
    pub contig_index: usize,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::assembly::{FastaContig, ReportSource};
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
                            contigs: vec![FastaContig {
                                name: "chr1".to_string(),
                                length: 248_956_422,
                                md5: "6aef897c3d6ff0c78aff06ac189178dd".to_string(),
                                sort_order: 0,
                                report_contig_id: Some(1),
                                aliases: vec![],
                            }],
                        },
                        FastaDistribution {
                            id: "hs38DH".to_string(),
                            display_name: "hs38DH (1KG)".to_string(),
                            source: ReferenceSource::OneThousandGenomes,
                            download_url: None,
                            tags: vec![],
                            contigs: vec![FastaContig {
                                name: "chr1".to_string(),
                                length: 248_956_422,
                                md5: "6aef897c3d6ff0c78aff06ac189178dd".to_string(),
                                sort_order: 0,
                                report_contig_id: Some(1),
                                aliases: vec![],
                            }],
                        },
                    ],
                }],
            }],
            standalone_distributions: vec![],
        }
    }

    #[test]
    fn test_index_by_md5() {
        let catalog = make_test_catalog();
        let index = catalog.build_index();

        let locations = index.find_by_md5("6aef897c3d6ff0c78aff06ac189178dd");

        // chr1 appears in both hg38_ucsc and hs38DH
        assert_eq!(locations.len(), 2);
        assert!(locations.iter().any(|l| l.distribution_id == "hg38_ucsc"));
        assert!(locations.iter().any(|l| l.distribution_id == "hs38DH"));
    }

    #[test]
    fn test_index_by_name_length() {
        let catalog = make_test_catalog();
        let index = catalog.build_index();

        let locations = index.find_by_name_length("chr1", 248_956_422);

        assert_eq!(locations.len(), 2);
    }

    #[test]
    fn test_get_distribution() {
        let catalog = make_test_catalog();

        let dist_ref = catalog.get_distribution("hg38_ucsc").unwrap();
        assert_eq!(dist_ref.distribution.id, "hg38_ucsc");
        assert!(dist_ref.assembly.is_some());
        assert!(dist_ref.version.is_some());
    }

    #[test]
    fn test_save_and_load_roundtrip() {
        let catalog = make_test_catalog();
        let temp = tempfile::NamedTempFile::with_suffix(".json").unwrap();

        catalog.save(temp.path()).unwrap();
        let loaded = HierarchicalCatalog::load(temp.path()).unwrap();

        assert_eq!(catalog.version, loaded.version);
        assert_eq!(catalog.assemblies.len(), loaded.assemblies.len());
        assert_eq!(
            catalog.assemblies[0].versions[0].fasta_distributions.len(),
            loaded.assemblies[0].versions[0].fasta_distributions.len()
        );
    }

    #[test]
    fn test_standalone_distribution() {
        let mut catalog = HierarchicalCatalog::new();
        catalog.standalone_distributions.push(FastaDistribution {
            id: "custom_ref".to_string(),
            display_name: "Custom Reference".to_string(),
            source: ReferenceSource::Custom("Local".to_string()),
            download_url: None,
            tags: vec![],
            contigs: vec![FastaContig::new("contig1", 1000, "abc123")],
        });

        let index = catalog.build_index();
        let locations = index.find_by_md5("abc123");

        assert_eq!(locations.len(), 1);
        assert!(locations[0].assembly_id.is_empty());
        assert_eq!(locations[0].distribution_id, "custom_ref");
    }

    fn make_catalog_for_inference() -> HierarchicalCatalog {
        HierarchicalCatalog {
            version: "1.0.0".to_string(),
            created_at: "2024-01-01T00:00:00Z".to_string(),
            assemblies: vec![
                HierarchicalAssembly {
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
                        fasta_distributions: vec![FastaDistribution {
                            id: "hg38_ucsc".to_string(),
                            display_name: "hg38 (UCSC)".to_string(),
                            source: ReferenceSource::Ucsc,
                            download_url: None,
                            tags: vec![],
                            contigs: vec![
                                FastaContig::new("chr1", 248_956_422, "md5_grch38_chr1"),
                                FastaContig::new("chr2", 242_193_529, "md5_grch38_chr2"),
                                FastaContig::new("chr3", 198_295_559, "md5_grch38_chr3"),
                            ],
                        }],
                    }],
                },
                HierarchicalAssembly {
                    id: "grch37".to_string(),
                    name: "GRCh37".to_string(),
                    organism: "Homo sapiens".to_string(),
                    versions: vec![AssemblyVersion {
                        id: "grch37_p13".to_string(),
                        version: "p13".to_string(),
                        source: ReportSource::Ncbi {
                            accession: "GCF_000001405.25".to_string(),
                            url: None,
                            date: None,
                        },
                        report_contigs: vec![],
                        fasta_distributions: vec![FastaDistribution {
                            id: "hg19_ucsc".to_string(),
                            display_name: "hg19 (UCSC)".to_string(),
                            source: ReferenceSource::Ucsc,
                            download_url: None,
                            tags: vec![],
                            contigs: vec![
                                FastaContig::new("chr1", 249_250_621, "md5_grch37_chr1"),
                                FastaContig::new("chr2", 243_199_373, "md5_grch37_chr2"),
                            ],
                        }],
                    }],
                },
            ],
            standalone_distributions: vec![],
        }
    }

    #[test]
    fn test_infer_base_assembly_exact_match() {
        let catalog = make_catalog_for_inference();

        // Input with all GRCh38 MD5s
        let input_contigs = vec![
            FastaContig::new("chr1", 248_956_422, "md5_grch38_chr1"),
            FastaContig::new("chr2", 242_193_529, "md5_grch38_chr2"),
            FastaContig::new("chr3", 198_295_559, "md5_grch38_chr3"),
        ];

        let inferred = catalog.infer_base_assembly(&input_contigs, 0.9).unwrap();

        assert_eq!(inferred.assembly_id, "grch38");
        assert_eq!(inferred.version_id, "grch38_p14");
        assert_eq!(inferred.matched_contigs, 3);
        assert!((inferred.match_rate - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_infer_base_assembly_partial_match() {
        let catalog = make_catalog_for_inference();

        // Input with 2 out of 3 GRCh38 MD5s - but match_rate is based on input contigs
        // So 2/2 = 100% of INPUT matched, even though only 2/3 of catalog matched
        let input_contigs = vec![
            FastaContig::new("chr1", 248_956_422, "md5_grch38_chr1"),
            FastaContig::new("chr2", 242_193_529, "md5_grch38_chr2"),
        ];

        // With 90% threshold, should match since 100% of input matches
        let inferred = catalog.infer_base_assembly(&input_contigs, 0.9).unwrap();
        assert_eq!(inferred.assembly_id, "grch38");
        assert!((inferred.match_rate - 1.0).abs() < 0.001); // 100% of input matched

        // Now test with input that has non-matching contigs
        let input_with_extra = vec![
            FastaContig::new("chr1", 248_956_422, "md5_grch38_chr1"),
            FastaContig::new("chr2", 242_193_529, "md5_grch38_chr2"),
            FastaContig::new("extra_contig", 5000, "unknown_md5"),
        ];

        // 2/3 = 66.7% of input matched, below 90% threshold
        let inferred = catalog.infer_base_assembly(&input_with_extra, 0.9);
        assert!(inferred.is_none());

        // With 50% threshold, should match
        let inferred = catalog.infer_base_assembly(&input_with_extra, 0.5).unwrap();
        assert_eq!(inferred.assembly_id, "grch38");
        assert!((inferred.match_rate - 0.6667).abs() < 0.01); // ~66.7% matched
    }

    #[test]
    fn test_infer_base_assembly_no_match() {
        let catalog = make_catalog_for_inference();

        // Input with completely different MD5s
        let input_contigs = vec![
            FastaContig::new("chr1", 100_000, "unknown_md5_1"),
            FastaContig::new("chr2", 200_000, "unknown_md5_2"),
        ];

        let inferred = catalog.infer_base_assembly(&input_contigs, 0.5);
        assert!(inferred.is_none());
    }

    #[test]
    fn test_infer_base_assembly_empty_input() {
        let catalog = make_catalog_for_inference();

        // Input with no MD5s
        let input_contigs = vec![
            FastaContig::new("chr1", 100_000, ""),
            FastaContig::new("chr2", 200_000, ""),
        ];

        let inferred = catalog.infer_base_assembly(&input_contigs, 0.5);
        assert!(inferred.is_none());
    }

    #[test]
    fn test_infer_base_assembly_distinguishes_versions() {
        let catalog = make_catalog_for_inference();

        // Input with GRCh37 MD5s
        let input_contigs = vec![
            FastaContig::new("chr1", 249_250_621, "md5_grch37_chr1"),
            FastaContig::new("chr2", 243_199_373, "md5_grch37_chr2"),
        ];

        let inferred = catalog.infer_base_assembly(&input_contigs, 0.9).unwrap();

        assert_eq!(inferred.assembly_id, "grch37");
        assert_eq!(inferred.version_id, "grch37_p13");
        assert_eq!(inferred.assembly_name, "GRCh37");
    }

    #[test]
    fn test_infer_base_assembly_default_threshold() {
        let catalog = make_catalog_for_inference();

        let input_contigs = vec![
            FastaContig::new("chr1", 248_956_422, "md5_grch38_chr1"),
            FastaContig::new("chr2", 242_193_529, "md5_grch38_chr2"),
            FastaContig::new("chr3", 198_295_559, "md5_grch38_chr3"),
        ];

        let inferred = catalog.infer_base_assembly_default(&input_contigs).unwrap();
        assert_eq!(inferred.assembly_id, "grch38");
    }
}
