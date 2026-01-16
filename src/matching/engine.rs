use crate::catalog::index::CandidateFinder;
use crate::catalog::store::ReferenceCatalog;
use crate::core::header::QueryHeader;
use crate::core::reference::KnownReference;
use crate::matching::diagnosis::MatchDiagnosis;
use crate::matching::scoring::MatchScore;

/// Result of matching a query against the catalog
#[derive(Debug, Clone)]
pub struct MatchResult {
    /// The matched reference
    pub reference: KnownReference,

    /// Match score details
    pub score: MatchScore,

    /// Detailed diagnosis
    pub diagnosis: MatchDiagnosis,
}

impl MatchResult {
    #[must_use]
    pub fn new(reference: &KnownReference, query: &QueryHeader) -> Self {
        let score = MatchScore::calculate(query, reference);
        let diagnosis = MatchDiagnosis::analyze(query, reference);

        Self {
            reference: reference.clone(),
            score,
            diagnosis,
        }
    }

    #[must_use]
    pub fn new_with_weights(
        reference: &KnownReference,
        query: &QueryHeader,
        weights: &ScoringWeights,
    ) -> Self {
        let score = MatchScore::calculate_with_weights(query, reference, weights);
        let diagnosis = MatchDiagnosis::analyze(query, reference);

        Self {
            reference: reference.clone(),
            score,
            diagnosis,
        }
    }
}

/// Default minimum score threshold for matches
pub const DEFAULT_MIN_SCORE: f64 = 0.1;

/// Configuration for the matching engine
#[derive(Debug, Clone)]
pub struct MatchingConfig {
    /// Minimum score threshold for including matches in results
    pub min_score: f64,
    /// Custom scoring weights
    pub scoring_weights: ScoringWeights,
}

impl Default for MatchingConfig {
    fn default() -> Self {
        Self {
            min_score: DEFAULT_MIN_SCORE,
            scoring_weights: ScoringWeights::default(),
        }
    }
}

/// Configurable weights for the scoring algorithm
///
/// The scoring algorithm uses three main components:
/// - `contig_match`: Per-contig match quality (exact, name+length, conflicts)
/// - coverage: What fraction of the reference is covered by good matches
/// - order: Whether contigs appear in the same order
///
/// Additionally, `conflict_penalty` controls how much credit MD5 conflicts receive
/// (name+length matches but MD5 differs, indicating different sequences).
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ScoringWeights {
    /// Weight for per-contig match score (default 0.70 = 70%)
    #[serde(default = "default_contig_match")]
    pub contig_match: f64,

    /// Weight for reference coverage (default 0.20 = 20%)
    #[serde(default = "default_coverage")]
    pub coverage: f64,

    /// Weight for contig ordering (default 0.10 = 10%)
    #[serde(default = "default_order")]
    pub order: f64,

    /// Credit given to MD5 conflicts (default 0.1 = 10% credit)
    /// Set to 0.0 for zero credit on conflicts, 1.0 to treat conflicts as matches
    #[serde(default = "default_conflict_penalty")]
    pub conflict_penalty: f64,
}

fn default_contig_match() -> f64 {
    0.70
}
fn default_coverage() -> f64 {
    0.20
}
fn default_order() -> f64 {
    0.10
}
fn default_conflict_penalty() -> f64 {
    0.1
}

impl Default for ScoringWeights {
    fn default() -> Self {
        Self {
            contig_match: 0.70,    // 70%
            coverage: 0.20,        // 20%
            order: 0.10,           // 10%
            conflict_penalty: 0.1, // 10% credit for conflicts
        }
    }
}

impl ScoringWeights {
    /// Normalize the main weights (`contig_match`, coverage, order) to sum to 1.0
    /// The `conflict_penalty` is kept as-is since it's a multiplier, not a weight.
    #[must_use]
    pub fn normalized(&self) -> Self {
        let total = self.contig_match + self.coverage + self.order;

        if total <= 0.0 {
            return Self::default();
        }

        Self {
            contig_match: self.contig_match / total,
            coverage: self.coverage / total,
            order: self.order / total,
            conflict_penalty: self.conflict_penalty.clamp(0.0, 1.0),
        }
    }
}

/// The main matching engine
pub struct MatchingEngine<'a> {
    catalog: &'a ReferenceCatalog,
    /// Configuration including scoring weights and thresholds
    config: MatchingConfig,
}

impl<'a> MatchingEngine<'a> {
    /// Create a new matching engine with default configuration
    #[must_use]
    pub fn new(catalog: &'a ReferenceCatalog) -> Self {
        Self {
            catalog,
            config: MatchingConfig::default(),
        }
    }

    /// Create a new matching engine with custom configuration
    #[must_use]
    pub fn with_config(catalog: &'a ReferenceCatalog, config: MatchingConfig) -> Self {
        Self { catalog, config }
    }

    /// Find the best matching references for a query
    #[must_use]
    pub fn find_matches(&self, query: &QueryHeader, limit: usize) -> Vec<MatchResult> {
        // Step 1: Try exact signature match
        if let Some(sig) = &query.signature {
            if let Some(reference) = self.catalog.find_by_signature(sig) {
                return vec![MatchResult::new_with_weights(
                    reference,
                    query,
                    &self.config.scoring_weights,
                )];
            }
        }

        // Step 2: Find candidates via index
        let finder = CandidateFinder::new(self.catalog);
        let candidate_indices = finder.find_top_candidates(query, limit * 2);

        // Step 3: Score and rank candidates with custom weights
        let mut results: Vec<MatchResult> = candidate_indices
            .into_iter()
            .map(|idx| {
                let reference = &self.catalog.references[idx];
                MatchResult::new_with_weights(reference, query, &self.config.scoring_weights)
            })
            .collect();

        // Sort by composite score descending
        results.sort_by(|a, b| {
            b.score
                .composite
                .partial_cmp(&a.score.composite)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        // Filter to meaningful matches and limit
        results
            .into_iter()
            .filter(|r| r.score.composite > self.config.min_score)
            .take(limit)
            .collect()
    }

    /// Find the single best match
    #[cfg(test)]
    #[must_use]
    pub fn find_best_match(&self, query: &QueryHeader) -> Option<MatchResult> {
        self.find_matches(query, 1).into_iter().next()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::contig::Contig;

    fn make_test_catalog() -> ReferenceCatalog {
        ReferenceCatalog::load_embedded().unwrap()
    }

    #[test]
    fn test_find_matches_grch38() {
        let catalog = make_test_catalog();
        let engine = MatchingEngine::new(&catalog);

        // Create query with GRCh38 chromosomes
        let contigs = vec![
            Contig::new("chr1", 248_956_422).with_md5("6aef897c3d6ff0c78aff06ac189178dd"),
            Contig::new("chr2", 242_193_529).with_md5("f98db672eb0993dcfdabafe2a882905c"),
            Contig::new("chr3", 198_295_559).with_md5("76635a41ea913a405ded820447d067b0"),
        ];
        let query = QueryHeader::new(contigs);

        let matches = engine.find_matches(&query, 5);
        assert!(!matches.is_empty());

        // Best match should be a GRCh38 reference
        let best = &matches[0];
        assert!(
            best.reference.display_name.contains("38")
                || best.reference.display_name.contains("hg38"),
            "Expected GRCh38 match, got: {}",
            best.reference.display_name
        );
    }

    #[test]
    fn test_find_matches_grch37() {
        let catalog = make_test_catalog();
        let engine = MatchingEngine::new(&catalog);

        // Create query with GRCh37 chromosomes
        let contigs = vec![
            Contig::new("chr1", 249_250_621).with_md5("1b22b98cdeb4a9304cb5d48026a85128"),
            Contig::new("chr2", 243_199_373).with_md5("a0d9851da00400dec1098a9255ac712e"),
        ];
        let query = QueryHeader::new(contigs);

        let matches = engine.find_matches(&query, 5);
        assert!(!matches.is_empty());

        // Best match should be a GRCh37/hg19 reference
        let best = &matches[0];
        assert!(
            best.reference.display_name.contains("37")
                || best.reference.display_name.contains("hg19")
                || best.reference.display_name.contains("b37"),
            "Expected GRCh37 match, got: {}",
            best.reference.display_name
        );
    }

    #[test]
    fn test_find_matches_ncbi_naming() {
        let catalog = make_test_catalog();
        let engine = MatchingEngine::new(&catalog);

        // Create query with NCBI naming but GRCh38 sequences
        let contigs = vec![
            Contig::new("1", 248_956_422).with_md5("6aef897c3d6ff0c78aff06ac189178dd"),
            Contig::new("2", 242_193_529).with_md5("f98db672eb0993dcfdabafe2a882905c"),
        ];
        let query = QueryHeader::new(contigs);

        let matches = engine.find_matches(&query, 5);
        assert!(!matches.is_empty());

        // Should still match GRCh38 references (score may be lower with only 2 contigs)
        let best = &matches[0];
        assert!(best.score.composite > 0.1); // Found a match above threshold
        assert!(
            best.reference.display_name.contains("38")
                || best.reference.display_name.contains("hs38"),
            "Expected GRCh38 match, got: {}",
            best.reference.display_name
        );
    }

    #[test]
    fn test_find_matches_no_md5() {
        let catalog = make_test_catalog();
        let engine = MatchingEngine::new(&catalog);

        // Create query without MD5s
        let contigs = vec![
            Contig::new("chr1", 248_956_422),
            Contig::new("chr2", 242_193_529),
            Contig::new("chr3", 198_295_559),
        ];
        let query = QueryHeader::new(contigs);

        // Should still find matches based on name+length
        let matches = engine.find_matches(&query, 5);
        assert!(!matches.is_empty());
    }

    #[test]
    fn test_find_best_match() {
        let catalog = make_test_catalog();
        let engine = MatchingEngine::new(&catalog);

        let contigs =
            vec![Contig::new("chr1", 248_956_422).with_md5("6aef897c3d6ff0c78aff06ac189178dd")];
        let query = QueryHeader::new(contigs);

        let best = engine.find_best_match(&query);
        assert!(best.is_some());
    }

    #[test]
    fn test_no_matches_for_invalid_contigs() {
        let catalog = make_test_catalog();
        let engine = MatchingEngine::new(&catalog);

        // Create query with invalid contigs
        let contigs = vec![
            Contig::new("invalid_contig", 12345),
            Contig::new("another_invalid", 67890),
        ];
        let query = QueryHeader::new(contigs);

        let matches = engine.find_matches(&query, 5);
        // Should return empty or low-scoring matches
        assert!(matches.is_empty() || matches[0].score.composite < 0.2);
    }

    // ========================================================================
    // UCSC-Style Patch Name Matching Integration Tests
    // ========================================================================

    #[test]
    fn test_ucsc_style_fix_patch_matching() {
        // Test that queries with UCSC-style fix-patch names match catalog references
        // The catalog should have these names as aliases (generated from NCBI reports)
        let catalog = make_test_catalog();
        let engine = MatchingEngine::new(&catalog);

        // Query with UCSC-style fix-patch name (chr1_KN196472v1_fix)
        // This should match grch38_v38 (or similar) which has KN196472.1 as a fix-patch
        let contigs = vec![
            // Primary chromosomes with MD5s for strong matching
            Contig::new("chr1", 248_956_422).with_md5("6aef897c3d6ff0c78aff06ac189178dd"),
            Contig::new("chr2", 242_193_529).with_md5("f98db672eb0993dcfdabafe2a882905c"),
            // UCSC-style fix-patch name (generated from KN196472.1)
            Contig::new("chr1_KN196472v1_fix", 186_494),
        ];
        let query = QueryHeader::new(contigs);

        let matches = engine.find_matches(&query, 5);
        assert!(
            !matches.is_empty(),
            "Should find matches for UCSC-style patch names"
        );

        // Best match should be a GRCh38 reference (since fix-patches are GRCh38-specific)
        let best = &matches[0];
        assert!(
            best.reference.display_name.contains("38")
                || best.reference.display_name.contains("hg38")
                || best.reference.display_name.contains("hs38"),
            "Expected GRCh38 match for fix-patch query, got: {}",
            best.reference.display_name
        );
    }

    #[test]
    fn test_ucsc_style_novel_patch_matching() {
        // Test that queries with UCSC-style novel-patch (alt) names match catalog references
        let catalog = make_test_catalog();
        let engine = MatchingEngine::new(&catalog);

        // Query with UCSC-style novel-patch name (chr1_KQ458382v1_alt)
        let contigs = vec![
            Contig::new("chr1", 248_956_422).with_md5("6aef897c3d6ff0c78aff06ac189178dd"),
            Contig::new("chr2", 242_193_529).with_md5("f98db672eb0993dcfdabafe2a882905c"),
            // UCSC-style novel-patch name (generated from KQ458382.1)
            Contig::new("chr1_KQ458382v1_alt", 141_019),
        ];
        let query = QueryHeader::new(contigs);

        let matches = engine.find_matches(&query, 5);
        assert!(
            !matches.is_empty(),
            "Should find matches for UCSC-style alt patch names"
        );

        // Best match should be a GRCh38 reference
        let best = &matches[0];
        assert!(
            best.reference.display_name.contains("38")
                || best.reference.display_name.contains("hg38")
                || best.reference.display_name.contains("hs38"),
            "Expected GRCh38 match for novel-patch query, got: {}",
            best.reference.display_name
        );
    }

    #[test]
    fn test_ucsc_style_y_chromosome_patch_matching() {
        // Test Y chromosome fix-patch matching
        let catalog = make_test_catalog();
        let engine = MatchingEngine::new(&catalog);

        let contigs = vec![
            Contig::new("chr1", 248_956_422).with_md5("6aef897c3d6ff0c78aff06ac189178dd"),
            Contig::new("chrY", 57_227_415).with_md5("ce3e31103314a704255f3cd90369ecce"),
            // Y chromosome fix-patch (generated from KN196487.1)
            Contig::new("chrY_KN196487v1_fix", 101_150),
        ];
        let query = QueryHeader::new(contigs);

        let matches = engine.find_matches(&query, 5);
        assert!(
            !matches.is_empty(),
            "Should find matches for Y chromosome fix-patch names"
        );
    }

    #[test]
    fn test_catalog_contains_ucsc_patch_aliases() {
        // Verify that the embedded catalog has UCSC-style names indexed for matching
        // Note: The matching tests above verify functional matching works. This test
        // verifies the catalog structure supports UCSC-style names.
        let catalog = make_test_catalog();

        // Check that some GRCh38 references have UCSC-style names indexed
        // Either in name_length_to_refs or alias_length_to_refs
        let has_fix_patch_index = catalog
            .name_length_to_refs
            .keys()
            .any(|(name, _)| name.contains("_fix") && name.starts_with("chr"))
            || catalog
                .alias_length_to_refs
                .keys()
                .any(|(name, _)| name.contains("_fix") && name.starts_with("chr"));

        let has_alt_patch_index = catalog.name_length_to_refs.keys().any(|(name, _)| {
            name.ends_with("_alt") && name.starts_with("chr") && name.contains('v')
        }) || catalog.alias_length_to_refs.keys().any(|(name, _)| {
            name.ends_with("_alt") && name.starts_with("chr") && name.contains('v')
        });

        // At minimum, the matching tests above prove UCSC patch names work.
        // This test is informational - if the catalog structure changes,
        // we want to ensure UCSC names are still indexed somewhere.
        if !has_fix_patch_index && !has_alt_patch_index {
            // Check if any reference has UCSC-style names in its name_length_set
            let any_ref_has_patches = catalog.references.iter().any(|r| {
                r.name_length_set.iter().any(|(name, _)| {
                    (name.contains("_fix") || name.contains("_alt"))
                        && name.starts_with("chr")
                        && name.contains('v')
                })
            });

            assert!(
                any_ref_has_patches,
                "Catalog should have UCSC-style patch names indexed for matching. \
                 The matching tests pass, so this may indicate a catalog structure change."
            );
        }
    }
}
