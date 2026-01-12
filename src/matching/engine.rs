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
    pub fn new(reference: &KnownReference, query: &QueryHeader) -> Self {
        let score = MatchScore::calculate(query, reference);
        let diagnosis = MatchDiagnosis::analyze(query, reference);

        Self {
            reference: reference.clone(),
            score,
            diagnosis,
        }
    }

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

/// Configurable weights for different scoring components
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ScoringWeights {
    /// Weight for MD5 Jaccard similarity
    pub md5_jaccard: f64,
    /// Weight for name+length Jaccard similarity
    pub name_length_jaccard: f64,
    /// Weight for MD5 query coverage
    pub md5_query_coverage: f64,
    /// Weight for order score
    pub order_score: f64,
}

impl Default for ScoringWeights {
    fn default() -> Self {
        Self {
            md5_jaccard: 0.4,         // 40%
            name_length_jaccard: 0.3, // 30%
            md5_query_coverage: 0.2,  // 20%
            order_score: 0.1,         // 10%
        }
    }
}

impl ScoringWeights {
    /// Normalize weights to sum to 1.0
    pub fn normalized(&self) -> Self {
        let total = self.md5_jaccard
            + self.name_length_jaccard
            + self.md5_query_coverage
            + self.order_score;

        if total <= 0.0 {
            return Self::default();
        }

        Self {
            md5_jaccard: self.md5_jaccard / total,
            name_length_jaccard: self.name_length_jaccard / total,
            md5_query_coverage: self.md5_query_coverage / total,
            order_score: self.order_score / total,
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
    pub fn new(catalog: &'a ReferenceCatalog) -> Self {
        Self {
            catalog,
            config: MatchingConfig::default(),
        }
    }

    /// Create a new matching engine with custom configuration
    pub fn with_config(catalog: &'a ReferenceCatalog, config: MatchingConfig) -> Self {
        Self { catalog, config }
    }

    /// Find the best matching references for a query
    pub fn find_matches(&self, query: &QueryHeader, limit: usize) -> Vec<MatchResult> {
        // Step 1: Try exact signature match
        if let Some(sig) = &query.signature {
            if let Some(reference) = self.catalog.find_by_signature(sig) {
                return vec![MatchResult::new(reference, query)];
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
}
