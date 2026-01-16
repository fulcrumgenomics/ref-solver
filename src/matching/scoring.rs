use std::collections::HashSet;

use crate::core::header::QueryHeader;
use crate::core::reference::KnownReference;
use crate::core::types::Confidence;

/// Safely convert usize to f64 for percentage calculations
///
/// This function explicitly handles the precision loss that occurs when converting
/// usize to f64 on 64-bit platforms. For typical bioinformatics use cases, the
/// numbers are well within the safe range of f64 mantissa precision.
#[inline]
fn count_to_f64(count: usize) -> f64 {
    #[allow(clippy::cast_precision_loss)]
    {
        count as f64
    }
}

/// Detailed similarity scores between a query and a reference
#[derive(Debug, Clone)]
pub struct MatchScore {
    /// Jaccard similarity of MD5 sets: |intersection| / |union|
    pub md5_jaccard: f64,

    /// Jaccard similarity of (`normalized_name`, length) pairs
    pub name_length_jaccard: f64,

    /// Fraction of query contigs matched by MD5
    pub md5_query_coverage: f64,

    /// Fraction of query contigs matched by name+length
    pub name_length_query_coverage: f64,

    /// Are matched contigs in the same order?
    pub order_preserved: bool,

    /// Fraction of contigs in correct relative order
    pub order_score: f64,

    /// Weighted composite score
    pub composite: f64,

    /// Confidence level derived from score
    pub confidence: Confidence,
}

impl MatchScore {
    /// Calculate match score between query and reference
    #[must_use]
    pub fn calculate(query: &QueryHeader, reference: &KnownReference) -> Self {
        let md5_jaccard = jaccard_similarity(&query.md5_set, &reference.md5_set);

        // For name+length matching, also consider alias matches
        // A query contig matches if: query name matches ref name OR query alias matches ref name
        let (name_length_jaccard, name_length_query_coverage) =
            calculate_name_length_similarity_with_aliases(query, reference);

        // Query coverage
        let md5_query_coverage = if query.md5_set.is_empty() {
            0.0
        } else {
            count_to_f64(query.md5_set.intersection(&reference.md5_set).count())
                / count_to_f64(query.md5_set.len())
        };

        // Order analysis
        let (order_preserved, order_score) = analyze_order(query, reference);

        // Composite score with weights
        // MD5 is most reliable when available
        let has_md5 = !query.md5_set.is_empty() && !reference.md5_set.is_empty();

        let composite = if has_md5 {
            0.50 * md5_jaccard
                + 0.25 * name_length_jaccard
                + 0.15 * md5_query_coverage
                + 0.10 * order_score
        } else {
            0.60 * name_length_jaccard + 0.25 * name_length_query_coverage + 0.15 * order_score
        };

        let confidence = Confidence::from_score(composite);

        Self {
            md5_jaccard,
            name_length_jaccard,
            md5_query_coverage,
            name_length_query_coverage,
            order_preserved,
            order_score,
            composite,
            confidence,
        }
    }

    /// Calculate match score with custom scoring weights
    #[must_use]
    pub fn calculate_with_weights(
        query: &QueryHeader,
        reference: &KnownReference,
        weights: &crate::matching::engine::ScoringWeights,
    ) -> Self {
        let md5_jaccard = jaccard_similarity(&query.md5_set, &reference.md5_set);

        // For name+length matching, also consider alias matches
        let (name_length_jaccard, name_length_query_coverage) =
            calculate_name_length_similarity_with_aliases(query, reference);

        // Query coverage
        let md5_query_coverage = if query.md5_set.is_empty() {
            0.0
        } else {
            count_to_f64(query.md5_set.intersection(&reference.md5_set).count())
                / count_to_f64(query.md5_set.len())
        };

        // Order analysis
        let (order_preserved, order_score) = analyze_order(query, reference);

        // Normalize weights to ensure they sum to 1.0
        let normalized_weights = weights.normalized();

        // Composite score with custom weights
        let has_md5 = !query.md5_set.is_empty() && !reference.md5_set.is_empty();

        let composite = if has_md5 {
            // When MD5s are available, use all four components
            normalized_weights.md5_jaccard * md5_jaccard
                + normalized_weights.name_length_jaccard * name_length_jaccard
                + normalized_weights.md5_query_coverage * md5_query_coverage
                + normalized_weights.order_score * order_score
        } else {
            // Without MD5s, redistribute weights proportionally among non-MD5 components
            // The MD5 weights (md5_jaccard + md5_query_coverage) get distributed to the others
            let total_non_md5_weight =
                normalized_weights.name_length_jaccard + normalized_weights.order_score;

            if total_non_md5_weight > 0.0 {
                // Normalize the non-MD5 weights to sum to 1.0
                let redistributed_name_length =
                    normalized_weights.name_length_jaccard / total_non_md5_weight;
                let redistributed_order = normalized_weights.order_score / total_non_md5_weight;

                // Use name_length_query_coverage as the coverage metric when MD5 isn't available
                // Split name_length weight between jaccard and coverage
                let name_weight = redistributed_name_length * 0.7; // 70% to jaccard
                let coverage_weight = redistributed_name_length * 0.3; // 30% to coverage

                name_weight * name_length_jaccard
                    + coverage_weight * name_length_query_coverage
                    + redistributed_order * order_score
            } else {
                // Fallback to equal weights (sum to 1.0)
                0.6 * name_length_jaccard + 0.25 * name_length_query_coverage + 0.15 * order_score
            }
        };

        let confidence = Confidence::from_score(composite);

        Self {
            md5_jaccard,
            name_length_jaccard,
            md5_query_coverage,
            name_length_query_coverage,
            order_preserved,
            order_score,
            composite,
            confidence,
        }
    }
}

/// Jaccard similarity: |A ∩ B| / |A ∪ B|
///
/// Returns 0.0 when both sets are empty (undefined mathematically, but 0.0 is safer
/// for matching to avoid false positives from two references with no MD5s).
fn jaccard_similarity<T: Eq + std::hash::Hash>(a: &HashSet<T>, b: &HashSet<T>) -> f64 {
    let intersection = a.intersection(b).count();
    let union = a.union(b).count();
    if union == 0 {
        // Both sets empty - return 0.0 to avoid false positive matches
        0.0
    } else {
        count_to_f64(intersection) / count_to_f64(union)
    }
}

/// Calculate name+length similarity accounting for alias-based matching.
///
/// A query contig matches a reference contig if:
/// - Query name+length matches reference name+length directly, OR
/// - Query alias+length matches reference name+length (reverse alias matching)
///
/// Returns (`jaccard_similarity`, `query_coverage`)
fn calculate_name_length_similarity_with_aliases(
    query: &QueryHeader,
    reference: &KnownReference,
) -> (f64, f64) {
    // Count query contigs that match (by name or alias)
    let mut matched_query_contigs = 0usize;
    let mut matched_ref_keys: HashSet<(String, u64)> = HashSet::new();

    for contig in &query.contigs {
        let name_key = (contig.name.clone(), contig.length);

        // Check direct name match
        if reference.name_length_set.contains(&name_key) {
            matched_query_contigs += 1;
            matched_ref_keys.insert(name_key);
            continue;
        }

        // Check alias matches (query alias -> ref name)
        for alias in &contig.aliases {
            let alias_key = (alias.clone(), contig.length);
            if reference.name_length_set.contains(&alias_key) {
                matched_query_contigs += 1;
                matched_ref_keys.insert(alias_key);
                break;
            }
        }
    }

    // Calculate Jaccard: |matched| / |union|
    // Union = query contigs + ref contigs - matched (to avoid double counting)
    let query_count = query.contigs.len();
    let ref_count = reference.contigs.len();
    let union_size = query_count + ref_count - matched_ref_keys.len();

    let jaccard = if union_size == 0 {
        0.0
    } else {
        count_to_f64(matched_ref_keys.len()) / count_to_f64(union_size)
    };

    let coverage = if query_count == 0 {
        0.0
    } else {
        count_to_f64(matched_query_contigs) / count_to_f64(query_count)
    };

    (jaccard, coverage)
}

/// Analyze if contigs are in the same order.
/// Considers both direct name matches and alias-based matches.
///
/// Returns (`order_preserved`, `order_score`) where:
/// - `order_preserved`: true if all matched contigs appear in the same order
/// - `order_score`: fraction of contigs in the longest increasing subsequence
///
/// When there are < 2 matching contigs, returns (true, 0.0) since order
/// cannot be meaningfully assessed with so few matches.
fn analyze_order(query: &QueryHeader, reference: &KnownReference) -> (bool, f64) {
    use std::collections::HashMap;

    // Build position map for reference using exact names
    let ref_positions: HashMap<(&str, u64), usize> = reference
        .contigs
        .iter()
        .enumerate()
        .map(|(i, c)| ((c.name.as_str(), c.length), i))
        .collect();

    // Get positions of query contigs in reference ordering
    // Check both direct name match and alias matches
    let mut positions: Vec<usize> = Vec::new();
    for contig in &query.contigs {
        // Try direct name match first
        let key = (contig.name.as_str(), contig.length);
        if let Some(&pos) = ref_positions.get(&key) {
            positions.push(pos);
            continue;
        }

        // Try alias matches
        for alias in &contig.aliases {
            let alias_key = (alias.as_str(), contig.length);
            if let Some(&pos) = ref_positions.get(&alias_key) {
                positions.push(pos);
                break;
            }
        }
    }

    // Need at least 2 matching contigs to assess ordering
    if positions.len() < 2 {
        // Return neutral values: order is trivially preserved, but score is 0
        // to avoid inflating composite scores with insufficient data
        return (true, 0.0);
    }

    // Check if sorted
    let is_sorted = positions.windows(2).all(|w| w[0] < w[1]);

    // Calculate longest increasing subsequence for order score
    let lis_len = longest_increasing_subsequence(&positions);
    let order_score = count_to_f64(lis_len) / count_to_f64(positions.len());

    (is_sorted, order_score)
}

/// Length of longest increasing subsequence (LIS).
///
/// Returns 0 for empty arrays, otherwise the length of the LIS (minimum 1).
fn longest_increasing_subsequence(arr: &[usize]) -> usize {
    if arr.is_empty() {
        return 0;
    }

    // dp[i] = length of LIS ending at index i (minimum 1 for any single element)
    let mut dp = vec![1usize; arr.len()];

    for i in 1..arr.len() {
        for j in 0..i {
            if arr[j] < arr[i] {
                dp[i] = dp[i].max(dp[j] + 1);
            }
        }
    }

    // Safety: dp is non-empty (arr.len() >= 1 after early return check)
    dp.into_iter().max().expect("dp is non-empty")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_jaccard_similarity() {
        let a: HashSet<i32> = [1, 2, 3].into_iter().collect();
        let b: HashSet<i32> = [2, 3, 4].into_iter().collect();

        let similarity = jaccard_similarity(&a, &b);
        // intersection = {2, 3} = 2, union = {1, 2, 3, 4} = 4
        assert!((similarity - 0.5).abs() < 0.001);

        // Empty sets should return 0.0 to avoid false positives
        let empty1: HashSet<i32> = HashSet::new();
        let empty2: HashSet<i32> = HashSet::new();
        assert!((jaccard_similarity(&empty1, &empty2) - 0.0).abs() < 0.001);

        // One empty, one non-empty should return 0.0
        assert!((jaccard_similarity(&a, &empty1) - 0.0).abs() < 0.001);

        // Identical sets should return 1.0
        let c: HashSet<i32> = [1, 2, 3].into_iter().collect();
        assert!((jaccard_similarity(&a, &c) - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_longest_increasing_subsequence() {
        assert_eq!(longest_increasing_subsequence(&[1, 2, 3, 4, 5]), 5);
        assert_eq!(longest_increasing_subsequence(&[5, 4, 3, 2, 1]), 1);
        assert_eq!(longest_increasing_subsequence(&[1, 3, 2, 4, 5]), 4);
        assert_eq!(longest_increasing_subsequence(&[]), 0);
    }

    #[test]
    fn test_composite_score_never_exceeds_one() {
        use crate::core::contig::Contig;
        use crate::core::reference::KnownReference;
        use crate::core::types::{Assembly, ReferenceSource};

        // Create a reference with some contigs
        let ref_contigs = vec![
            Contig::new("chr1", 1000),
            Contig::new("chr2", 2000),
            Contig::new("chr3", 3000),
        ];
        let reference = KnownReference::new(
            "test_ref",
            "Test Reference",
            Assembly::Grch38,
            ReferenceSource::Custom("test".to_string()),
        )
        .with_contigs(ref_contigs);

        // Test with query that has NO MD5s (triggers the redistribution code path)
        let query_no_md5 = QueryHeader::new(vec![
            Contig::new("chr1", 1000),
            Contig::new("chr2", 2000),
            Contig::new("chr3", 3000),
        ]);

        let score_no_md5 = MatchScore::calculate(&query_no_md5, &reference);
        assert!(
            score_no_md5.composite <= 1.0,
            "Composite score without MD5 should not exceed 1.0, got {}",
            score_no_md5.composite
        );
        assert!(
            score_no_md5.composite >= 0.0,
            "Composite score should not be negative, got {}",
            score_no_md5.composite
        );

        // Test with query that HAS MD5s
        let query_with_md5 = QueryHeader::new(vec![
            Contig::new("chr1", 1000).with_md5("abc123"),
            Contig::new("chr2", 2000).with_md5("def456"),
            Contig::new("chr3", 3000).with_md5("ghi789"),
        ]);

        let score_with_md5 = MatchScore::calculate(&query_with_md5, &reference);
        assert!(
            score_with_md5.composite <= 1.0,
            "Composite score with MD5 should not exceed 1.0, got {}",
            score_with_md5.composite
        );
        assert!(
            score_with_md5.composite >= 0.0,
            "Composite score should not be negative, got {}",
            score_with_md5.composite
        );
    }

    #[test]
    fn test_composite_score_with_custom_weights_never_exceeds_one() {
        use crate::core::contig::Contig;
        use crate::core::reference::KnownReference;
        use crate::core::types::{Assembly, ReferenceSource};
        use crate::matching::engine::ScoringWeights;

        let ref_contigs = vec![Contig::new("chr1", 1000), Contig::new("chr2", 2000)];
        let reference = KnownReference::new(
            "test_ref",
            "Test Reference",
            Assembly::Grch38,
            ReferenceSource::Custom("test".to_string()),
        )
        .with_contigs(ref_contigs);

        let query = QueryHeader::new(vec![Contig::new("chr1", 1000), Contig::new("chr2", 2000)]);

        // Test with various weight configurations
        let weight_configs = vec![
            ScoringWeights {
                md5_jaccard: 0.4,
                name_length_jaccard: 0.3,
                md5_query_coverage: 0.2,
                order_score: 0.1,
            },
            ScoringWeights {
                md5_jaccard: 0.0,
                name_length_jaccard: 0.5,
                md5_query_coverage: 0.0,
                order_score: 0.5,
            },
            ScoringWeights {
                md5_jaccard: 1.0,
                name_length_jaccard: 1.0,
                md5_query_coverage: 1.0,
                order_score: 1.0,
            },
        ];

        for weights in weight_configs {
            let score = MatchScore::calculate_with_weights(&query, &reference, &weights);
            assert!(
                score.composite <= 1.0,
                "Composite score should not exceed 1.0 with weights {:?}, got {}",
                weights,
                score.composite
            );
            assert!(
                score.composite >= 0.0,
                "Composite score should not be negative with weights {:?}, got {}",
                weights,
                score.composite
            );
        }
    }

    // Tests for all 4 alias matching combinations:
    // 1. Query name → Ref name (direct match)
    // 2. Query name → Ref alias (UCSC query chr1 matches NCBI ref with chr1 alias)
    // 3. Query alias → Ref name (NCBI query with chr1 alias matches UCSC ref chr1)
    // 4. Query alias → Ref alias (both have aliases that match)

    #[test]
    fn test_alias_match_query_name_to_ref_name() {
        use crate::core::contig::Contig;
        use crate::core::reference::KnownReference;
        use crate::core::types::{Assembly, ReferenceSource};

        // Direct match: query chr1 → ref chr1
        let ref_contigs = vec![
            Contig::new("chr1", 248_956_422),
            Contig::new("chr2", 242_193_529),
        ];
        let reference = KnownReference::new(
            "test_ref",
            "Test Reference",
            Assembly::Grch38,
            ReferenceSource::Custom("test".to_string()),
        )
        .with_contigs(ref_contigs);

        let query = QueryHeader::new(vec![
            Contig::new("chr1", 248_956_422),
            Contig::new("chr2", 242_193_529),
        ]);

        let score = MatchScore::calculate(&query, &reference);
        assert!(
            (score.name_length_jaccard - 1.0).abs() < 0.001,
            "Direct name match should have jaccard = 1.0, got {}",
            score.name_length_jaccard
        );
    }

    #[test]
    fn test_alias_match_query_name_to_ref_alias() {
        use crate::core::contig::Contig;
        use crate::core::reference::KnownReference;
        use crate::core::types::{Assembly, ReferenceSource};

        // Query uses UCSC names (chr1), ref uses NCBI names with chr1 as alias
        let ref_contigs = vec![
            Contig::new("NC_000001.11", 248_956_422)
                .with_aliases(vec!["chr1".to_string(), "1".to_string()]),
            Contig::new("NC_000002.12", 242_193_529)
                .with_aliases(vec!["chr2".to_string(), "2".to_string()]),
        ];
        let reference = KnownReference::new(
            "test_ncbi_ref",
            "Test NCBI Reference",
            Assembly::Grch38,
            ReferenceSource::Custom("test".to_string()),
        )
        .with_contigs(ref_contigs);

        // Query with UCSC names - NO aliases
        let query = QueryHeader::new(vec![
            Contig::new("chr1", 248_956_422),
            Contig::new("chr2", 242_193_529),
        ]);

        let score = MatchScore::calculate(&query, &reference);
        assert!(
            (score.name_length_jaccard - 1.0).abs() < 0.001,
            "Query name matching ref alias should have jaccard = 1.0, got {}",
            score.name_length_jaccard
        );
    }

    #[test]
    fn test_alias_match_query_alias_to_ref_name() {
        use crate::core::contig::Contig;
        use crate::core::reference::KnownReference;
        use crate::core::types::{Assembly, ReferenceSource};

        // Ref uses UCSC names (chr1), query uses NCBI names with chr1 as alias
        let ref_contigs = vec![
            Contig::new("chr1", 248_956_422),
            Contig::new("chr2", 242_193_529),
        ];
        let reference = KnownReference::new(
            "test_ucsc_ref",
            "Test UCSC Reference",
            Assembly::Grch38,
            ReferenceSource::Custom("test".to_string()),
        )
        .with_contigs(ref_contigs);

        // Query with NCBI names that have UCSC aliases
        let query = QueryHeader::new(vec![
            Contig::new("NC_000001.11", 248_956_422)
                .with_aliases(vec!["chr1".to_string(), "1".to_string()]),
            Contig::new("NC_000002.12", 242_193_529)
                .with_aliases(vec!["chr2".to_string(), "2".to_string()]),
        ]);

        let score = MatchScore::calculate(&query, &reference);
        assert!(
            (score.name_length_jaccard - 1.0).abs() < 0.001,
            "Query alias matching ref name should have jaccard = 1.0, got {}",
            score.name_length_jaccard
        );
    }

    #[test]
    fn test_alias_match_query_alias_to_ref_alias() {
        use crate::core::contig::Contig;
        use crate::core::reference::KnownReference;
        use crate::core::types::{Assembly, ReferenceSource};

        // Both query and ref use NCBI names with UCSC aliases
        // They should match via their common alias
        let ref_contigs = vec![
            Contig::new("NC_000001.11", 248_956_422)
                .with_aliases(vec!["chr1".to_string(), "CM000663.2".to_string()]),
            Contig::new("NC_000002.12", 242_193_529)
                .with_aliases(vec!["chr2".to_string(), "CM000664.2".to_string()]),
        ];
        let reference = KnownReference::new(
            "test_ncbi_ref_1",
            "Test NCBI Reference 1",
            Assembly::Grch38,
            ReferenceSource::Custom("test".to_string()),
        )
        .with_contigs(ref_contigs);

        // Query with different primary names but overlapping aliases
        let query = QueryHeader::new(vec![
            Contig::new("CM000663.2", 248_956_422)
                .with_aliases(vec!["chr1".to_string(), "NC_000001.11".to_string()]),
            Contig::new("CM000664.2", 242_193_529)
                .with_aliases(vec!["chr2".to_string(), "NC_000002.12".to_string()]),
        ]);

        let score = MatchScore::calculate(&query, &reference);
        assert!(
            (score.name_length_jaccard - 1.0).abs() < 0.001,
            "Query alias matching ref alias should have jaccard = 1.0, got {}",
            score.name_length_jaccard
        );
    }
}
