use std::collections::HashSet;

use crate::core::contig::Contig;
use crate::core::header::QueryHeader;
use crate::core::reference::KnownReference;
use crate::core::types::Confidence;

/// Classification of how a query contig matches a reference
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContigMatchType {
    /// Name+length match AND MD5 matches (or both are primary chromosomes with same length)
    Exact,
    /// Name+length match, but MD5 not available on one or both sides (neutral)
    NameLengthNoMd5,
    /// Name+length match, but MD5 differs (different sequence - heavy penalty!)
    Md5Conflict,
    /// No name+length match found
    Unmatched,
}

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
    /// Weighted composite score (the primary ranking metric)
    pub composite: f64,

    /// Confidence level derived from score
    pub confidence: Confidence,

    // === New per-contig match classification metrics ===
    /// Number of contigs with exact match (name+length+MD5 all match)
    pub exact_matches: usize,

    /// Number of contigs with name+length match but no MD5 available (neutral)
    pub name_length_matches: usize,

    /// Number of contigs with name+length match but MD5 conflicts (different sequence!)
    pub md5_conflicts: usize,

    /// Number of contigs with no match found
    pub unmatched: usize,

    // === Derived scores ===
    /// Per-contig match score: (exact + neutral + 0.1*conflicts) / total
    pub match_quality: f64,

    /// Coverage score: `good_matches` / `reference_contigs`
    pub coverage_score: f64,

    /// Fraction of contigs in correct relative order
    pub order_score: f64,

    /// Are matched contigs in the same order?
    pub order_preserved: bool,

    // === Legacy metrics (for backward compatibility) ===
    /// Jaccard similarity of MD5 sets: |intersection| / |union|
    pub md5_jaccard: f64,

    /// Jaccard similarity of (`normalized_name`, length) pairs
    pub name_length_jaccard: f64,

    /// Fraction of query contigs matched by MD5
    pub md5_query_coverage: f64,

    /// Fraction of query contigs matched by name+length
    pub name_length_query_coverage: f64,
}

impl MatchScore {
    /// Calculate match score between query and reference
    ///
    /// Uses the new contig-based scoring algorithm:
    /// - Classifies each query contig as Exact, `NameLengthNoMd5`, `Md5Conflict`, or Unmatched
    /// - Exact and `NameLengthNoMd5` get full credit (1.0)
    /// - `Md5Conflict` gets heavy penalty (0.1 credit - name is right but sequence is wrong)
    /// - Unmatched gets no credit
    ///
    /// Composite = 70% `match_quality` + 20% `coverage_score` + 10% `order_score`
    #[must_use]
    pub fn calculate(query: &QueryHeader, reference: &KnownReference) -> Self {
        // Classify each query contig
        let mut exact_matches = 0usize;
        let mut name_length_matches = 0usize; // No MD5 available
        let mut md5_conflicts = 0usize;
        let mut unmatched = 0usize;

        for contig in &query.contigs {
            match classify_contig_match(contig, reference) {
                ContigMatchType::Exact => exact_matches += 1,
                ContigMatchType::NameLengthNoMd5 => name_length_matches += 1,
                ContigMatchType::Md5Conflict => md5_conflicts += 1,
                ContigMatchType::Unmatched => unmatched += 1,
            }
        }

        let total_query = query.contigs.len();

        // Contig match score: full credit for exact/neutral, 10% for conflicts, 0 for unmatched
        // Key principle: MD5 absence is neutral (full credit), MD5 conflict is penalized
        let weighted_matches = count_to_f64(exact_matches)
            + count_to_f64(name_length_matches)
            + (count_to_f64(md5_conflicts) * 0.1);

        let match_quality = if total_query > 0 {
            weighted_matches / count_to_f64(total_query)
        } else {
            0.0
        };

        // Coverage: what fraction of reference is covered by good matches (exact + neutral)
        let good_matches = exact_matches + name_length_matches;
        let coverage_score = if reference.contigs.is_empty() {
            0.0
        } else {
            (count_to_f64(good_matches) / count_to_f64(reference.contigs.len())).min(1.0)
        };

        // Order analysis (existing function)
        let (order_preserved, order_score) = analyze_order(query, reference);

        // New composite formula: 70% match quality, 20% coverage, 10% order
        // Clamp to [0.0, 1.0] to handle floating point precision issues
        let composite = ((0.70 * match_quality) + (0.20 * coverage_score) + (0.10 * order_score))
            .clamp(0.0, 1.0);

        let confidence = Confidence::from_score(composite);

        // Compute legacy metrics for backward compatibility
        let md5_jaccard = jaccard_similarity(&query.md5_set, &reference.md5_set);
        let (name_length_jaccard, name_length_query_coverage) =
            calculate_name_length_similarity_with_aliases(query, reference);
        let md5_query_coverage = if query.md5_set.is_empty() {
            0.0
        } else {
            count_to_f64(query.md5_set.intersection(&reference.md5_set).count())
                / count_to_f64(query.md5_set.len())
        };

        Self {
            composite,
            confidence,
            exact_matches,
            name_length_matches,
            md5_conflicts,
            unmatched,
            match_quality,
            coverage_score,
            order_score,
            order_preserved,
            md5_jaccard,
            name_length_jaccard,
            md5_query_coverage,
            name_length_query_coverage,
        }
    }

    /// Calculate match score with custom scoring weights
    ///
    /// Uses the new contig-based scoring algorithm with configurable weights:
    /// - `contig_match_weight`: Weight for per-contig match score (default 70%)
    /// - `coverage_weight`: Weight for reference coverage (default 20%)
    /// - `order_weight`: Weight for contig ordering (default 10%)
    /// - `conflict_penalty`: Credit given to MD5 conflicts (default 0.1 = 10%)
    #[must_use]
    pub fn calculate_with_weights(
        query: &QueryHeader,
        reference: &KnownReference,
        weights: &crate::matching::engine::ScoringWeights,
    ) -> Self {
        // Classify each query contig
        let mut exact_matches = 0usize;
        let mut name_length_matches = 0usize;
        let mut md5_conflicts = 0usize;
        let mut unmatched = 0usize;

        for contig in &query.contigs {
            match classify_contig_match(contig, reference) {
                ContigMatchType::Exact => exact_matches += 1,
                ContigMatchType::NameLengthNoMd5 => name_length_matches += 1,
                ContigMatchType::Md5Conflict => md5_conflicts += 1,
                ContigMatchType::Unmatched => unmatched += 1,
            }
        }

        let total_query = query.contigs.len();

        // Normalize weights
        let normalized_weights = weights.normalized();

        // Contig match score using configurable conflict penalty
        let weighted_matches = count_to_f64(exact_matches)
            + count_to_f64(name_length_matches)
            + (count_to_f64(md5_conflicts) * normalized_weights.conflict_penalty);

        let match_quality = if total_query > 0 {
            weighted_matches / count_to_f64(total_query)
        } else {
            0.0
        };

        // Coverage: what fraction of reference is covered by good matches
        let good_matches = exact_matches + name_length_matches;
        let coverage_score = if reference.contigs.is_empty() {
            0.0
        } else {
            (count_to_f64(good_matches) / count_to_f64(reference.contigs.len())).min(1.0)
        };

        // Order analysis
        let (order_preserved, order_score) = analyze_order(query, reference);

        // Composite with custom weights
        // Clamp to [0.0, 1.0] to handle floating point precision issues
        let composite = ((normalized_weights.contig_match * match_quality)
            + (normalized_weights.coverage * coverage_score)
            + (normalized_weights.order * order_score))
            .clamp(0.0, 1.0);

        let confidence = Confidence::from_score(composite);

        // Compute legacy metrics for backward compatibility
        let md5_jaccard = jaccard_similarity(&query.md5_set, &reference.md5_set);
        let (name_length_jaccard, name_length_query_coverage) =
            calculate_name_length_similarity_with_aliases(query, reference);
        let md5_query_coverage = if query.md5_set.is_empty() {
            0.0
        } else {
            count_to_f64(query.md5_set.intersection(&reference.md5_set).count())
                / count_to_f64(query.md5_set.len())
        };

        Self {
            composite,
            confidence,
            exact_matches,
            name_length_matches,
            md5_conflicts,
            unmatched,
            match_quality,
            coverage_score,
            order_score,
            order_preserved,
            md5_jaccard,
            name_length_jaccard,
            md5_query_coverage,
            name_length_query_coverage,
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

    // Build position map for reference using exact names AND aliases
    // This allows order matching when query uses different naming convention
    let mut ref_positions: HashMap<(&str, u64), usize> = HashMap::new();
    for (i, contig) in reference.contigs.iter().enumerate() {
        // Add primary name
        ref_positions.insert((contig.name.as_str(), contig.length), i);
        // Add all aliases (so query "chr1" can find reference "1" with alias "chr1")
        for alias in &contig.aliases {
            ref_positions.insert((alias.as_str(), contig.length), i);
        }
    }

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

/// Classify how a query contig matches against a reference.
///
/// Returns the match type indicating whether the contig:
/// - Has an exact match (name+length+MD5 all match)
/// - Has a name+length match but no MD5 to compare (neutral)
/// - Has a name+length match but MD5 differs (conflict - different sequence!)
/// - Has no match
fn classify_contig_match(query_contig: &Contig, reference: &KnownReference) -> ContigMatchType {
    // First, find if there's a name+length match (including aliases)
    let matched_ref_contig = find_matching_reference_contig(query_contig, reference);

    match matched_ref_contig {
        None => ContigMatchType::Unmatched,
        Some(ref_contig) => {
            // Name+length matched, now check MD5
            match (&query_contig.md5, &ref_contig.md5) {
                // Both have MD5 - compare them
                (Some(q_md5), Some(r_md5)) => {
                    if q_md5.eq_ignore_ascii_case(r_md5) {
                        ContigMatchType::Exact
                    } else {
                        // MD5 mismatch - same name/length but different sequence!
                        ContigMatchType::Md5Conflict
                    }
                }
                // Either or both missing MD5 - neutral (can't verify but not a conflict)
                _ => ContigMatchType::NameLengthNoMd5,
            }
        }
    }
}

/// Find the reference contig that matches the query contig by name+length.
/// Checks both direct name match and alias matches.
fn find_matching_reference_contig<'a>(
    query_contig: &Contig,
    reference: &'a KnownReference,
) -> Option<&'a Contig> {
    // Check if query name+length matches any reference contig
    let query_key = (query_contig.name.clone(), query_contig.length);
    if reference.name_length_set.contains(&query_key) {
        // Find the actual contig (could be by name or alias)
        for ref_contig in &reference.contigs {
            if ref_contig.name == query_contig.name && ref_contig.length == query_contig.length {
                return Some(ref_contig);
            }
            // Check if query name matches a reference alias
            for alias in &ref_contig.aliases {
                if alias == &query_contig.name && ref_contig.length == query_contig.length {
                    return Some(ref_contig);
                }
            }
        }
    }

    // Check if any query alias matches a reference name or alias
    for query_alias in &query_contig.aliases {
        let alias_key = (query_alias.clone(), query_contig.length);
        if reference.name_length_set.contains(&alias_key) {
            // Find the actual contig
            for ref_contig in &reference.contigs {
                if ref_contig.name == *query_alias && ref_contig.length == query_contig.length {
                    return Some(ref_contig);
                }
                for ref_alias in &ref_contig.aliases {
                    if ref_alias == query_alias && ref_contig.length == query_contig.length {
                        return Some(ref_contig);
                    }
                }
            }
        }
    }

    None
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
                contig_match: 0.7,
                coverage: 0.2,
                order: 0.1,
                conflict_penalty: 0.1,
            },
            ScoringWeights {
                contig_match: 0.5,
                coverage: 0.3,
                order: 0.2,
                conflict_penalty: 0.0,
            },
            ScoringWeights {
                contig_match: 1.0,
                coverage: 1.0,
                order: 1.0,
                conflict_penalty: 0.5,
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

    // ========================================================================
    // New Scoring Algorithm Tests
    // ========================================================================

    #[test]
    fn test_perfect_name_length_match_no_md5_gets_high_score() {
        use crate::core::contig::Contig;
        use crate::core::reference::KnownReference;
        use crate::core::types::{Assembly, ReferenceSource};

        // Reference with 3 contigs, no MD5
        let ref_contigs = vec![
            Contig::new("chr1", 248_956_422),
            Contig::new("chr2", 242_193_529),
            Contig::new("chr3", 198_295_559),
        ];
        let reference = KnownReference::new(
            "test_ref",
            "Test Reference",
            Assembly::Grch38,
            ReferenceSource::Custom("test".to_string()),
        )
        .with_contigs(ref_contigs);

        // Query matches perfectly, also no MD5
        let query = QueryHeader::new(vec![
            Contig::new("chr1", 248_956_422),
            Contig::new("chr2", 242_193_529),
            Contig::new("chr3", 198_295_559),
        ]);

        let score = MatchScore::calculate(&query, &reference);

        // All contigs should be classified as NameLengthNoMd5 (neutral = full credit)
        assert_eq!(
            score.exact_matches, 0,
            "No exact matches (no MD5 available)"
        );
        assert_eq!(
            score.name_length_matches, 3,
            "All 3 contigs should be name+length matches"
        );
        assert_eq!(score.md5_conflicts, 0, "No conflicts");
        assert_eq!(score.unmatched, 0, "No unmatched");

        // Contig match score should be 1.0 (3 neutral matches / 3 total)
        assert!(
            (score.match_quality - 1.0).abs() < 0.001,
            "Contig match score should be 1.0, got {}",
            score.match_quality
        );

        // Coverage should be 1.0 (3 good matches / 3 reference contigs)
        assert!(
            (score.coverage_score - 1.0).abs() < 0.001,
            "Coverage score should be 1.0, got {}",
            score.coverage_score
        );

        // Composite should be high: 0.7*1.0 + 0.2*1.0 + 0.1*order = 0.9 + order
        // Order score should be 1.0 for 3 contigs in order
        assert!(
            score.composite > 0.9,
            "100% name+length match with no MD5 should score > 90%, got {:.1}%",
            score.composite * 100.0
        );
    }

    #[test]
    fn test_md5_conflict_gets_heavy_penalty() {
        use crate::core::contig::Contig;
        use crate::core::reference::KnownReference;
        use crate::core::types::{Assembly, ReferenceSource};

        // Reference with MD5
        let ref_contigs = vec![
            Contig::new("chr1", 248_956_422).with_md5("abc123"),
            Contig::new("chr2", 242_193_529).with_md5("def456"),
        ];
        let reference = KnownReference::new(
            "test_ref",
            "Test Reference",
            Assembly::Grch38,
            ReferenceSource::Custom("test".to_string()),
        )
        .with_contigs(ref_contigs);

        // Query with DIFFERENT MD5s (conflict!)
        let query = QueryHeader::new(vec![
            Contig::new("chr1", 248_956_422).with_md5("DIFFERENT1"),
            Contig::new("chr2", 242_193_529).with_md5("DIFFERENT2"),
        ]);

        let score = MatchScore::calculate(&query, &reference);

        // All contigs should be classified as Md5Conflict
        assert_eq!(score.exact_matches, 0, "No exact matches");
        assert_eq!(score.name_length_matches, 0, "No name+length-only matches");
        assert_eq!(score.md5_conflicts, 2, "Both contigs have MD5 conflicts");
        assert_eq!(score.unmatched, 0, "No unmatched");

        // Contig match score should be 0.1 * 2 / 2 = 0.1 (10% credit for conflicts)
        assert!(
            (score.match_quality - 0.1).abs() < 0.001,
            "Contig match score for conflicts should be 0.1, got {}",
            score.match_quality
        );

        // Composite should be low due to conflicts
        assert!(
            score.composite < 0.2,
            "MD5 conflicts should result in low score, got {:.1}%",
            score.composite * 100.0
        );
    }

    #[test]
    fn test_mixed_match_types_weighted_correctly() {
        use crate::core::contig::Contig;
        use crate::core::reference::KnownReference;
        use crate::core::types::{Assembly, ReferenceSource};

        // Reference: 3 contigs with MD5
        let ref_contigs = vec![
            Contig::new("chr1", 248_956_422).with_md5("correct_md5_1"),
            Contig::new("chr2", 242_193_529).with_md5("correct_md5_2"),
            Contig::new("chr3", 198_295_559), // No MD5
        ];
        let reference = KnownReference::new(
            "test_ref",
            "Test Reference",
            Assembly::Grch38,
            ReferenceSource::Custom("test".to_string()),
        )
        .with_contigs(ref_contigs);

        // Query: 4 contigs
        // - chr1: exact match (MD5 matches)
        // - chr2: conflict (MD5 differs)
        // - chr3: name+length only (no MD5 on either side)
        // - chr4: unmatched
        let query = QueryHeader::new(vec![
            Contig::new("chr1", 248_956_422).with_md5("correct_md5_1"),
            Contig::new("chr2", 242_193_529).with_md5("WRONG_MD5"),
            Contig::new("chr3", 198_295_559),
            Contig::new("chr4", 100_000), // Unmatched
        ]);

        let score = MatchScore::calculate(&query, &reference);

        assert_eq!(score.exact_matches, 1, "1 exact match (chr1)");
        assert_eq!(score.name_length_matches, 1, "1 name+length match (chr3)");
        assert_eq!(score.md5_conflicts, 1, "1 conflict (chr2)");
        assert_eq!(score.unmatched, 1, "1 unmatched (chr4)");

        // Contig match score = (1*1.0 + 1*1.0 + 1*0.1 + 0) / 4 = 2.1 / 4 = 0.525
        let expected_match_score = (1.0 + 1.0 + 0.1) / 4.0;
        assert!(
            (score.match_quality - expected_match_score).abs() < 0.001,
            "Expected match_quality = {}, got {}",
            expected_match_score,
            score.match_quality
        );
    }

    #[test]
    fn test_classify_contig_match_exact() {
        use crate::core::contig::Contig;
        use crate::core::reference::KnownReference;
        use crate::core::types::{Assembly, ReferenceSource};

        let ref_contigs = vec![Contig::new("chr1", 1000).with_md5("abc123")];
        let reference = KnownReference::new(
            "test",
            "Test",
            Assembly::Grch38,
            ReferenceSource::Custom("test".to_string()),
        )
        .with_contigs(ref_contigs);

        let query_contig = Contig::new("chr1", 1000).with_md5("abc123");
        assert_eq!(
            classify_contig_match(&query_contig, &reference),
            ContigMatchType::Exact
        );
    }

    #[test]
    fn test_classify_contig_match_no_md5() {
        use crate::core::contig::Contig;
        use crate::core::reference::KnownReference;
        use crate::core::types::{Assembly, ReferenceSource};

        let ref_contigs = vec![Contig::new("chr1", 1000)]; // No MD5
        let reference = KnownReference::new(
            "test",
            "Test",
            Assembly::Grch38,
            ReferenceSource::Custom("test".to_string()),
        )
        .with_contigs(ref_contigs);

        let query_contig = Contig::new("chr1", 1000); // No MD5
        assert_eq!(
            classify_contig_match(&query_contig, &reference),
            ContigMatchType::NameLengthNoMd5
        );
    }

    #[test]
    fn test_classify_contig_match_md5_conflict() {
        use crate::core::contig::Contig;
        use crate::core::reference::KnownReference;
        use crate::core::types::{Assembly, ReferenceSource};

        let ref_contigs = vec![Contig::new("chr1", 1000).with_md5("abc123")];
        let reference = KnownReference::new(
            "test",
            "Test",
            Assembly::Grch38,
            ReferenceSource::Custom("test".to_string()),
        )
        .with_contigs(ref_contigs);

        let query_contig = Contig::new("chr1", 1000).with_md5("DIFFERENT");
        assert_eq!(
            classify_contig_match(&query_contig, &reference),
            ContigMatchType::Md5Conflict
        );
    }

    #[test]
    fn test_classify_contig_match_unmatched() {
        use crate::core::contig::Contig;
        use crate::core::reference::KnownReference;
        use crate::core::types::{Assembly, ReferenceSource};

        let ref_contigs = vec![Contig::new("chr1", 1000)];
        let reference = KnownReference::new(
            "test",
            "Test",
            Assembly::Grch38,
            ReferenceSource::Custom("test".to_string()),
        )
        .with_contigs(ref_contigs);

        // Different name
        let query_contig = Contig::new("chr2", 1000);
        assert_eq!(
            classify_contig_match(&query_contig, &reference),
            ContigMatchType::Unmatched
        );

        // Different length
        let query_contig2 = Contig::new("chr1", 2000);
        assert_eq!(
            classify_contig_match(&query_contig2, &reference),
            ContigMatchType::Unmatched
        );
    }
}
