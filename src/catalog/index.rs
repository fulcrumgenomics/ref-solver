use std::collections::{HashMap, HashSet};

use crate::core::header::QueryHeader;

use super::store::ReferenceCatalog;

/// Finds candidate references that might match a query header
pub struct CandidateFinder<'a> {
    catalog: &'a ReferenceCatalog,
}

impl<'a> CandidateFinder<'a> {
    pub fn new(catalog: &'a ReferenceCatalog) -> Self {
        Self { catalog }
    }

    /// Find candidate references based on MD5 overlap
    /// Returns indices sorted by number of overlapping MD5s (descending)
    pub fn find_candidates_by_md5(&self, query: &QueryHeader) -> Vec<(usize, usize)> {
        let mut ref_counts: HashMap<usize, usize> = HashMap::new();

        for md5 in &query.md5_set {
            if let Some(indices) = self.catalog.md5_to_refs.get(md5) {
                for &idx in indices {
                    *ref_counts.entry(idx).or_default() += 1;
                }
            }
        }

        let mut candidates: Vec<_> = ref_counts.into_iter().collect();
        candidates.sort_by(|a, b| b.1.cmp(&a.1)); // Sort by count descending
        candidates
    }

    /// Find candidates by (name, length) pairs when MD5s aren't available
    pub fn find_candidates_by_name_length(&self, query: &QueryHeader) -> Vec<(usize, usize)> {
        let mut ref_counts: HashMap<usize, usize> = HashMap::new();

        // Check primary names
        for (name, length) in &query.name_length_set {
            if let Some(indices) = self
                .catalog
                .name_length_to_refs
                .get(&(name.clone(), *length))
            {
                for &idx in indices {
                    *ref_counts.entry(idx).or_default() += 1;
                }
            }
        }

        // Also check query aliases against catalog names (reverse alias matching)
        // This handles cases where query has aliases like NC_000001.11 that match
        // catalog names directly
        for (alias, length) in &query.alias_length_set {
            if let Some(indices) = self
                .catalog
                .name_length_to_refs
                .get(&(alias.clone(), *length))
            {
                for &idx in indices {
                    *ref_counts.entry(idx).or_default() += 1;
                }
            }
        }

        let mut candidates: Vec<_> = ref_counts.into_iter().collect();
        candidates.sort_by(|a, b| b.1.cmp(&a.1));
        candidates
    }

    /// Get top N candidates combining MD5 and name/length matching
    pub fn find_top_candidates(&self, query: &QueryHeader, limit: usize) -> Vec<usize> {
        let mut seen: HashSet<usize> = HashSet::new();
        let mut result = Vec::new();

        // First add MD5-based candidates (higher priority)
        for (idx, _) in self.find_candidates_by_md5(query) {
            if seen.insert(idx) {
                result.push(idx);
                if result.len() >= limit {
                    return result;
                }
            }
        }

        // Then add name/length candidates
        for (idx, _) in self.find_candidates_by_name_length(query) {
            if seen.insert(idx) {
                result.push(idx);
                if result.len() >= limit {
                    return result;
                }
            }
        }

        result
    }
}
