use std::collections::{HashMap, HashSet};

use crate::core::header::QueryHeader;

use super::store::ReferenceCatalog;

/// Finds candidate references that might match a query header
pub struct CandidateFinder<'a> {
    catalog: &'a ReferenceCatalog,
}

impl<'a> CandidateFinder<'a> {
    #[must_use]
    pub fn new(catalog: &'a ReferenceCatalog) -> Self {
        Self { catalog }
    }

    /// Find candidate references based on MD5 overlap
    /// Returns indices sorted by number of overlapping MD5s (descending)
    #[must_use]
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
    #[must_use]
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
    #[must_use]
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::contig::Contig;
    use crate::core::reference::KnownReference;
    use crate::core::types::{Assembly, ReferenceSource};

    #[test]
    fn test_find_candidates_via_alias() {
        // Create a catalog with an NCBI reference that has UCSC aliases
        let mut catalog = ReferenceCatalog::new();

        let ref_contigs = vec![
            Contig::new("NC_000001.11", 248_956_422)
                .with_md5("6aef897c3d6ff0c78aff06ac189178dd")
                .with_aliases(vec!["chr1".to_string(), "1".to_string()]),
            Contig::new("NC_000002.12", 242_193_529)
                .with_md5("f98db672eb0993dcfdabafe2a882905c")
                .with_aliases(vec!["chr2".to_string(), "2".to_string()]),
        ];

        let reference = KnownReference::new(
            "test_ncbi_ref",
            "Test NCBI Reference",
            Assembly::Grch38,
            ReferenceSource::Custom("test".to_string()),
        )
        .with_contigs(ref_contigs);

        catalog.add_reference(reference);

        // Create a query with UCSC names (no aliases, no MD5s)
        let query = QueryHeader::new(vec![
            Contig::new("chr1", 248_956_422),
            Contig::new("chr2", 242_193_529),
        ]);

        // The CandidateFinder should find the NCBI reference via alias matching
        let finder = CandidateFinder::new(&catalog);
        let candidates = finder.find_candidates_by_name_length(&query);

        assert!(
            !candidates.is_empty(),
            "Should find NCBI reference as candidate via UCSC aliases"
        );

        // The reference should have matched 2 contigs via alias
        let (idx, count) = candidates[0];
        assert_eq!(idx, 0, "Should match the first (and only) reference");
        assert_eq!(count, 2, "Should match both chr1 and chr2 via aliases");
    }

    #[test]
    fn test_find_candidates_real_catalog() {
        // Test with the real embedded catalog
        let catalog = ReferenceCatalog::load_embedded().unwrap();

        // Create a query with UCSC names (matching hg38)
        let query = QueryHeader::new(vec![
            Contig::new("chr1", 248_956_422),
            Contig::new("chr2", 242_193_529),
        ]);

        let finder = CandidateFinder::new(&catalog);
        let candidates = finder.find_candidates_by_name_length(&query);

        // Should find some candidates
        assert!(
            !candidates.is_empty(),
            "Should find candidates for UCSC chr1/chr2 query"
        );

        // Check if grch38_v38 (p12) is among the candidates
        let grch38_v38_idx = catalog
            .references
            .iter()
            .position(|r| r.id.0 == "grch38_v38");

        if let Some(expected_idx) = grch38_v38_idx {
            let found = candidates.iter().any(|(idx, _)| *idx == expected_idx);
            assert!(
                found,
                "grch38_v38 (with UCSC aliases) should be found as candidate for UCSC query"
            );
        }
    }
}
