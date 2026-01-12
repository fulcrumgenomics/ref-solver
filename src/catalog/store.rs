use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;
use thiserror::Error;

use crate::core::reference::KnownReference;
use crate::core::types::ReferenceId;

#[derive(Error, Debug)]
pub enum CatalogError {
    #[error("Failed to read catalog: {0}")]
    ReadError(#[from] std::io::Error),

    #[error("Failed to parse catalog: {0}")]
    ParseError(#[from] serde_json::Error),
}

/// Catalog version for compatibility checking
pub const CATALOG_VERSION: &str = "1.0.0";

/// Serializable catalog format
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CatalogData {
    pub version: String,
    pub created_at: String,
    pub references: Vec<KnownReference>,
}

/// The main reference catalog with indexes
#[derive(Debug)]
pub struct ReferenceCatalog {
    /// All known references
    pub references: Vec<KnownReference>,

    /// Index: reference ID -> index in references vec
    id_to_index: HashMap<ReferenceId, usize>,

    /// Index: MD5 -> indices of references containing this MD5
    pub md5_to_refs: HashMap<String, Vec<usize>>,

    /// Index: (exact_name, length) -> indices of references
    pub name_length_to_refs: HashMap<(String, u64), Vec<usize>>,

    /// Index: (alias, length) -> indices of references
    /// Separate from name_length_to_refs to distinguish primary names from aliases
    pub alias_length_to_refs: HashMap<(String, u64), Vec<usize>>,

    /// Index: signature -> reference index (for exact matches)
    signature_to_ref: HashMap<String, usize>,
}

impl ReferenceCatalog {
    /// Create an empty catalog
    pub fn new() -> Self {
        Self {
            references: Vec::new(),
            id_to_index: HashMap::new(),
            md5_to_refs: HashMap::new(),
            name_length_to_refs: HashMap::new(),
            alias_length_to_refs: HashMap::new(),
            signature_to_ref: HashMap::new(),
        }
    }

    /// Load the embedded default catalog
    pub fn load_embedded() -> Result<Self, CatalogError> {
        // Embedded at compile time via build.rs
        const EMBEDDED_CATALOG: &str = include_str!("../../catalogs/human_references.json");
        Self::from_json(EMBEDDED_CATALOG)
    }

    /// Load catalog from a JSON file
    pub fn load_from_file(path: &Path) -> Result<Self, CatalogError> {
        let content = std::fs::read_to_string(path)?;
        Self::from_json(&content)
    }

    /// Parse catalog from JSON string
    pub fn from_json(json: &str) -> Result<Self, CatalogError> {
        let data: CatalogData = serde_json::from_str(json)?;

        // Version check (warn but don't fail)
        if data.version != CATALOG_VERSION {
            eprintln!(
                "Warning: Catalog version mismatch (expected {}, found {})",
                CATALOG_VERSION, data.version
            );
        }

        let mut catalog = Self::new();
        for mut reference in data.references {
            reference.rebuild_indexes();
            catalog.add_reference(reference);
        }

        Ok(catalog)
    }

    /// Add a reference to the catalog
    pub fn add_reference(&mut self, reference: KnownReference) {
        let index = self.references.len();

        // Index by ID
        self.id_to_index.insert(reference.id.clone(), index);

        // Index by MD5s
        for md5 in &reference.md5_set {
            self.md5_to_refs.entry(md5.clone()).or_default().push(index);
        }

        // Index by (name, length) pairs
        for (name, length) in &reference.name_length_set {
            self.name_length_to_refs
                .entry((name.clone(), *length))
                .or_default()
                .push(index);
        }

        // Index by (alias, length) pairs
        for contig in &reference.contigs {
            for alias in &contig.aliases {
                self.alias_length_to_refs
                    .entry((alias.clone(), contig.length))
                    .or_default()
                    .push(index);
            }
        }

        // Index by signature
        if let Some(sig) = &reference.signature {
            self.signature_to_ref.insert(sig.clone(), index);
        }

        self.references.push(reference);
    }

    /// Get a reference by ID
    pub fn get(&self, id: &ReferenceId) -> Option<&KnownReference> {
        self.id_to_index.get(id).map(|&idx| &self.references[idx])
    }

    /// Find exact match by signature
    pub fn find_by_signature(&self, signature: &str) -> Option<&KnownReference> {
        self.signature_to_ref
            .get(signature)
            .map(|&idx| &self.references[idx])
    }

    /// Export catalog to JSON
    pub fn to_json(&self) -> Result<String, CatalogError> {
        let data = CatalogData {
            version: CATALOG_VERSION.to_string(),
            created_at: chrono::Utc::now().to_rfc3339(),
            references: self.references.clone(),
        };
        Ok(serde_json::to_string_pretty(&data)?)
    }

    /// Number of references in catalog
    pub fn len(&self) -> usize {
        self.references.len()
    }

    /// Check if catalog is empty
    pub fn is_empty(&self) -> bool {
        self.references.is_empty()
    }
}

impl Default for ReferenceCatalog {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::contig::Contig;
    use crate::core::types::{Assembly, ReferenceSource};

    #[test]
    fn test_load_embedded_catalog() {
        let catalog = ReferenceCatalog::load_embedded().unwrap();
        assert!(!catalog.is_empty());
    }

    #[test]
    fn test_catalog_get_by_id() {
        let catalog = ReferenceCatalog::load_embedded().unwrap();

        let hg38 = catalog.get(&ReferenceId::new("hg38_ucsc"));
        assert!(hg38.is_some());
        let hg38 = hg38.unwrap();
        assert_eq!(hg38.display_name, "hg38 (UCSC)");
        assert!(!hg38.contigs.is_empty());
    }

    #[test]
    fn test_catalog_get_nonexistent() {
        let catalog = ReferenceCatalog::load_embedded().unwrap();
        let result = catalog.get(&ReferenceId::new("nonexistent_ref"));
        assert!(result.is_none());
    }

    #[test]
    fn test_catalog_to_json() {
        let catalog = ReferenceCatalog::load_embedded().unwrap();
        let json = catalog.to_json().unwrap();

        assert!(json.contains("\"version\""));
        assert!(json.contains("\"references\""));
        assert!(json.contains("hg38_ucsc"));
    }

    #[test]
    fn test_add_reference() {
        let mut catalog = ReferenceCatalog::new();
        assert_eq!(catalog.len(), 0);

        let contigs = vec![Contig::new("chr1", 100).with_md5("abc123")];
        let mut reference = KnownReference::new(
            "test_ref",
            "Test Reference",
            Assembly::Grch38,
            ReferenceSource::Custom("test".to_string()),
        );
        reference = reference.with_contigs(contigs);

        catalog.add_reference(reference);
        assert_eq!(catalog.len(), 1);

        let retrieved = catalog.get(&ReferenceId::new("test_ref"));
        assert!(retrieved.is_some());
        assert_eq!(retrieved.unwrap().display_name, "Test Reference");
    }
}
