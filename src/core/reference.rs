use serde::{Deserialize, Serialize};
use std::collections::HashSet;

use crate::core::contig::{detect_naming_convention, Contig, SequenceRole};
use crate::core::types::{Assembly, NamingConvention, ReferenceId, ReferenceSource};
use crate::utils::validation::compute_signature as compute_sig;

/// A known reference genome in the catalog
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KnownReference {
    /// Unique identifier
    pub id: ReferenceId,

    /// Human-readable display name
    pub display_name: String,

    /// Assembly version
    pub assembly: Assembly,

    /// Source organization
    pub source: ReferenceSource,

    /// Naming convention used
    pub naming_convention: NamingConvention,

    /// Download URL for the reference FASTA
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub download_url: Option<String>,

    /// Path to NCBI assembly report (if applicable)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub assembly_report_url: Option<String>,

    /// All contigs in this reference
    pub contigs: Vec<Contig>,

    /// Description/notes about this reference
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,

    /// Tags for filtering (e.g., "`with_decoy`", "`no_alt`", "`analysis_set`")
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub tags: Vec<String>,

    /// Contigs that appear in the assembly report but not in the FASTA/dict
    /// (e.g., MT in CHM13 which uses standard rCRS mitochondria)
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub contigs_missing_from_fasta: Vec<String>,

    // === Pre-computed for fast matching (populated on load) ===
    /// Set of all MD5 checksums in this reference
    #[serde(skip)]
    pub md5_set: HashSet<String>,

    /// Set of all (`exact_name`, length) pairs for matching
    #[serde(skip)]
    pub name_length_set: HashSet<(String, u64)>,

    /// Signature for exact matching (hash of sorted MD5s)
    #[serde(skip)]
    pub signature: Option<String>,
}

impl KnownReference {
    pub fn new(
        id: impl Into<String>,
        display_name: impl Into<String>,
        assembly: Assembly,
        source: ReferenceSource,
    ) -> Self {
        Self {
            id: ReferenceId::new(id),
            display_name: display_name.into(),
            assembly,
            source,
            naming_convention: NamingConvention::Mixed,
            download_url: None,
            assembly_report_url: None,
            contigs: Vec::new(),
            description: None,
            tags: Vec::new(),
            contigs_missing_from_fasta: Vec::new(),
            md5_set: HashSet::new(),
            name_length_set: HashSet::new(),
            signature: None,
        }
    }

    #[must_use]
    pub fn with_contigs(mut self, contigs: Vec<Contig>) -> Self {
        self.naming_convention = detect_naming_convention(&contigs);
        self.contigs = contigs;
        self.rebuild_indexes();
        self
    }

    /// Rebuild the internal indexes after modifying contigs
    pub fn rebuild_indexes(&mut self) {
        self.md5_set.clear();
        self.name_length_set.clear();

        for contig in &self.contigs {
            if let Some(md5) = &contig.md5 {
                self.md5_set.insert(md5.clone());
            }
            // Use exact name for matching (no normalization)
            self.name_length_set
                .insert((contig.name.clone(), contig.length));
        }

        // Compute signature from sorted MD5s
        self.signature = self.compute_signature();
    }

    /// Compute a signature for exact matching
    /// Uses sorted MD5s concatenated and hashed
    fn compute_signature(&self) -> Option<String> {
        let sig = compute_sig(&self.md5_set);
        if sig.is_empty() {
            None
        } else {
            Some(sig)
        }
    }

    /// Check if this reference has decoy sequences
    #[must_use]
    pub fn has_decoy(&self) -> bool {
        self.contigs.iter().any(super::contig::Contig::is_decoy)
    }

    /// Check if this reference has ALT contigs
    #[must_use]
    pub fn has_alt(&self) -> bool {
        self.contigs.iter().any(super::contig::Contig::is_alt)
    }

    /// Count contigs by sequence role
    #[must_use]
    pub fn role_counts(&self) -> RoleCounts {
        let mut counts = RoleCounts::default();
        for contig in &self.contigs {
            match contig.sequence_role {
                SequenceRole::AssembledMolecule => counts.assembled_molecule += 1,
                SequenceRole::AltScaffold => counts.alt_scaffold += 1,
                SequenceRole::FixPatch => counts.fix_patch += 1,
                SequenceRole::NovelPatch => counts.novel_patch += 1,
                SequenceRole::UnlocalizedScaffold => counts.unlocalized_scaffold += 1,
                SequenceRole::UnplacedScaffold => counts.unplaced_scaffold += 1,
                SequenceRole::Unknown => counts.unknown += 1,
            }
        }
        counts
    }
}

/// Counts of contigs by sequence role
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct RoleCounts {
    pub assembled_molecule: usize,
    pub alt_scaffold: usize,
    pub fix_patch: usize,
    pub novel_patch: usize,
    pub unlocalized_scaffold: usize,
    pub unplaced_scaffold: usize,
    pub unknown: usize,
}
