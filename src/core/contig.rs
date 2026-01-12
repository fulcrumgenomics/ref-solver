use serde::{Deserialize, Serialize};

use crate::core::types::NamingConvention;

/// Sequence role from NCBI assembly report
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "kebab-case")]
pub enum SequenceRole {
    /// Primary chromosome (1-22, X, Y, MT)
    AssembledMolecule,
    /// Alternate locus scaffold
    AltScaffold,
    /// Fix patch (error correction)
    FixPatch,
    /// Novel patch (new sequence)
    NovelPatch,
    /// Unlocalized scaffold (known chromosome, unknown location)
    UnlocalizedScaffold,
    /// Unplaced scaffold (unknown chromosome)
    UnplacedScaffold,
    /// Role not specified or unknown
    #[default]
    Unknown,
}

impl SequenceRole {
    /// Parse a sequence role from string representation (e.g. from NCBI assembly report)
    pub fn parse(s: &str) -> Self {
        match s.to_lowercase().as_str() {
            "assembled-molecule" => SequenceRole::AssembledMolecule,
            "alt-scaffold" => SequenceRole::AltScaffold,
            "fix-patch" => SequenceRole::FixPatch,
            "novel-patch" => SequenceRole::NovelPatch,
            "unlocalized-scaffold" => SequenceRole::UnlocalizedScaffold,
            "unplaced-scaffold" => SequenceRole::UnplacedScaffold,
            _ => SequenceRole::Unknown,
        }
    }
}

/// A single contig/sequence in a reference genome
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Contig {
    /// Sequence name (SN tag in SAM)
    pub name: String,

    /// Sequence length (LN tag in SAM)
    pub length: u64,

    /// MD5 checksum of the sequence (M5 tag in SAM)
    /// Lowercase hex, 32 characters
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub md5: Option<String>,

    /// Assembly identifier (AS tag in SAM)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub assembly: Option<String>,

    /// URI where sequence can be retrieved (UR tag in SAM)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub uri: Option<String>,

    /// Species (SP tag in SAM)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub species: Option<String>,

    /// Known alternative names for this contig
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub aliases: Vec<String>,

    /// Sequence role from NCBI assembly report
    #[serde(default, skip_serializing_if = "is_unknown_role")]
    pub sequence_role: SequenceRole,
}

fn is_unknown_role(role: &SequenceRole) -> bool {
    matches!(role, SequenceRole::Unknown)
}

impl Contig {
    pub fn new(name: impl Into<String>, length: u64) -> Self {
        Self {
            name: name.into(),
            length,
            md5: None,
            assembly: None,
            uri: None,
            species: None,
            aliases: Vec::new(),
            sequence_role: SequenceRole::Unknown,
        }
    }

    #[cfg(test)]
    pub fn with_md5(mut self, md5: impl Into<String>) -> Self {
        self.md5 = Some(md5.into());
        self
    }

    /// Check if this contig is a primary chromosome (1-22, X, Y)
    /// Matches both UCSC (chr1) and NCBI (1) naming conventions exactly
    pub fn is_primary_chromosome(&self) -> bool {
        // NCBI style: 1-22, X, Y
        // UCSC style: chr1-chr22, chrX, chrY
        matches!(
            self.name.as_str(),
            "1" | "2"
                | "3"
                | "4"
                | "5"
                | "6"
                | "7"
                | "8"
                | "9"
                | "10"
                | "11"
                | "12"
                | "13"
                | "14"
                | "15"
                | "16"
                | "17"
                | "18"
                | "19"
                | "20"
                | "21"
                | "22"
                | "X"
                | "Y"
                | "chr1"
                | "chr2"
                | "chr3"
                | "chr4"
                | "chr5"
                | "chr6"
                | "chr7"
                | "chr8"
                | "chr9"
                | "chr10"
                | "chr11"
                | "chr12"
                | "chr13"
                | "chr14"
                | "chr15"
                | "chr16"
                | "chr17"
                | "chr18"
                | "chr19"
                | "chr20"
                | "chr21"
                | "chr22"
                | "chrX"
                | "chrY"
        )
    }

    /// Check if this is a mitochondrial contig
    /// Matches common mitochondrial names from various reference builds
    pub fn is_mitochondrial(&self) -> bool {
        let name_lower = self.name.to_lowercase();
        matches!(
            name_lower.as_str(),
            "mt" | "m" | "chrm" | "chrmt" | "mito" | "mitochondrion" | "rcrs" | "nc_012920.1"
        ) || name_lower.contains("mitochon")
    }

    /// Check if this is an ALT contig (GRCh38)
    pub fn is_alt(&self) -> bool {
        self.name.ends_with("_alt") || self.name.contains("_alt_")
    }

    /// Check if this is a decoy contig
    pub fn is_decoy(&self) -> bool {
        self.name.contains("decoy")
            || self.name == "hs37d5"
            || self.name.starts_with("chrUn_")
            || self.name.contains("_random")
    }
}

// NOTE: normalize_contig_name() was removed.
// Name equivalence is now defined ONLY through explicit aliases (AN tag in SAM/dict,
// or NCBI assembly report columns). Matching uses exact names.

/// Detect the naming convention used by a set of contigs
pub fn detect_naming_convention(contigs: &[Contig]) -> NamingConvention {
    let mut has_chr_prefix = false;
    let mut has_no_prefix = false;

    for contig in contigs {
        if contig.is_primary_chromosome() {
            if contig.name.starts_with("chr") {
                has_chr_prefix = true;
            } else {
                has_no_prefix = true;
            }
        }
    }

    match (has_chr_prefix, has_no_prefix) {
        (true, false) => NamingConvention::Ucsc,
        (false, true) => NamingConvention::Ncbi,
        _ => NamingConvention::Mixed,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_primary_chromosome() {
        assert!(Contig::new("chr1", 100).is_primary_chromosome());
        assert!(Contig::new("1", 100).is_primary_chromosome());
        assert!(Contig::new("chrX", 100).is_primary_chromosome());
        assert!(Contig::new("Y", 100).is_primary_chromosome());
        assert!(!Contig::new("chrM", 100).is_primary_chromosome());
        assert!(!Contig::new("chr1_random", 100).is_primary_chromosome());
    }

    #[test]
    fn test_is_mitochondrial() {
        // Standard names
        assert!(Contig::new("chrM", 100).is_mitochondrial());
        assert!(Contig::new("MT", 100).is_mitochondrial());
        assert!(Contig::new("chrMT", 100).is_mitochondrial());
        assert!(Contig::new("M", 100).is_mitochondrial());
        // Extended names from older references
        assert!(Contig::new("mito", 100).is_mitochondrial());
        assert!(Contig::new("Mitochondrion", 100).is_mitochondrial());
        assert!(Contig::new("rCRS", 100).is_mitochondrial());
        assert!(Contig::new("NC_012920.1", 100).is_mitochondrial());
        // Substring match
        assert!(Contig::new("mitochondrial_genome", 100).is_mitochondrial());
        // Non-mitochondrial
        assert!(!Contig::new("chr1", 100).is_mitochondrial());
        assert!(!Contig::new("chrX", 100).is_mitochondrial());
    }
}
