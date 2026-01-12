//! Parser for NCBI assembly report files.
//!
//! NCBI assembly reports contain rich metadata about contigs including multiple
//! naming conventions. The key columns are:
//!
//! - Sequence-Name: The primary name (e.g., "1", "X", "MT")
//! - GenBank-Accn: GenBank accession (e.g., "CM000663.2")
//! - RefSeq-Accn: RefSeq accession (e.g., "NC_000001.11")
//! - UCSC-style-name: UCSC-style name (e.g., "chr1")
//! - Sequence-Length: Length in base pairs
//!
//! All non-empty names from these columns become aliases for matching.

use std::collections::HashMap;

use crate::core::contig::{Contig, SequenceRole};
use crate::parsing::sam::ParseError;
use crate::utils::validation::check_contig_limit;

/// A parsed contig from an NCBI assembly report with all naming variants
#[derive(Debug, Clone)]
pub struct NcbiContigEntry {
    /// Primary sequence name (from Sequence-Name column)
    pub sequence_name: String,
    /// Sequence length
    pub length: u64,
    /// GenBank accession
    pub genbank_accn: Option<String>,
    /// RefSeq accession
    pub refseq_accn: Option<String>,
    /// UCSC-style name
    pub ucsc_name: Option<String>,
    /// Sequence role (e.g., "assembled-molecule", "unlocalized-scaffold")
    pub role: Option<String>,
}

impl NcbiContigEntry {
    /// Get all unique, non-empty names for this contig as aliases
    pub fn all_names(&self) -> Vec<String> {
        let mut names = vec![self.sequence_name.clone()];

        // Add other naming variants as aliases
        if let Some(ref name) = self.genbank_accn {
            if !name.is_empty() && name != "na" && !names.contains(name) {
                names.push(name.clone());
            }
        }
        if let Some(ref name) = self.refseq_accn {
            if !name.is_empty() && name != "na" && !names.contains(name) {
                names.push(name.clone());
            }
        }
        if let Some(ref name) = self.ucsc_name {
            if !name.is_empty() && name != "na" && !names.contains(name) {
                names.push(name.clone());
            }
        }

        names
    }

    /// Convert to a Contig with aliases populated
    pub fn to_contig(&self) -> Contig {
        let mut contig = Contig::new(&self.sequence_name, self.length);

        // All names except the primary become aliases
        let all = self.all_names();
        if all.len() > 1 {
            contig.aliases = all.into_iter().skip(1).collect();
        }

        // Set sequence role if available
        if let Some(ref role_str) = self.role {
            contig.sequence_role = SequenceRole::parse(role_str);
        }

        contig
    }
}

/// Parse NCBI assembly report from text
pub fn parse_ncbi_report_text(text: &str) -> Result<Vec<NcbiContigEntry>, ParseError> {
    let mut entries = Vec::new();
    // Use lowercase keys for case-insensitive matching
    let mut header_map: HashMap<String, usize> = HashMap::new();
    let mut found_header = false;

    for line in text.lines() {
        // Skip comment lines except the header
        if line.starts_with('#') {
            // The header line starts with "# " and contains column names
            // Use case-insensitive detection
            let line_lower = line.to_lowercase();
            if line_lower.contains("sequence-name") {
                let header_line = line.trim_start_matches('#').trim();
                for (idx, col) in header_line.split('\t').enumerate() {
                    // Store column names in lowercase for case-insensitive lookup
                    header_map.insert(col.trim().to_lowercase(), idx);
                }
                found_header = true;
            }
            continue;
        }

        if line.trim().is_empty() {
            continue;
        }

        if !found_header {
            return Err(ParseError::InvalidFormat(
                "NCBI assembly report header not found".to_string(),
            ));
        }

        let fields: Vec<&str> = line.split('\t').collect();

        // Get required fields (case-insensitive)
        let seq_name_idx = header_map
            .get("sequence-name")
            .ok_or_else(|| ParseError::InvalidFormat("Missing Sequence-Name column".to_string()))?;
        let length_idx = header_map.get("sequence-length").ok_or_else(|| {
            ParseError::InvalidFormat("Missing Sequence-Length column".to_string())
        })?;

        if fields.len() <= *seq_name_idx || fields.len() <= *length_idx {
            continue; // Skip malformed lines
        }

        let sequence_name = fields[*seq_name_idx].trim().to_string();
        let length: u64 = fields[*length_idx].trim().parse().map_err(|_| {
            ParseError::InvalidFormat(format!(
                "Invalid length for '{}': {}",
                sequence_name, fields[*length_idx]
            ))
        })?;

        // Get optional fields
        let get_optional = |name: &str| -> Option<String> {
            header_map
                .get(name)
                .and_then(|&idx| {
                    fields.get(idx).map(|s| {
                        let s = s.trim();
                        if s.is_empty() || s == "na" {
                            None
                        } else {
                            Some(s.to_string())
                        }
                    })
                })
                .flatten()
        };

        // Check contig limit for DOS protection
        if check_contig_limit(entries.len()).is_some() {
            return Err(ParseError::TooManyContigs(entries.len()));
        }

        entries.push(NcbiContigEntry {
            sequence_name,
            length,
            // Column names are lowercase in header_map for case-insensitive matching
            genbank_accn: get_optional("genbank-accn"),
            refseq_accn: get_optional("refseq-accn"),
            ucsc_name: get_optional("ucsc-style-name"),
            role: get_optional("sequence-role"),
        });
    }

    if entries.is_empty() {
        return Err(ParseError::InvalidFormat(
            "No contigs found in NCBI assembly report".to_string(),
        ));
    }

    Ok(entries)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_ncbi_report() {
        // Simplified NCBI assembly report format
        let report = r#"# Assembly name:  GRCh38.p14
# Organism name:  Homo sapiens
# Sequence-Name	Sequence-Role	Assigned-Molecule	Assigned-Molecule-Location/Type	GenBank-Accn	Relationship	RefSeq-Accn	Assembly-Unit	Sequence-Length	UCSC-style-name
1	assembled-molecule	1	Chromosome	CM000663.2	=	NC_000001.11	Primary Assembly	248956422	chr1
2	assembled-molecule	2	Chromosome	CM000664.2	=	NC_000002.12	Primary Assembly	242193529	chr2
MT	assembled-molecule	MT	Mitochondrion	J01415.2	=	NC_012920.1	non-nuclear	16569	chrM
"#;

        let entries = parse_ncbi_report_text(report).unwrap();
        assert_eq!(entries.len(), 3);

        // Check chr1
        let chr1 = &entries[0];
        assert_eq!(chr1.sequence_name, "1");
        assert_eq!(chr1.length, 248956422);
        assert_eq!(chr1.genbank_accn, Some("CM000663.2".to_string()));
        assert_eq!(chr1.refseq_accn, Some("NC_000001.11".to_string()));
        assert_eq!(chr1.ucsc_name, Some("chr1".to_string()));

        // Check all names returns aliases correctly
        let names = chr1.all_names();
        assert_eq!(names.len(), 4); // 1, CM000663.2, NC_000001.11, chr1
        assert!(names.contains(&"1".to_string()));
        assert!(names.contains(&"chr1".to_string()));
        assert!(names.contains(&"NC_000001.11".to_string()));

        // Check MT
        let mt = &entries[2];
        assert_eq!(mt.sequence_name, "MT");
        assert_eq!(mt.length, 16569);
        assert_eq!(mt.ucsc_name, Some("chrM".to_string()));
    }

    #[test]
    fn test_ncbi_entry_to_contig() {
        let entry = NcbiContigEntry {
            sequence_name: "1".to_string(),
            length: 248956422,
            genbank_accn: Some("CM000663.2".to_string()),
            refseq_accn: Some("NC_000001.11".to_string()),
            ucsc_name: Some("chr1".to_string()),
            role: Some("assembled-molecule".to_string()),
        };

        let contig = entry.to_contig();
        assert_eq!(contig.name, "1");
        assert_eq!(contig.length, 248956422);
        assert_eq!(contig.aliases.len(), 3); // CM000663.2, NC_000001.11, chr1
        assert!(contig.aliases.contains(&"chr1".to_string()));
        assert!(contig.aliases.contains(&"NC_000001.11".to_string()));
    }

    #[test]
    fn test_parse_ncbi_report_no_header() {
        let report = "1\tassembled-molecule\t1\t248956422\n";
        let result = parse_ncbi_report_text(report);
        assert!(result.is_err());
    }
}
