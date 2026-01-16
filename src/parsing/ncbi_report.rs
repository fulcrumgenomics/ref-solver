//! Parser for NCBI assembly report files.
//!
//! NCBI assembly reports contain rich metadata about contigs including multiple
//! naming conventions. The key columns are:
//!
//! - Sequence-Name: The primary name (e.g., "1", "X", "MT")
//! - GenBank-Accn: `GenBank` accession (e.g., "CM000663.2")
//! - RefSeq-Accn: `RefSeq` accession (e.g., "`NC_000001.11`")
//! - UCSC-style-name: UCSC-style name (e.g., "chr1")
//! - Sequence-Length: Length in base pairs
//!
//! All non-empty names from these columns become aliases for matching.
//!
//! ## UCSC-style Name Generation for Patches
//!
//! For `GRCh38` assembly reports prior to p13, the UCSC-style-name column shows "na"
//! for fix-patches and novel-patches. However, UCSC does assign names to these
//! patches following a specific convention:
//!
//! - **Format**: `chr{chromosome}_{accession}v{version}_{suffix}`
//! - **Suffix**: `_fix` for fix-patches, `_alt` for novel-patches
//! - **Example**: `GenBank` accession `KN196472.1` on chromosome 1 as a fix-patch
//!   becomes `chr1_KN196472v1_fix`
//!
//! This module can optionally generate these UCSC-style names when they are missing
//! from the assembly report. This is controlled by the `generate_ucsc_names` parameter.
//!
//! ### Sources and References
//!
//! - UCSC FAQ on chromosome naming: <https://genome.ucsc.edu/FAQ/FAQdownloads.html>
//! - UCSC Patches blog post: <https://genome-blog.soe.ucsc.edu/blog/2019/02/22/patches/>
//! - UCSC hg38.p12 chrom.sizes: <https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p12/>
//! - GRC Patches documentation: <https://www.ncbi.nlm.nih.gov/grc/help/patches/>
//!
//! ### Verification
//!
//! The naming convention has been verified against UCSC's official chromosome size
//! files for GRCh38.p12 and cross-referenced with NCBI assembly reports for p12-p14.

use std::collections::HashMap;

use crate::core::contig::{Contig, SequenceRole};
use crate::parsing::sam::ParseError;
use crate::utils::validation::check_contig_limit;

/// Patch type for NCBI assembly report entries
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PatchType {
    /// A fix-patch corrects errors in the primary assembly
    Fix,
    /// A novel-patch adds new alternate sequence (labeled `_alt` by UCSC)
    Novel,
}

/// A parsed contig from an NCBI assembly report with all naming variants
#[derive(Debug, Clone)]
pub struct NcbiContigEntry {
    /// Primary sequence name (from Sequence-Name column)
    pub sequence_name: String,
    /// Sequence length
    pub length: u64,
    /// `GenBank` accession
    pub genbank_accn: Option<String>,
    /// `RefSeq` accession
    pub refseq_accn: Option<String>,
    /// UCSC-style name
    pub ucsc_name: Option<String>,
    /// Sequence role (e.g., "assembled-molecule", "unlocalized-scaffold")
    pub role: Option<String>,
    /// Assigned molecule (chromosome number, e.g., "1", "X", "Y", "MT")
    pub assigned_molecule: Option<String>,
    /// Patch type if this is a patch contig
    pub patch_type: Option<PatchType>,
}

/// Generate a UCSC-style name for a patch contig.
///
/// This function implements the UCSC naming convention for fix-patches and novel-patches:
/// - Format: `chr{chromosome}_{accession}v{version}_{suffix}`
/// - Where `suffix` is `fix` for fix-patches and `alt` for novel-patches
///
/// # Arguments
///
/// * `genbank_accession` - The `GenBank` accession (e.g., "KN196472.1")
/// * `chromosome` - The assigned chromosome (e.g., "1", "X", "Y")
/// * `patch_type` - Whether this is a fix-patch or novel-patch
///
/// # Returns
///
/// The generated UCSC-style name, or `None` if the accession format is invalid.
///
/// # Examples
///
/// ```
/// use ref_solver::parsing::ncbi_report::{generate_ucsc_patch_name, PatchType};
///
/// // Fix-patch example
/// let name = generate_ucsc_patch_name("KN196472.1", "1", PatchType::Fix);
/// assert_eq!(name, Some("chr1_KN196472v1_fix".to_string()));
///
/// // Novel-patch (alt) example
/// let name = generate_ucsc_patch_name("KQ458382.1", "1", PatchType::Novel);
/// assert_eq!(name, Some("chr1_KQ458382v1_alt".to_string()));
///
/// // Y chromosome fix-patch
/// let name = generate_ucsc_patch_name("KN196487.1", "Y", PatchType::Fix);
/// assert_eq!(name, Some("chrY_KN196487v1_fix".to_string()));
/// ```
#[must_use]
pub fn generate_ucsc_patch_name(
    genbank_accession: &str,
    chromosome: &str,
    patch_type: PatchType,
) -> Option<String> {
    // Parse the accession: expect format like "KN196472.1" -> ("KN196472", "1")
    let parts: Vec<&str> = genbank_accession.split('.').collect();
    if parts.len() != 2 {
        return None;
    }

    let accession_base = parts[0];
    let version = parts[1];

    // Validate that version is numeric
    if !version.chars().all(|c| c.is_ascii_digit()) {
        return None;
    }

    // Validate accession base is alphanumeric
    if !accession_base
        .chars()
        .all(|c| c.is_ascii_alphanumeric() || c == '_')
    {
        return None;
    }

    // Determine suffix based on patch type
    let suffix = match patch_type {
        PatchType::Fix => "fix",
        PatchType::Novel => "alt",
    };

    // Generate UCSC-style name: chr{chromosome}_{accession}v{version}_{suffix}
    Some(format!(
        "chr{chromosome}_{accession_base}v{version}_{suffix}"
    ))
}

impl NcbiContigEntry {
    /// Get all unique, non-empty names for this contig as aliases.
    ///
    /// # Arguments
    ///
    /// * `generate_ucsc_names` - If `true` and this is a patch contig without a
    ///   UCSC-style name in the assembly report, generate one using the UCSC
    ///   naming convention. This is useful for assembly reports prior to p13
    ///   where patch UCSC names were not included.
    ///
    /// # UCSC Name Generation
    ///
    /// When `generate_ucsc_names` is `true` and:
    /// - This contig is a fix-patch or novel-patch (determined from sequence role)
    /// - The UCSC-style-name column is "na" or missing
    /// - The `GenBank` accession and assigned molecule are available
    ///
    /// Then a UCSC-style name will be generated following the convention:
    /// `chr{chromosome}_{accession}v{version}_{suffix}`
    ///
    /// Where `suffix` is `_fix` for fix-patches and `_alt` for novel-patches.
    ///
    /// ## Sources
    ///
    /// - UCSC FAQ: <https://genome.ucsc.edu/FAQ/FAQdownloads.html>
    /// - UCSC Patches blog: <https://genome-blog.soe.ucsc.edu/blog/2019/02/22/patches/>
    /// - Verified against: <https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p12/>
    #[must_use]
    pub fn all_names_with_options(&self, generate_ucsc_names: bool) -> Vec<String> {
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

        // Handle UCSC name - either from assembly report or generated
        let effective_ucsc_name = if let Some(ref name) = self.ucsc_name {
            if !name.is_empty() && name != "na" {
                Some(name.clone())
            } else {
                None
            }
        } else {
            None
        };

        if let Some(ucsc_name) = effective_ucsc_name {
            if !names.contains(&ucsc_name) {
                names.push(ucsc_name);
            }
        } else if generate_ucsc_names {
            // Try to generate UCSC name for patches
            if let (Some(ref genbank), Some(ref chromosome), Some(patch_type)) =
                (&self.genbank_accn, &self.assigned_molecule, self.patch_type)
            {
                if !genbank.is_empty() && genbank != "na" && !chromosome.is_empty() {
                    if let Some(generated_name) =
                        generate_ucsc_patch_name(genbank, chromosome, patch_type)
                    {
                        if !names.contains(&generated_name) {
                            names.push(generated_name);
                        }
                    }
                }
            }
        }

        names
    }

    /// Convert to a Contig with aliases populated.
    ///
    /// This is equivalent to calling `to_contig_with_options(true)` - UCSC name
    /// generation is enabled by default.
    #[must_use]
    pub fn to_contig(&self) -> Contig {
        self.to_contig_with_options(true)
    }

    /// Convert to a Contig with aliases populated.
    ///
    /// # Arguments
    ///
    /// * `generate_ucsc_names` - If `true`, generate UCSC-style names for patches
    ///   that don't have them in the assembly report. See [`all_names_with_options`]
    ///   for details on the naming convention.
    #[must_use]
    pub fn to_contig_with_options(&self, generate_ucsc_names: bool) -> Contig {
        let mut contig = Contig::new(&self.sequence_name, self.length);

        // All names except the primary become aliases
        let all = self.all_names_with_options(generate_ucsc_names);
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
///
/// # Errors
///
/// Returns `ParseError::InvalidFormat` if the header is missing, required columns
/// are not found, or field values cannot be parsed.
#[allow(clippy::too_many_lines)]
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

        // Get optional fields (returns None for empty or "na")
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

        // Get raw optional fields (keeps "na" and empty values)
        let get_raw_optional = |name: &str| -> Option<String> {
            header_map
                .get(name)
                .and_then(|&idx| {
                    fields.get(idx).map(|s| {
                        let s = s.trim();
                        if s.is_empty() {
                            None
                        } else {
                            Some(s.to_string())
                        }
                    })
                })
                .flatten()
        };

        // Get sequence role to determine patch type
        let role = get_raw_optional("sequence-role");
        let patch_type = role.as_ref().and_then(|r| {
            let r_lower = r.to_lowercase();
            if r_lower == "fix-patch" {
                Some(PatchType::Fix)
            } else if r_lower == "novel-patch" {
                Some(PatchType::Novel)
            } else {
                None
            }
        });

        // Get assigned molecule (chromosome number)
        let assigned_molecule = get_optional("assigned-molecule");

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
            ucsc_name: get_raw_optional("ucsc-style-name"), // Keep "na" to detect missing names
            role,
            assigned_molecule,
            patch_type,
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
        let report = r"# Assembly name:  GRCh38.p14
# Organism name:  Homo sapiens
# Sequence-Name	Sequence-Role	Assigned-Molecule	Assigned-Molecule-Location/Type	GenBank-Accn	Relationship	RefSeq-Accn	Assembly-Unit	Sequence-Length	UCSC-style-name
1	assembled-molecule	1	Chromosome	CM000663.2	=	NC_000001.11	Primary Assembly	248956422	chr1
2	assembled-molecule	2	Chromosome	CM000664.2	=	NC_000002.12	Primary Assembly	242193529	chr2
MT	assembled-molecule	MT	Mitochondrion	J01415.2	=	NC_012920.1	non-nuclear	16569	chrM
";

        let entries = parse_ncbi_report_text(report).unwrap();
        assert_eq!(entries.len(), 3);

        // Check chr1
        let chr1 = &entries[0];
        assert_eq!(chr1.sequence_name, "1");
        assert_eq!(chr1.length, 248_956_422);
        assert_eq!(chr1.genbank_accn, Some("CM000663.2".to_string()));
        assert_eq!(chr1.refseq_accn, Some("NC_000001.11".to_string()));
        assert_eq!(chr1.ucsc_name, Some("chr1".to_string()));
        assert_eq!(chr1.assigned_molecule, Some("1".to_string()));
        assert!(chr1.patch_type.is_none()); // Not a patch

        // Check all names returns aliases correctly
        let names = chr1.all_names_with_options(true);
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
            length: 248_956_422,
            genbank_accn: Some("CM000663.2".to_string()),
            refseq_accn: Some("NC_000001.11".to_string()),
            ucsc_name: Some("chr1".to_string()),
            role: Some("assembled-molecule".to_string()),
            assigned_molecule: Some("1".to_string()),
            patch_type: None,
        };

        let contig = entry.to_contig();
        assert_eq!(contig.name, "1");
        assert_eq!(contig.length, 248_956_422);
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

    // ========================================================================
    // UCSC Name Generation Tests
    // ========================================================================

    #[test]
    fn test_generate_ucsc_patch_name_fix() {
        // Fix-patch on chromosome 1
        let name = generate_ucsc_patch_name("KN196472.1", "1", PatchType::Fix);
        assert_eq!(name, Some("chr1_KN196472v1_fix".to_string()));

        // Fix-patch on chromosome 8
        let name = generate_ucsc_patch_name("KZ208915.1", "8", PatchType::Fix);
        assert_eq!(name, Some("chr8_KZ208915v1_fix".to_string()));

        // Fix-patch on Y chromosome
        let name = generate_ucsc_patch_name("KN196487.1", "Y", PatchType::Fix);
        assert_eq!(name, Some("chrY_KN196487v1_fix".to_string()));
    }

    #[test]
    fn test_generate_ucsc_patch_name_novel() {
        // Novel-patch (alt) on chromosome 1
        let name = generate_ucsc_patch_name("KQ458382.1", "1", PatchType::Novel);
        assert_eq!(name, Some("chr1_KQ458382v1_alt".to_string()));

        // Novel-patch on X chromosome
        let name = generate_ucsc_patch_name("KV766199.1", "X", PatchType::Novel);
        assert_eq!(name, Some("chrX_KV766199v1_alt".to_string()));
    }

    #[test]
    fn test_generate_ucsc_patch_name_version_2() {
        // Accession with version 2
        let name = generate_ucsc_patch_name("GL000256.2", "6", PatchType::Novel);
        assert_eq!(name, Some("chr6_GL000256v2_alt".to_string()));
    }

    #[test]
    fn test_generate_ucsc_patch_name_invalid_format() {
        // Missing version
        let name = generate_ucsc_patch_name("KN196472", "1", PatchType::Fix);
        assert_eq!(name, None);

        // Too many dots
        let name = generate_ucsc_patch_name("KN196472.1.2", "1", PatchType::Fix);
        assert_eq!(name, None);

        // Non-numeric version
        let name = generate_ucsc_patch_name("KN196472.a", "1", PatchType::Fix);
        assert_eq!(name, None);

        // Empty accession
        let name = generate_ucsc_patch_name("", "1", PatchType::Fix);
        assert_eq!(name, None);
    }

    #[test]
    fn test_parse_ncbi_report_with_patches() {
        // Test parsing assembly report with fix-patch and novel-patch entries
        let report = r"# Assembly name:  GRCh38.p12
# Organism name:  Homo sapiens
# Sequence-Name	Sequence-Role	Assigned-Molecule	Assigned-Molecule-Location/Type	GenBank-Accn	Relationship	RefSeq-Accn	Assembly-Unit	Sequence-Length	UCSC-style-name
1	assembled-molecule	1	Chromosome	CM000663.2	=	NC_000001.11	Primary Assembly	248956422	chr1
HG986_PATCH	fix-patch	1	Chromosome	KN196472.1	=	NW_009646194.1	PATCHES	186494	na
HSCHR1_3_CTG3	novel-patch	1	Chromosome	KQ458382.1	=	NW_014040925.1	PATCHES	141019	na
";

        let entries = parse_ncbi_report_text(report).unwrap();
        assert_eq!(entries.len(), 3);

        // Check the fix-patch
        let fix_patch = &entries[1];
        assert_eq!(fix_patch.sequence_name, "HG986_PATCH");
        assert_eq!(fix_patch.patch_type, Some(PatchType::Fix));
        assert_eq!(fix_patch.assigned_molecule, Some("1".to_string()));
        assert_eq!(fix_patch.genbank_accn, Some("KN196472.1".to_string()));
        assert_eq!(fix_patch.ucsc_name, Some("na".to_string())); // "na" is kept

        // Check the novel-patch
        let novel_patch = &entries[2];
        assert_eq!(novel_patch.sequence_name, "HSCHR1_3_CTG3");
        assert_eq!(novel_patch.patch_type, Some(PatchType::Novel));
        assert_eq!(novel_patch.assigned_molecule, Some("1".to_string()));
    }

    #[test]
    fn test_all_names_with_ucsc_generation_enabled() {
        // Fix-patch without UCSC name in assembly report
        let entry = NcbiContigEntry {
            sequence_name: "HG986_PATCH".to_string(),
            length: 186_494,
            genbank_accn: Some("KN196472.1".to_string()),
            refseq_accn: Some("NW_009646194.1".to_string()),
            ucsc_name: Some("na".to_string()), // "na" means no UCSC name
            role: Some("fix-patch".to_string()),
            assigned_molecule: Some("1".to_string()),
            patch_type: Some(PatchType::Fix),
        };

        // With UCSC name generation enabled (default)
        let names = entry.all_names_with_options(true);
        assert!(
            names.contains(&"chr1_KN196472v1_fix".to_string()),
            "Generated UCSC name should be present: {names:?}"
        );
        assert!(names.contains(&"HG986_PATCH".to_string()));
        assert!(names.contains(&"KN196472.1".to_string()));
        assert!(names.contains(&"NW_009646194.1".to_string()));
    }

    #[test]
    fn test_all_names_with_ucsc_generation_disabled() {
        // Fix-patch without UCSC name in assembly report
        let entry = NcbiContigEntry {
            sequence_name: "HG986_PATCH".to_string(),
            length: 186_494,
            genbank_accn: Some("KN196472.1".to_string()),
            refseq_accn: Some("NW_009646194.1".to_string()),
            ucsc_name: Some("na".to_string()),
            role: Some("fix-patch".to_string()),
            assigned_molecule: Some("1".to_string()),
            patch_type: Some(PatchType::Fix),
        };

        // With UCSC name generation disabled
        let names = entry.all_names_with_options(false);
        assert!(
            !names.contains(&"chr1_KN196472v1_fix".to_string()),
            "Generated UCSC name should NOT be present: {names:?}"
        );
        assert!(names.contains(&"HG986_PATCH".to_string()));
        assert!(names.contains(&"KN196472.1".to_string()));
    }

    #[test]
    fn test_all_names_with_existing_ucsc_name() {
        // Entry that already has a UCSC name in the assembly report
        let entry = NcbiContigEntry {
            sequence_name: "1".to_string(),
            length: 248_956_422,
            genbank_accn: Some("CM000663.2".to_string()),
            refseq_accn: Some("NC_000001.11".to_string()),
            ucsc_name: Some("chr1".to_string()), // Has UCSC name
            role: Some("assembled-molecule".to_string()),
            assigned_molecule: Some("1".to_string()),
            patch_type: None,
        };

        // Both options should return the same result - use existing UCSC name
        let names_enabled = entry.all_names_with_options(true);
        let names_disabled = entry.all_names_with_options(false);

        assert!(names_enabled.contains(&"chr1".to_string()));
        assert!(names_disabled.contains(&"chr1".to_string()));
        assert_eq!(names_enabled.len(), names_disabled.len());
    }

    #[test]
    fn test_to_contig_with_ucsc_generation() {
        // Fix-patch without UCSC name
        let entry = NcbiContigEntry {
            sequence_name: "HG986_PATCH".to_string(),
            length: 186_494,
            genbank_accn: Some("KN196472.1".to_string()),
            refseq_accn: Some("NW_009646194.1".to_string()),
            ucsc_name: Some("na".to_string()),
            role: Some("fix-patch".to_string()),
            assigned_molecule: Some("1".to_string()),
            patch_type: Some(PatchType::Fix),
        };

        // With UCSC generation enabled
        let contig = entry.to_contig_with_options(true);
        assert!(
            contig.aliases.contains(&"chr1_KN196472v1_fix".to_string()),
            "Contig aliases should include generated UCSC name: {:?}",
            contig.aliases
        );

        // With UCSC generation disabled
        let contig = entry.to_contig_with_options(false);
        assert!(
            !contig.aliases.contains(&"chr1_KN196472v1_fix".to_string()),
            "Contig aliases should NOT include generated UCSC name: {:?}",
            contig.aliases
        );
    }

    #[test]
    fn test_novel_patch_ucsc_generation() {
        // Novel-patch should get "_alt" suffix
        let entry = NcbiContigEntry {
            sequence_name: "HSCHR1_3_CTG3".to_string(),
            length: 141_019,
            genbank_accn: Some("KQ458382.1".to_string()),
            refseq_accn: Some("NW_014040925.1".to_string()),
            ucsc_name: Some("na".to_string()),
            role: Some("novel-patch".to_string()),
            assigned_molecule: Some("1".to_string()),
            patch_type: Some(PatchType::Novel),
        };

        let names = entry.all_names_with_options(true);
        assert!(
            names.contains(&"chr1_KQ458382v1_alt".to_string()),
            "Novel patch should have _alt suffix: {names:?}"
        );
    }

    #[test]
    fn test_ucsc_generation_for_different_chromosomes() {
        // Test various chromosome values
        let test_cases = vec![
            ("1", "KN196472.1", PatchType::Fix, "chr1_KN196472v1_fix"),
            ("X", "KV766199.1", PatchType::Novel, "chrX_KV766199v1_alt"),
            ("Y", "KN196487.1", PatchType::Fix, "chrY_KN196487v1_fix"),
            ("22", "KZ208920.1", PatchType::Fix, "chr22_KZ208920v1_fix"),
        ];

        for (chrom, accession, patch_type, expected) in test_cases {
            let entry = NcbiContigEntry {
                sequence_name: "TEST_PATCH".to_string(),
                length: 1000,
                genbank_accn: Some(accession.to_string()),
                refseq_accn: None,
                ucsc_name: Some("na".to_string()),
                role: None,
                assigned_molecule: Some(chrom.to_string()),
                patch_type: Some(patch_type),
            };

            let names = entry.all_names_with_options(true);
            assert!(
                names.contains(&expected.to_string()),
                "Expected {expected} for chromosome {chrom}, accession {accession}: {names:?}"
            );
        }
    }
}
