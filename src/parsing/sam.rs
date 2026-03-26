use std::io::BufReader;
use std::path::Path;
use thiserror::Error;
use tracing::warn;

use crate::core::contig::Contig;
use crate::core::header::QueryHeader;
use crate::utils::validation::{check_contig_limit, normalize_md5};

#[derive(Error, Debug)]
pub enum ParseError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Invalid SAM header format: {0}")]
    InvalidFormat(String),

    #[error("noodles error: {0}")]
    Noodles(String),

    #[error("Unsupported file format: {0}")]
    UnsupportedFormat(String),

    #[error("Too many contigs: {0} exceeds maximum allowed (100000)")]
    TooManyContigs(usize),
}

/// Parse a SAM/BAM/CRAM file and extract the header
///
/// # Errors
///
/// Returns `ParseError::Io` if the file cannot be read, `ParseError::Noodles` if
/// parsing fails, `ParseError::UnsupportedFormat` for unknown extensions,
/// `ParseError::InvalidFormat` if no contigs are found, or
/// `ParseError::TooManyContigs` if the limit is exceeded.
pub fn parse_file(path: &Path) -> Result<QueryHeader, ParseError> {
    let extension = path
        .extension()
        .and_then(|e| e.to_str())
        .map(str::to_lowercase);

    match extension.as_deref() {
        Some("sam") => parse_sam_file(path),
        Some("bam") => parse_bam_file(path),
        Some("cram") => parse_cram_file(path),
        Some(ext) => Err(ParseError::UnsupportedFormat(ext.to_string())),
        None => {
            // Try to detect from content - default to SAM
            parse_sam_file(path)
        }
    }
}

/// Parse a SAM file (text format)
fn parse_sam_file(path: &Path) -> Result<QueryHeader, ParseError> {
    use noodles::sam;

    let mut reader = std::fs::File::open(path)
        .map(BufReader::new)
        .map(sam::io::Reader::new)?;

    let header = reader
        .read_header()
        .map_err(|e| ParseError::Noodles(e.to_string()))?;

    header_to_query(&header, Some(path))
}

/// Parse a BAM file (binary format)
fn parse_bam_file(path: &Path) -> Result<QueryHeader, ParseError> {
    use noodles::bam;

    let mut reader = std::fs::File::open(path).map(bam::io::Reader::new)?;

    let header = reader
        .read_header()
        .map_err(|e| ParseError::Noodles(e.to_string()))?;

    header_to_query(&header, Some(path))
}

/// Parse a CRAM file
fn parse_cram_file(path: &Path) -> Result<QueryHeader, ParseError> {
    use noodles::cram;

    let mut reader = std::fs::File::open(path).map(cram::io::Reader::new)?;

    // Read file definition
    reader
        .read_file_definition()
        .map_err(|e| ParseError::Noodles(e.to_string()))?;

    let header = reader
        .read_file_header()
        .map_err(|e| ParseError::Noodles(e.to_string()))?;

    header_to_query(&header, Some(path))
}

/// Convert noodles header to `QueryHeader`
fn header_to_query(
    header: &noodles::sam::Header,
    source: Option<&Path>,
) -> Result<QueryHeader, ParseError> {
    use noodles::sam::header::record::value::map::tag::Other;

    let mut contigs = Vec::new();

    for (name, map) in header.reference_sequences() {
        let name_str = name.to_string();
        let length = map.length().get() as u64;

        let mut contig = Contig::new(name_str, length);

        // Extract M5 (MD5) tag from other_fields
        // The M5 tag contains the MD5 checksum of the sequence
        if let Ok(m5_tag) = Other::try_from(*b"M5") {
            if let Some(md5_value) = map.other_fields().get(&m5_tag) {
                let md5_str = md5_value.to_string();
                // Validate and normalize MD5 using centralized helper
                if let Some(normalized) = normalize_md5(&md5_str) {
                    contig.md5 = Some(normalized);
                } else {
                    warn!(
                        contig = %contig.name,
                        md5 = %md5_str,
                        "Invalid MD5 checksum format, ignoring"
                    );
                }
            }
        }

        // Extract AS (Assembly) tag
        if let Ok(as_tag) = Other::try_from(*b"AS") {
            if let Some(assembly_value) = map.other_fields().get(&as_tag) {
                contig.assembly = Some(assembly_value.to_string());
            }
        }

        // Extract UR (URI) tag
        if let Ok(ur_tag) = Other::try_from(*b"UR") {
            if let Some(uri_value) = map.other_fields().get(&ur_tag) {
                contig.uri = Some(uri_value.to_string());
            }
        }

        // Extract SP (Species) tag
        if let Ok(sp_tag) = Other::try_from(*b"SP") {
            if let Some(species_value) = map.other_fields().get(&sp_tag) {
                contig.species = Some(species_value.to_string());
            }
        }

        // Extract AN (Alternate Names) tag - comma-separated list of aliases
        if let Ok(an_tag) = Other::try_from(*b"AN") {
            if let Some(aliases_value) = map.other_fields().get(&an_tag) {
                let aliases: Vec<String> = aliases_value
                    .to_string()
                    .split(',')
                    .map(|s| s.trim().to_string())
                    .filter(|s| !s.is_empty())
                    .collect();
                if !aliases.is_empty() {
                    contig.aliases = aliases;
                }
            }
        }

        // Check contig limit for DOS protection
        if check_contig_limit(contigs.len()).is_some() {
            return Err(ParseError::TooManyContigs(contigs.len()));
        }

        contigs.push(contig);
    }

    let mut query = QueryHeader::new(contigs);
    if let Some(path) = source {
        query = query.with_source(path.display().to_string());
    }

    Ok(query)
}

/// Normalize SAM header lines that use spaces instead of tabs.
///
/// Browsers and copy-paste often convert tabs to spaces. This function detects
/// SAM header lines (`@XX` prefix) where fields are space-separated instead of
/// tab-separated and converts the spaces to tabs.
///
/// Only normalizes lines that start with a SAM header record type (`@HD`, `@SQ`,
/// `@RG`, `@PG`) followed by a space and a TAG: pattern. Lines that already
/// contain tabs are left unchanged. `@CO` (comment) lines are not normalized
/// because their content is free-form text where spaces are meaningful.
///
/// Returns the (possibly normalized) text and a boolean indicating whether any
/// normalization was performed.
#[must_use]
pub fn normalize_sam_whitespace(text: &str) -> (String, bool) {
    let mut normalized = String::with_capacity(text.len());
    let mut any_normalized = false;

    for line in text.lines() {
        if !normalized.is_empty() {
            normalized.push('\n');
        }

        if needs_space_to_tab_normalization(line) {
            let mut first = true;
            for field in line.split_whitespace() {
                if first {
                    normalized.push_str(field);
                    first = false;
                } else {
                    normalized.push('\t');
                    normalized.push_str(field);
                }
            }
            any_normalized = true;
        } else {
            normalized.push_str(line);
        }
    }

    // Preserve trailing newline if present
    if text.ends_with('\n') {
        normalized.push('\n');
    }

    (normalized, any_normalized)
}

/// Check if a line is a SAM header line that uses spaces instead of tabs.
///
/// Returns true only when the line starts with a recognized SAM record type
/// followed by whitespace and a TAG: pattern, and the line contains NO tab
/// characters.
fn needs_space_to_tab_normalization(line: &str) -> bool {
    if line.contains('\t') {
        return false;
    }

    // @CO lines are free-form comments — do not normalize
    let sam_prefixes = ["@HD ", "@SQ ", "@RG ", "@PG "];
    if !sam_prefixes.iter().any(|p| line.starts_with(p)) {
        return false;
    }

    // Must have at least one TAG:VALUE pattern after the record type
    line.split_whitespace().skip(1).any(|field| {
        field.len() >= 3
            && field.as_bytes().get(2) == Some(&b':')
            && field.as_bytes()[0].is_ascii_uppercase()
            && field.as_bytes()[1].is_ascii_uppercase()
    })
}

/// Parse header from raw text (stdin or pasted)
///
/// # Errors
///
/// Returns `ParseError::InvalidFormat` if the text has invalid format, missing
/// required fields, or no contigs are found, or `ParseError::TooManyContigs`
/// if the limit is exceeded.
pub fn parse_header_text(text: &str) -> Result<QueryHeader, ParseError> {
    let (normalized_text, was_normalized) = normalize_sam_whitespace(text);
    if was_normalized {
        warn!("SAM header contained spaces instead of tabs between fields; auto-corrected");
    }
    let text = &normalized_text;
    let mut contigs = Vec::new();

    for line in text.lines() {
        if !line.starts_with("@SQ") {
            continue;
        }

        let mut name: Option<String> = None;
        let mut length: Option<u64> = None;
        let mut md5_raw: Option<String> = None;
        let mut assembly: Option<String> = None;
        let mut uri: Option<String> = None;
        let mut species: Option<String> = None;
        let mut aliases: Vec<String> = Vec::new();

        for field in line.split('\t').skip(1) {
            if let Some((tag, value)) = field.split_once(':') {
                match tag {
                    "SN" => name = Some(value.to_string()),
                    "LN" => length = value.parse().ok(),
                    "M5" => md5_raw = Some(value.to_string()),
                    "AS" => assembly = Some(value.to_string()),
                    "UR" => uri = Some(value.to_string()),
                    "SP" => species = Some(value.to_string()),
                    "AN" => {
                        // Alternate names (aliases), comma-separated
                        aliases = value
                            .split(',')
                            .map(|s| s.trim().to_string())
                            .filter(|s| !s.is_empty())
                            .collect();
                    }
                    _ => {}
                }
            }
        }

        if let (Some(ref name_str), Some(length)) = (&name, length) {
            // Check contig limit for DOS protection
            if check_contig_limit(contigs.len()).is_some() {
                return Err(ParseError::TooManyContigs(contigs.len()));
            }

            // Validate and normalize MD5, warn if invalid
            let md5 = if let Some(ref raw) = md5_raw {
                if let Some(normalized) = normalize_md5(raw) {
                    Some(normalized)
                } else {
                    warn!(
                        contig = %name_str,
                        md5 = %raw,
                        "Invalid MD5 checksum format, ignoring"
                    );
                    None
                }
            } else {
                None
            };

            let mut contig = Contig::new(name_str.clone(), length);
            contig.md5 = md5;
            contig.assembly = assembly;
            contig.uri = uri;
            contig.species = species;
            contig.aliases = aliases;
            contigs.push(contig);
        }
    }

    if contigs.is_empty() {
        return Err(ParseError::InvalidFormat(
            "No @SQ lines found in header".to_string(),
        ));
    }

    Ok(QueryHeader::new(contigs))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_header_text() {
        let header = r"@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:248956422	M5:6aef897c3d6ff0c78aff06ac189178dd
@SQ	SN:chr2	LN:242193529	M5:f98db672eb0993dcfdabafe2a882905c
@SQ	SN:chrM	LN:16569
@RG	ID:sample1
";

        let query = parse_header_text(header).unwrap();
        assert_eq!(query.contigs.len(), 3);

        assert_eq!(query.contigs[0].name, "chr1");
        assert_eq!(query.contigs[0].length, 248_956_422);
        assert_eq!(
            query.contigs[0].md5,
            Some("6aef897c3d6ff0c78aff06ac189178dd".to_string())
        );

        assert_eq!(query.contigs[1].name, "chr2");
        assert_eq!(query.contigs[2].name, "chrM");
        assert!(query.contigs[2].md5.is_none());
    }

    #[test]
    fn test_parse_header_text_no_sq() {
        let header = "@HD\tVN:1.6\n@RG\tID:sample1\n";
        let result = parse_header_text(header);
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_header_text_with_aliases() {
        let header = r"@HD	VN:1.6
@SQ	SN:chr1	LN:248956422	M5:6aef897c3d6ff0c78aff06ac189178dd	AN:1,NC_000001.11
@SQ	SN:chrM	LN:16569	AN:MT,chrMT,NC_012920.1
";

        let query = parse_header_text(header).unwrap();
        assert_eq!(query.contigs.len(), 2);

        // Check chr1 aliases
        assert_eq!(query.contigs[0].name, "chr1");
        assert_eq!(
            query.contigs[0].aliases,
            vec!["1".to_string(), "NC_000001.11".to_string()]
        );

        // Check chrM aliases
        assert_eq!(query.contigs[1].name, "chrM");
        assert_eq!(
            query.contigs[1].aliases,
            vec![
                "MT".to_string(),
                "chrMT".to_string(),
                "NC_012920.1".to_string()
            ]
        );
    }

    #[test]
    fn test_normalize_sam_whitespace_spaces_to_tabs() {
        let input = "@SQ SN:chr1 LN:248956422 M5:6aef897c3d6ff0c78aff06ac189178dd\n";
        let (normalized, was_normalized) = normalize_sam_whitespace(input);
        assert!(was_normalized);
        assert_eq!(
            normalized,
            "@SQ\tSN:chr1\tLN:248956422\tM5:6aef897c3d6ff0c78aff06ac189178dd\n"
        );
    }

    #[test]
    fn test_normalize_sam_whitespace_already_tabs() {
        let input = "@SQ\tSN:chr1\tLN:248956422\n";
        let (normalized, was_normalized) = normalize_sam_whitespace(input);
        assert!(!was_normalized);
        assert_eq!(normalized, input);
    }

    #[test]
    fn test_normalize_sam_whitespace_mixed_lines() {
        let input =
            "@HD VN:1.6 SO:coordinate\n@SQ SN:chr1 LN:248956422\n@SQ SN:chr2 LN:242193529\n";
        let (normalized, was_normalized) = normalize_sam_whitespace(input);
        assert!(was_normalized);
        assert!(normalized.contains("@HD\tVN:1.6\tSO:coordinate"));
        assert!(normalized.contains("@SQ\tSN:chr1\tLN:248956422"));
        assert!(normalized.contains("@SQ\tSN:chr2\tLN:242193529"));
    }

    #[test]
    fn test_normalize_sam_whitespace_multiple_spaces() {
        let input = "@SQ  SN:chr1  LN:248956422\n";
        let (normalized, was_normalized) = normalize_sam_whitespace(input);
        assert!(was_normalized);
        assert_eq!(normalized, "@SQ\tSN:chr1\tLN:248956422\n");
    }

    #[test]
    fn test_normalize_sam_whitespace_preserves_non_header_lines() {
        let input = "some random text with spaces\n@SQ SN:chr1 LN:100\n";
        let (normalized, was_normalized) = normalize_sam_whitespace(input);
        assert!(was_normalized);
        assert!(normalized.starts_with("some random text with spaces\n"));
        assert!(normalized.contains("@SQ\tSN:chr1\tLN:100"));
    }

    #[test]
    fn test_normalize_sam_whitespace_tabs_and_spaces_mixed() {
        // Line has some tabs and some spaces — leave it alone
        let input = "@SQ\tSN:chr1 LN:248956422\n";
        let (normalized, was_normalized) = normalize_sam_whitespace(input);
        assert!(!was_normalized);
        assert_eq!(normalized, input);
    }

    #[test]
    fn test_normalize_sam_whitespace_skips_comment_lines() {
        let input = "@CO This is a comment with VN:1.0 mentioned\n@SQ SN:chr1 LN:100\n";
        let (normalized, was_normalized) = normalize_sam_whitespace(input);
        assert!(was_normalized);
        // Comment line should be preserved as-is
        assert!(normalized.starts_with("@CO This is a comment with VN:1.0 mentioned\n"));
        assert!(normalized.contains("@SQ\tSN:chr1\tLN:100"));
    }

    #[test]
    fn test_parse_header_text_with_spaces() {
        let header = "@SQ SN:chr1 LN:248956422 M5:6aef897c3d6ff0c78aff06ac189178dd\n\
                      @SQ SN:chr2 LN:242193529\n";
        let query = parse_header_text(header).unwrap();
        assert_eq!(query.contigs.len(), 2);
        assert_eq!(query.contigs[0].name, "chr1");
        assert_eq!(query.contigs[0].length, 248_956_422);
        assert_eq!(
            query.contigs[0].md5,
            Some("6aef897c3d6ff0c78aff06ac189178dd".to_string())
        );
        assert_eq!(query.contigs[1].name, "chr2");
        assert_eq!(query.contigs[1].length, 242_193_529);
    }
}
