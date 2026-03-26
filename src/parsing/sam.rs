use std::borrow::Cow;
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

/// Parse a BAM header from any reader (no file path required).
///
/// This enables parsing from in-memory buffers (e.g. `Cursor<Vec<u8>>`)
/// without writing to a temporary file. Only the header is read;
/// the reader does not need to contain complete record data.
///
/// # Errors
///
/// Returns `ParseError::Noodles` if the BAM header cannot be parsed,
/// `ParseError::InvalidFormat` if no contigs are found, or
/// `ParseError::TooManyContigs` if the limit is exceeded.
// TODO: remove #[allow(dead_code)] when wired into web server (Task 2)
#[allow(dead_code)]
pub fn parse_bam_from_reader<R: std::io::Read>(reader: R) -> Result<QueryHeader, ParseError> {
    use noodles::bam;

    let mut reader = bam::io::Reader::new(reader);

    let header = reader
        .read_header()
        .map_err(|e| ParseError::Noodles(e.to_string()))?;

    header_to_query(&header, None)
}

/// Parse a CRAM header from any reader (no file path required).
///
/// This enables parsing from in-memory buffers (e.g. `Cursor<Vec<u8>>`)
/// without writing to a temporary file. Only the file definition and
/// header are read.
///
/// # Errors
///
/// Returns `ParseError::Noodles` if the CRAM header cannot be parsed,
/// `ParseError::InvalidFormat` if no contigs are found, or
/// `ParseError::TooManyContigs` if the limit is exceeded.
// TODO: remove #[allow(dead_code)] when wired into web server (Task 2)
#[allow(dead_code)]
pub fn parse_cram_from_reader<R: std::io::Read>(reader: R) -> Result<QueryHeader, ParseError> {
    use noodles::cram;

    let mut reader = cram::io::Reader::new(reader);

    reader
        .read_file_definition()
        .map_err(|e| ParseError::Noodles(e.to_string()))?;

    let header = reader
        .read_file_header()
        .map_err(|e| ParseError::Noodles(e.to_string()))?;

    header_to_query(&header, None)
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
/// normalization was performed. Uses `Cow` to avoid allocation when no
/// normalization is needed (the common case for well-formed input).
#[must_use]
pub fn normalize_sam_whitespace(text: &str) -> (Cow<'_, str>, bool) {
    // Fast path: check if any line needs normalization before allocating
    if !text.lines().any(needs_space_to_tab_normalization) {
        return (Cow::Borrowed(text), false);
    }

    let mut normalized = String::with_capacity(text.len());

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
        } else {
            normalized.push_str(line);
        }
    }

    // Preserve trailing newline if present
    if text.ends_with('\n') {
        normalized.push('\n');
    }

    (Cow::Owned(normalized), true)
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
    let (normalized_text, _) = normalize_sam_whitespace(text);
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

#[cfg(test)]
mod reader_tests {
    use super::*;
    use std::io::Cursor;

    /// Helper: create a minimal BAM file in memory with the given contigs.
    fn create_test_bam(contigs: &[(&str, usize)]) -> Vec<u8> {
        use noodles::bam;
        use noodles::sam;
        use noodles::sam::header::record::value::map::ReferenceSequence;
        use noodles::sam::header::record::value::Map;
        use std::num::NonZeroUsize;

        let mut header = sam::Header::builder();
        for &(name, length) in contigs {
            header = header.add_reference_sequence(
                name,
                Map::<ReferenceSequence>::new(NonZeroUsize::new(length).unwrap()),
            );
        }
        let header = header.build();

        let mut buf = Vec::new();
        {
            let mut writer = bam::io::Writer::new(&mut buf);
            writer.write_header(&header).unwrap();
        }
        buf
    }

    #[test]
    fn test_parse_bam_from_reader_basic() {
        let bam_bytes = create_test_bam(&[("chr1", 248_956_422), ("chr2", 242_193_529)]);
        let cursor = Cursor::new(&bam_bytes);
        let query = parse_bam_from_reader(cursor).unwrap();
        assert_eq!(query.contigs.len(), 2);
        assert_eq!(query.contigs[0].name, "chr1");
        assert_eq!(query.contigs[0].length, 248_956_422);
        assert_eq!(query.contigs[1].name, "chr2");
        assert_eq!(query.contigs[1].length, 242_193_529);
    }

    #[test]
    fn test_parse_bam_from_reader_truncated_after_header() {
        let mut bam_bytes = create_test_bam(&[("chr1", 248_956_422)]);
        bam_bytes.extend_from_slice(&[0u8; 100]);
        let cursor = Cursor::new(&bam_bytes);
        let query = parse_bam_from_reader(cursor).unwrap();
        assert_eq!(query.contigs.len(), 1);
        assert_eq!(query.contigs[0].name, "chr1");
    }

    #[test]
    fn test_parse_bam_from_reader_empty() {
        let cursor = Cursor::new(Vec::<u8>::new());
        let result = parse_bam_from_reader(cursor);
        assert!(result.is_err());
    }

    /// Helper: create a minimal CRAM file in memory with the given contigs.
    ///
    /// Each contig is `(name, length, md5)`.  The CRAM writer requires MD5
    /// checksums on every reference sequence, so we must supply them here.
    fn create_test_cram(contigs: &[(&str, usize, &str)]) -> Vec<u8> {
        use noodles::cram;
        use noodles::sam;
        use noodles::sam::header::record::value::map::reference_sequence::tag;
        use noodles::sam::header::record::value::map::ReferenceSequence;
        use noodles::sam::header::record::value::Map;
        use std::num::NonZeroUsize;

        let mut header = sam::Header::builder();
        for &(name, length, md5) in contigs {
            let map = Map::<ReferenceSequence>::builder()
                .set_length(NonZeroUsize::new(length).unwrap())
                .insert(tag::MD5_CHECKSUM, md5)
                .build()
                .unwrap();
            header = header.add_reference_sequence(name, map);
        }
        let header = header.build();

        let mut buf = Vec::new();
        {
            let mut writer = cram::io::Writer::new(&mut buf);
            writer.write_file_definition().unwrap();
            writer.write_file_header(&header).unwrap();
        }
        buf
    }

    #[test]
    fn test_parse_cram_from_reader_basic() {
        let cram_bytes = create_test_cram(&[
            ("chr1", 248_956_422, "6aef897c3d6ff0c78aff06ac189178dd"),
            ("chrX", 156_040_895, "01234567890abcdef01234567890abcd"),
        ]);
        let cursor = Cursor::new(&cram_bytes);
        let query = parse_cram_from_reader(cursor).unwrap();
        assert_eq!(query.contigs.len(), 2);
        assert_eq!(query.contigs[0].name, "chr1");
        assert_eq!(query.contigs[0].length, 248_956_422);
        assert_eq!(
            query.contigs[0].md5,
            Some("6aef897c3d6ff0c78aff06ac189178dd".to_string())
        );
        assert_eq!(query.contigs[1].name, "chrX");
        assert_eq!(query.contigs[1].length, 156_040_895);
        assert_eq!(
            query.contigs[1].md5,
            Some("01234567890abcdef01234567890abcd".to_string())
        );
    }

    #[test]
    fn test_parse_cram_from_reader_truncated_after_header() {
        let mut cram_bytes =
            create_test_cram(&[("chr1", 248_956_422, "6aef897c3d6ff0c78aff06ac189178dd")]);
        // Append garbage bytes to simulate a truncated upload that ends mid-record.
        // parse_cram_from_reader only reads the header, so this must still succeed.
        cram_bytes.extend_from_slice(&[0u8; 100]);
        let cursor = Cursor::new(&cram_bytes);
        let query = parse_cram_from_reader(cursor).unwrap();
        assert_eq!(query.contigs.len(), 1);
        assert_eq!(query.contigs[0].name, "chr1");
    }

    #[test]
    fn test_parse_cram_from_reader_empty() {
        let cursor = Cursor::new(Vec::<u8>::new());
        let result = parse_cram_from_reader(cursor);
        assert!(result.is_err());
    }
}
