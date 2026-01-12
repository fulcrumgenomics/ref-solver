//! Parser for VCF header contig lines.
//!
//! VCF files have contig definitions in the header as:
//! `##contig=<ID=chr1,length=248956422>`
//!
//! Additional fields like `md5` and `assembly` may also be present.
//!
//! Uses noodles for file parsing, with manual fallback for text parsing
//! to ensure all fields (including MD5) are properly extracted.

use std::path::Path;

use crate::core::contig::Contig;
use crate::core::header::QueryHeader;
use crate::parsing::sam::ParseError;
use crate::utils::validation::{check_contig_limit, normalize_md5};

/// Parse VCF file and extract contig definitions from header
pub fn parse_vcf_file(path: &Path) -> Result<QueryHeader, ParseError> {
    let content = std::fs::read_to_string(path)?;
    parse_vcf_header_text(&content)
}

/// Parse VCF header text and extract contig definitions
pub fn parse_vcf_header_text(text: &str) -> Result<QueryHeader, ParseError> {
    let mut contigs = Vec::new();

    for line in text.lines() {
        // VCF contig lines start with ##contig=
        if !line.starts_with("##contig=") {
            // Stop at the header line (starts with #CHROM)
            if line.starts_with("#CHROM") {
                break;
            }
            continue;
        }

        if let Some(contig) = parse_contig_line(line)? {
            // Check contig limit for DOS protection
            if check_contig_limit(contigs.len()).is_some() {
                return Err(ParseError::TooManyContigs(contigs.len()));
            }
            contigs.push(contig);
        }
    }

    if contigs.is_empty() {
        return Err(ParseError::InvalidFormat(
            "No ##contig lines found in VCF header".to_string(),
        ));
    }

    Ok(QueryHeader::new(contigs))
}

/// Parse a single ##contig=<...> line
fn parse_contig_line(line: &str) -> Result<Option<Contig>, ParseError> {
    // Format: ##contig=<ID=chr1,length=248956422,md5=abc123,...>
    let content = line
        .strip_prefix("##contig=<")
        .and_then(|s| s.strip_suffix('>'))
        .ok_or_else(|| {
            ParseError::InvalidFormat(format!("Invalid contig line format: {}", line))
        })?;

    let mut name: Option<String> = None;
    let mut length: Option<u64> = None;
    let mut md5: Option<String> = None;
    let mut assembly: Option<String> = None;

    // Parse key=value pairs, handling quoted values
    for part in split_contig_fields(content) {
        if let Some((key, value)) = part.split_once('=') {
            let key = key.trim();
            // Remove quotes from value if present
            let value = value.trim().trim_matches('"');

            match key.to_lowercase().as_str() {
                "id" => name = Some(value.to_string()),
                "length" => length = value.parse().ok(),
                "md5" => {
                    // Validate and normalize MD5 using centralized helper
                    md5 = normalize_md5(value);
                }
                "assembly" => assembly = Some(value.to_string()),
                _ => {}
            }
        }
    }

    match (name, length) {
        (Some(name), Some(length)) => {
            let mut contig = Contig::new(name, length);
            contig.md5 = md5;
            contig.assembly = assembly;
            Ok(Some(contig))
        }
        (Some(name), None) => Err(ParseError::InvalidFormat(format!(
            "Contig '{}' missing length",
            name
        ))),
        _ => Ok(None), // Skip malformed lines without ID
    }
}

/// Split contig fields, handling commas inside quoted values.
///
/// This is UTF-8 safe because:
/// - Commas are single-byte ASCII (0x2C)
/// - `char_indices()` yields byte positions at character boundaries
/// - After a comma at position `i`, `i + 1` is always a valid boundary
fn split_contig_fields(content: &str) -> Vec<&str> {
    let mut fields = Vec::new();
    let mut start = 0;
    let mut in_quotes = false;

    for (i, c) in content.char_indices() {
        match c {
            '"' => in_quotes = !in_quotes,
            ',' if !in_quotes => {
                fields.push(&content[start..i]);
                // Safe: comma is 1 byte, so i + 1 is a valid char boundary
                start = i + 1;
            }
            _ => {}
        }
    }

    // Don't forget the last field
    if start <= content.len() {
        fields.push(&content[start..]);
    }

    fields
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_vcf_header() {
        let vcf = r#"##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529,md5=f98db672eb0993dcfdabafe2a882905c>
##contig=<ID=chrM,length=16569,assembly=GRCh38>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
"#;

        let query = parse_vcf_header_text(vcf).unwrap();
        assert_eq!(query.contigs.len(), 3);

        assert_eq!(query.contigs[0].name, "chr1");
        assert_eq!(query.contigs[0].length, 248956422);
        assert!(query.contigs[0].md5.is_none());

        assert_eq!(query.contigs[1].name, "chr2");
        assert_eq!(query.contigs[1].length, 242193529);
        assert_eq!(
            query.contigs[1].md5,
            Some("f98db672eb0993dcfdabafe2a882905c".to_string())
        );

        assert_eq!(query.contigs[2].name, "chrM");
        assert_eq!(query.contigs[2].length, 16569);
        assert_eq!(query.contigs[2].assembly, Some("GRCh38".to_string()));
    }

    #[test]
    fn test_parse_vcf_no_contigs() {
        let vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\n";
        let result = parse_vcf_header_text(vcf);
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_contig_line() {
        let line = "##contig=<ID=chr1,length=248956422,md5=6aef897c3d6ff0c78aff06ac189178dd>";
        let contig = parse_contig_line(line).unwrap().unwrap();

        assert_eq!(contig.name, "chr1");
        assert_eq!(contig.length, 248956422);
        assert_eq!(
            contig.md5,
            Some("6aef897c3d6ff0c78aff06ac189178dd".to_string())
        );
    }

    #[test]
    fn test_parse_vcf_quoted_values() {
        // Test that quoted values are handled properly
        let vcf = r#"##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422,assembly="GRCh38.p14">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
"#;
        let query = parse_vcf_header_text(vcf).unwrap();
        assert_eq!(query.contigs.len(), 1);
        assert_eq!(query.contigs[0].name, "chr1");
        assert_eq!(query.contigs[0].assembly, Some("GRCh38.p14".to_string()));
    }

    #[test]
    fn test_split_contig_fields() {
        let fields = split_contig_fields(r#"ID=chr1,length=123,desc="foo,bar""#);
        assert_eq!(fields.len(), 3);
        assert_eq!(fields[0], "ID=chr1");
        assert_eq!(fields[1], "length=123");
        assert_eq!(fields[2], r#"desc="foo,bar""#);
    }

    #[test]
    fn test_split_contig_fields_utf8() {
        // Verify UTF-8 handling with multi-byte characters
        let fields = split_contig_fields("ID=chrα,length=123,desc=日本語");
        assert_eq!(fields.len(), 3);
        assert_eq!(fields[0], "ID=chrα");
        assert_eq!(fields[1], "length=123");
        assert_eq!(fields[2], "desc=日本語");
    }

    #[test]
    fn test_split_contig_fields_empty() {
        // Edge case: empty string
        let fields = split_contig_fields("");
        assert_eq!(fields.len(), 1);
        assert_eq!(fields[0], "");

        // Edge case: single field
        let fields = split_contig_fields("ID=chr1");
        assert_eq!(fields.len(), 1);
        assert_eq!(fields[0], "ID=chr1");
    }
}
