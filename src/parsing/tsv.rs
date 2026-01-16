use std::path::Path;

use crate::core::contig::Contig;
use crate::core::header::QueryHeader;
use crate::parsing::sam::ParseError;
use crate::utils::validation::{check_contig_limit, normalize_md5};

/// Parse a TSV/CSV file with columns: name, length, [md5]
///
/// # Errors
///
/// Returns `ParseError::Io` if the file cannot be read, or other parse errors
/// if the content is invalid.
pub fn parse_tsv_file(path: &Path, delimiter: char) -> Result<QueryHeader, ParseError> {
    let content = std::fs::read_to_string(path)?;
    parse_tsv_text(&content, delimiter)
}

/// Parse TSV/CSV text with columns: name, length, [md5]
///
/// # Errors
///
/// Returns `ParseError::InvalidFormat` if lines have fewer than 2 fields,
/// contain invalid length values, or no contigs are found, or
/// `ParseError::TooManyContigs` if the limit is exceeded.
pub fn parse_tsv_text(text: &str, delimiter: char) -> Result<QueryHeader, ParseError> {
    let mut contigs = Vec::new();
    let mut first_data_line = true;

    for (i, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split(delimiter).collect();

        // Check if first non-empty/non-comment line is a header
        if first_data_line {
            first_data_line = false;
            let first = fields.first().map(|s| s.to_lowercase()).unwrap_or_default();
            if first == "name" || first == "sn" || first == "contig" || first == "chrom" {
                continue;
            }
        }

        // Line numbers in errors are 1-based for user friendliness
        let line_num = i + 1;

        if fields.len() < 2 {
            return Err(ParseError::InvalidFormat(format!(
                "Line {line_num} has fewer than 2 fields"
            )));
        }

        let name = fields[0].trim().to_string();
        let length: u64 = fields[1].trim().parse().map_err(|_| {
            ParseError::InvalidFormat(format!(
                "Invalid length on line {}: '{}'",
                line_num, fields[1]
            ))
        })?;

        let mut contig = Contig::new(name, length);

        // Optional MD5 in third column
        if fields.len() > 2 {
            // Validate and normalize MD5 using centralized helper
            contig.md5 = normalize_md5(fields[2].trim());
        }

        // Check contig limit for DOS protection
        if check_contig_limit(contigs.len()).is_some() {
            return Err(ParseError::TooManyContigs(contigs.len()));
        }

        contigs.push(contig);
    }

    if contigs.is_empty() {
        return Err(ParseError::InvalidFormat(
            "No contigs found in file".to_string(),
        ));
    }

    Ok(QueryHeader::new(contigs))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_tsv_text() {
        let tsv = r"name	length	md5
chr1	248956422	6aef897c3d6ff0c78aff06ac189178dd
chr2	242193529	f98db672eb0993dcfdabafe2a882905c
chrM	16569
";

        let query = parse_tsv_text(tsv, '\t').unwrap();
        assert_eq!(query.contigs.len(), 3);
        assert_eq!(query.contigs[0].name, "chr1");
        assert_eq!(query.contigs[0].length, 248_956_422);
        assert!(query.contigs[0].md5.is_some());
        assert!(query.contigs[2].md5.is_none());
    }

    #[test]
    fn test_parse_csv_text() {
        let csv = r"chrom,length,md5
chr1,248956422,6aef897c3d6ff0c78aff06ac189178dd
chr2,242193529,f98db672eb0993dcfdabafe2a882905c
";

        let query = parse_tsv_text(csv, ',').unwrap();
        assert_eq!(query.contigs.len(), 2);
    }

    #[test]
    fn test_parse_tsv_no_header() {
        let tsv = "chr1\t248956422\nchr2\t242193529\n";
        let query = parse_tsv_text(tsv, '\t').unwrap();
        assert_eq!(query.contigs.len(), 2);
    }

    #[test]
    fn test_parse_tsv_comments_before_header() {
        // Test that header detection works even with comments before it
        let tsv = r"# This is a comment
# Another comment

name	length	md5
chr1	248956422	6aef897c3d6ff0c78aff06ac189178dd
chr2	242193529	f98db672eb0993dcfdabafe2a882905c
";
        let query = parse_tsv_text(tsv, '\t').unwrap();
        assert_eq!(query.contigs.len(), 2);
        assert_eq!(query.contigs[0].name, "chr1");
    }
}
