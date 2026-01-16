//! Parser for FASTA index (.fai) files using noodles.
//!
//! FAI format provides name and length for each contig, but no MD5 or aliases.
//! Format: `name\tlength\toffset\tline_bases\tline_width`

use std::io::BufReader;
use std::path::Path;

use crate::core::contig::Contig;
use crate::core::header::QueryHeader;
use crate::parsing::sam::ParseError;
use crate::utils::validation::check_contig_limit;

/// Parse a FASTA index (.fai) file using noodles
///
/// # Errors
///
/// Returns `ParseError::Io` if the file cannot be read, `ParseError::Noodles` if
/// parsing fails, `ParseError::InvalidFormat` if no contigs are found, or
/// `ParseError::TooManyContigs` if the limit is exceeded.
pub fn parse_fai_file(path: &Path) -> Result<QueryHeader, ParseError> {
    use noodles::fasta;

    let reader = std::fs::File::open(path).map(BufReader::new)?;

    let index = fasta::fai::io::Reader::new(reader)
        .read_index()
        .map_err(|e| ParseError::Noodles(format!("Failed to parse FAI file: {e}")))?;

    index_to_query(&index)
}

/// Convert noodles FAI index to `QueryHeader`
fn index_to_query(index: &noodles::fasta::fai::Index) -> Result<QueryHeader, ParseError> {
    let mut contigs = Vec::new();

    for record in index.as_ref() {
        // Check contig limit for DOS protection
        if check_contig_limit(contigs.len()).is_some() {
            return Err(ParseError::TooManyContigs(contigs.len()));
        }

        let name = String::from_utf8_lossy(record.name()).to_string();
        let length = record.length();

        contigs.push(Contig::new(name, length));
    }

    if contigs.is_empty() {
        return Err(ParseError::InvalidFormat(
            "No contigs found in FAI file".to_string(),
        ));
    }

    Ok(QueryHeader::new(contigs))
}

/// Parse FAI from text (fallback for raw text input)
///
/// # Errors
///
/// Returns `ParseError::InvalidFormat` if the text has invalid format or no contigs,
/// or `ParseError::TooManyContigs` if the limit is exceeded.
pub fn parse_fai_text(text: &str) -> Result<QueryHeader, ParseError> {
    let mut contigs = Vec::new();

    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 2 {
            continue;
        }

        // Check contig limit for DOS protection
        if check_contig_limit(contigs.len()).is_some() {
            return Err(ParseError::TooManyContigs(contigs.len()));
        }

        let name = fields[0].to_string();
        let length: u64 = fields[1].parse().map_err(|_| {
            ParseError::InvalidFormat(format!(
                "Invalid length for contig '{}': {}",
                name, fields[1]
            ))
        })?;

        contigs.push(Contig::new(name, length));
    }

    if contigs.is_empty() {
        return Err(ParseError::InvalidFormat(
            "No contigs found in FAI file".to_string(),
        ));
    }

    Ok(QueryHeader::new(contigs))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Parsed contig from FAI file with offset information
    #[derive(Debug, Clone)]
    pub struct FaiEntry {
        pub name: String,
        pub length: u64,
        pub offset: u64,
        pub line_bases: u32,
        pub line_width: u32,
    }

    /// Parse FAI file with full entry information
    pub fn parse_fai_entries(text: &str) -> Result<Vec<FaiEntry>, ParseError> {
        let mut entries = Vec::new();

        for line in text.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 5 {
                return Err(ParseError::InvalidFormat(format!(
                    "FAI line has {} fields, expected 5: {}",
                    fields.len(),
                    line
                )));
            }

            // Check limit
            if check_contig_limit(entries.len()).is_some() {
                return Err(ParseError::TooManyContigs(entries.len()));
            }

            let name = fields[0].to_string();
            let length: u64 = fields[1]
                .parse()
                .map_err(|_| ParseError::InvalidFormat(format!("Invalid length: {}", fields[1])))?;
            let offset: u64 = fields[2]
                .parse()
                .map_err(|_| ParseError::InvalidFormat(format!("Invalid offset: {}", fields[2])))?;
            let line_bases: u32 = fields[3].parse().map_err(|_| {
                ParseError::InvalidFormat(format!("Invalid line_bases: {}", fields[3]))
            })?;
            let line_width: u32 = fields[4].parse().map_err(|_| {
                ParseError::InvalidFormat(format!("Invalid line_width: {}", fields[4]))
            })?;

            entries.push(FaiEntry {
                name,
                length,
                offset,
                line_bases,
                line_width,
            });
        }

        if entries.is_empty() {
            return Err(ParseError::InvalidFormat(
                "No entries found in FAI file".to_string(),
            ));
        }

        Ok(entries)
    }

    #[test]
    fn test_parse_fai_text() {
        let fai = r"chr1	248956422	112	70	71
chr2	242193529	253404903	70	71
chrM	16569	3099922541	70	71
";

        let query = parse_fai_text(fai).unwrap();
        assert_eq!(query.contigs.len(), 3);

        assert_eq!(query.contigs[0].name, "chr1");
        assert_eq!(query.contigs[0].length, 248_956_422);
        assert!(query.contigs[0].md5.is_none()); // FAI doesn't have MD5

        assert_eq!(query.contigs[1].name, "chr2");
        assert_eq!(query.contigs[1].length, 242_193_529);

        assert_eq!(query.contigs[2].name, "chrM");
        assert_eq!(query.contigs[2].length, 16569);
    }

    #[test]
    fn test_parse_fai_entries() {
        let fai = "chr1\t248956422\t112\t70\t71\n";

        let entries = parse_fai_entries(fai).unwrap();
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].name, "chr1");
        assert_eq!(entries[0].length, 248_956_422);
        assert_eq!(entries[0].offset, 112);
        assert_eq!(entries[0].line_bases, 70);
        assert_eq!(entries[0].line_width, 71);
    }

    #[test]
    fn test_parse_fai_empty() {
        let result = parse_fai_text("");
        assert!(result.is_err());
    }
}
