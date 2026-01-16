//! Parser for FASTA files using noodles.
//!
//! Extracts contig names and lengths from FASTA files.
//! Supports both uncompressed and gzip/bgzip compressed files.
//!
//! Supported extensions:
//! - `.fa`, `.fasta`, `.fna` (uncompressed)
//! - `.fa.gz`, `.fasta.gz`, `.fna.gz` (gzip compressed)
//! - `.fa.bgz`, `.fasta.bgz`, `.fna.bgz` (bgzip compressed)

use std::ffi::OsStr;
use std::io::{BufRead, BufReader};
use std::path::Path;

use flate2::read::GzDecoder;
use noodles::fasta;

use crate::core::contig::Contig;
use crate::core::header::QueryHeader;
use crate::parsing::sam::ParseError;
use crate::utils::validation::check_contig_limit;

/// Check if the path has a FASTA extension
pub fn is_fasta_file(path: &Path) -> bool {
    let path_str = path.to_string_lossy().to_lowercase();

    // Check for gzipped FASTA
    if path_str.ends_with(".fa.gz")
        || path_str.ends_with(".fasta.gz")
        || path_str.ends_with(".fna.gz")
        || path_str.ends_with(".fa.bgz")
        || path_str.ends_with(".fasta.bgz")
        || path_str.ends_with(".fna.bgz")
    {
        return true;
    }

    // Check for uncompressed FASTA
    matches!(
        path.extension()
            .and_then(OsStr::to_str)
            .map(str::to_lowercase)
            .as_deref(),
        Some("fa" | "fasta" | "fna")
    )
}

/// Check if the path is a gzipped file
#[allow(clippy::case_sensitive_file_extension_comparisons)] // Already lowercased
fn is_gzipped(path: &Path) -> bool {
    let path_str = path.to_string_lossy().to_lowercase();
    path_str.ends_with(".gz") || path_str.ends_with(".bgz")
}

/// Parse a FASTA file and extract contig names and lengths.
///
/// This reads through the entire FASTA file to determine sequence lengths.
/// For large files, consider using an existing .fai index instead.
///
/// # Errors
///
/// Returns `ParseError::Io` if the file cannot be read, `ParseError::Noodles` if
/// parsing fails, `ParseError::InvalidFormat` if no contigs are found, or
/// `ParseError::TooManyContigs` if the limit is exceeded.
pub fn parse_fasta_file(path: &Path) -> Result<QueryHeader, ParseError> {
    if is_gzipped(path) {
        parse_fasta_gzipped(path)
    } else {
        parse_fasta_uncompressed(path)
    }
}

/// Parse a FASTA file and compute MD5 checksums for each sequence.
///
/// This reads through the entire FASTA file and computes MD5 on the
/// uppercase sequence (standard for sequence checksums).
/// This is slower than `parse_fasta_file` but provides MD5 for matching.
///
/// # Errors
///
/// Returns `ParseError::Io` if the file cannot be read, `ParseError::Noodles` if
/// parsing fails, `ParseError::InvalidFormat` if no contigs are found, or
/// `ParseError::TooManyContigs` if the limit is exceeded.
pub fn parse_fasta_file_with_md5(path: &Path) -> Result<QueryHeader, ParseError> {
    if is_gzipped(path) {
        parse_fasta_gzipped_with_md5(path)
    } else {
        parse_fasta_uncompressed_with_md5(path)
    }
}

/// Parse an uncompressed FASTA file with MD5 computation
fn parse_fasta_uncompressed_with_md5(path: &Path) -> Result<QueryHeader, ParseError> {
    let file = std::fs::File::open(path)?;
    let reader = BufReader::new(file);
    let mut fasta_reader = fasta::io::Reader::new(reader);

    parse_fasta_reader_with_md5(&mut fasta_reader)
}

/// Parse a gzip-compressed FASTA file with MD5 computation
fn parse_fasta_gzipped_with_md5(path: &Path) -> Result<QueryHeader, ParseError> {
    let file = std::fs::File::open(path)?;
    let decoder = GzDecoder::new(file);
    let reader = BufReader::new(decoder);
    let mut fasta_reader = fasta::io::Reader::new(reader);

    parse_fasta_reader_with_md5(&mut fasta_reader)
}

/// Parse from a noodles FASTA reader with MD5 computation
fn parse_fasta_reader_with_md5<R: BufRead>(
    reader: &mut fasta::io::Reader<R>,
) -> Result<QueryHeader, ParseError> {
    let mut contigs = Vec::new();

    for result in reader.records() {
        let record = result
            .map_err(|e| ParseError::Noodles(format!("Failed to parse FASTA record: {e}")))?;

        // Check contig limit for DOS protection
        if check_contig_limit(contigs.len()).is_some() {
            return Err(ParseError::TooManyContigs(contigs.len()));
        }

        let name = String::from_utf8_lossy(record.name()).to_string();
        let sequence = record.sequence();
        let length = sequence.len() as u64;

        // Compute MD5 on uppercase sequence (standard convention)
        let uppercase: Vec<u8> = sequence
            .as_ref()
            .iter()
            .map(u8::to_ascii_uppercase)
            .collect();
        let md5 = format!("{:x}", md5::compute(&uppercase));

        let mut contig = Contig::new(name, length);
        contig.md5 = Some(md5);
        contigs.push(contig);
    }

    if contigs.is_empty() {
        return Err(ParseError::InvalidFormat(
            "No sequences found in FASTA file".to_string(),
        ));
    }

    Ok(QueryHeader::new(contigs))
}

/// Parse an uncompressed FASTA file
fn parse_fasta_uncompressed(path: &Path) -> Result<QueryHeader, ParseError> {
    let file = std::fs::File::open(path)?;
    let reader = BufReader::new(file);
    let mut fasta_reader = fasta::io::Reader::new(reader);

    parse_fasta_reader(&mut fasta_reader)
}

/// Parse a gzip-compressed FASTA file
fn parse_fasta_gzipped(path: &Path) -> Result<QueryHeader, ParseError> {
    let file = std::fs::File::open(path)?;
    let decoder = GzDecoder::new(file);
    let reader = BufReader::new(decoder);
    let mut fasta_reader = fasta::io::Reader::new(reader);

    parse_fasta_reader(&mut fasta_reader)
}

/// Parse from a noodles FASTA reader
fn parse_fasta_reader<R: BufRead>(
    reader: &mut fasta::io::Reader<R>,
) -> Result<QueryHeader, ParseError> {
    let mut contigs = Vec::new();

    for result in reader.records() {
        let record = result
            .map_err(|e| ParseError::Noodles(format!("Failed to parse FASTA record: {e}")))?;

        // Check contig limit for DOS protection
        if check_contig_limit(contigs.len()).is_some() {
            return Err(ParseError::TooManyContigs(contigs.len()));
        }

        let name = String::from_utf8_lossy(record.name()).to_string();
        let length = record.sequence().len() as u64;

        contigs.push(Contig::new(name, length));
    }

    if contigs.is_empty() {
        return Err(ParseError::InvalidFormat(
            "No sequences found in FASTA file".to_string(),
        ));
    }

    Ok(QueryHeader::new(contigs))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_is_fasta_file() {
        assert!(is_fasta_file(Path::new("test.fa")));
        assert!(is_fasta_file(Path::new("test.fasta")));
        assert!(is_fasta_file(Path::new("test.fna")));
        assert!(is_fasta_file(Path::new("test.fa.gz")));
        assert!(is_fasta_file(Path::new("test.fasta.gz")));
        assert!(is_fasta_file(Path::new("test.fna.bgz")));
        assert!(is_fasta_file(Path::new("/path/to/Reference.FA")));

        assert!(!is_fasta_file(Path::new("test.bam")));
        assert!(!is_fasta_file(Path::new("test.sam")));
        assert!(!is_fasta_file(Path::new("test.fai")));
    }

    #[test]
    fn test_parse_fasta_file() {
        let fasta_content = b">chr1 description\nACGTACGT\nACGT\n>chr2\nGGGG\n";

        let mut temp = NamedTempFile::with_suffix(".fa").unwrap();
        temp.write_all(fasta_content).unwrap();
        temp.flush().unwrap();

        let query = parse_fasta_file(temp.path()).unwrap();
        assert_eq!(query.contigs.len(), 2);
        assert_eq!(query.contigs[0].name, "chr1");
        assert_eq!(query.contigs[0].length, 12); // 8 + 4 bases
        assert_eq!(query.contigs[1].name, "chr2");
        assert_eq!(query.contigs[1].length, 4);
    }

    #[test]
    fn test_parse_empty_fasta() {
        let mut temp = NamedTempFile::with_suffix(".fa").unwrap();
        temp.write_all(b"").unwrap();
        temp.flush().unwrap();

        let result = parse_fasta_file(temp.path());
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_fasta_with_md5() {
        // Simple sequence to verify MD5 computation
        // "ACGT" uppercase -> MD5 = f1f8f4bf413b16ad135722aa4591043e
        let fasta_content = b">chr1\nACGT\n";

        let mut temp = NamedTempFile::with_suffix(".fa").unwrap();
        temp.write_all(fasta_content).unwrap();
        temp.flush().unwrap();

        let query = parse_fasta_file_with_md5(temp.path()).unwrap();
        assert_eq!(query.contigs.len(), 1);
        assert_eq!(query.contigs[0].name, "chr1");
        assert_eq!(query.contigs[0].length, 4);
        assert_eq!(
            query.contigs[0].md5,
            Some("f1f8f4bf413b16ad135722aa4591043e".to_string())
        );
    }

    #[test]
    fn test_parse_fasta_with_md5_lowercase() {
        // MD5 should be computed on uppercase, so "acgt" should give same result as "ACGT"
        let fasta_content = b">chr1\nacgt\n";

        let mut temp = NamedTempFile::with_suffix(".fa").unwrap();
        temp.write_all(fasta_content).unwrap();
        temp.flush().unwrap();

        let query = parse_fasta_file_with_md5(temp.path()).unwrap();
        assert_eq!(
            query.contigs[0].md5,
            Some("f1f8f4bf413b16ad135722aa4591043e".to_string())
        );
    }
}
