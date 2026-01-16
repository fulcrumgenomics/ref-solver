use crate::core::header::QueryHeader;
use std::path::Path;

/// Supported file formats for reference identification
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FileFormat {
    /// SAM/BAM/CRAM files or plain text headers
    Sam,
    /// BAM binary format
    Bam,
    /// CRAM binary format
    Cram,
    /// Picard sequence dictionary files
    Dict,
    /// VCF files with contig headers
    Vcf,
    /// NCBI assembly report files
    NcbiReport,
    /// TSV/CSV tabular files
    Tsv,
    /// FASTA index (.fai) files
    Fai,
    /// FASTA sequence files
    Fasta,
    /// Automatically detect format
    Auto,
}

/// Errors that can occur during format detection
#[derive(Debug, PartialEq, thiserror::Error)]
pub enum FormatError {
    #[error("Unable to detect file format from content and filename")]
    UnknownFormat,
    #[error("File appears to be binary but cannot determine specific format")]
    UnsupportedBinary,
}

/// Errors that can occur during parsing
#[derive(Debug, thiserror::Error)]
pub enum ParseError {
    #[error("Failed to parse {format:?} content: {message}")]
    ParseFailed { format: FileFormat, message: String },
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
}

impl FileFormat {
    /// Get the display name for this format
    #[must_use]
    #[allow(clippy::trivially_copy_pass_by_ref)] // Idiomatic method signature
    pub fn display_name(&self) -> &'static str {
        match self {
            FileFormat::Sam => "SAM/BAM Header",
            FileFormat::Bam => "BAM File",
            FileFormat::Cram => "CRAM File",
            FileFormat::Dict => "Sequence Dictionary",
            FileFormat::Vcf => "VCF File",
            FileFormat::NcbiReport => "NCBI Assembly Report",
            FileFormat::Tsv => "TSV/CSV Table",
            FileFormat::Fai => "FASTA Index",
            FileFormat::Fasta => "FASTA File",
            FileFormat::Auto => "Auto-detect",
        }
    }
}

/// Detect file format from content and optional filename
///
/// # Errors
///
/// Returns `FormatError::UnknownFormat` if the format cannot be detected, or
/// `FormatError::UnsupportedBinary` if the file appears to be binary but the
/// specific format cannot be determined.
pub fn detect_format(content: &str, filename: Option<&str>) -> Result<FileFormat, FormatError> {
    // First try filename-based detection if available
    if let Some(name) = filename {
        if let Some(format) = detect_format_from_filename(name) {
            // For binary formats, trust filename-based detection without content validation
            if matches!(format, FileFormat::Bam | FileFormat::Cram) {
                return Ok(format);
            }
            // For text formats, validate that content matches expected format
            if validate_format_content(content, &format) {
                return Ok(format);
            }
        }
    }

    // Fall back to content-based detection
    detect_format_from_content(content)
}

/// Detect format based on filename and extension
fn detect_format_from_filename(filename: &str) -> Option<FileFormat> {
    let path = Path::new(filename);
    let lower_name = filename.to_lowercase();

    // Check for compressed formats first (multi-extension patterns)
    if lower_name.ends_with(".vcf.gz") {
        return Some(FileFormat::Vcf);
    }
    if lower_name.ends_with(".fa.gz")
        || lower_name.ends_with(".fasta.gz")
        || lower_name.ends_with(".fna.gz")
        || lower_name.ends_with(".fa.bgz")
        || lower_name.ends_with(".fasta.bgz")
        || lower_name.ends_with(".fna.bgz")
    {
        return Some(FileFormat::Fasta);
    }

    let extension = path.extension()?.to_str()?.to_lowercase();

    match extension.as_str() {
        "sam" => Some(FileFormat::Sam),
        "bam" => Some(FileFormat::Bam),
        "cram" => Some(FileFormat::Cram),
        "dict" => Some(FileFormat::Dict),
        "vcf" => Some(FileFormat::Vcf),
        "fai" => Some(FileFormat::Fai),
        "fa" | "fasta" | "fna" => Some(FileFormat::Fasta),
        "tsv" | "csv" => Some(FileFormat::Tsv),
        "txt" => {
            // Disambiguate .txt files based on filename patterns
            if lower_name.contains("assembly") || lower_name.contains("report") {
                Some(FileFormat::NcbiReport)
            } else if lower_name.ends_with(".dict.txt") {
                Some(FileFormat::Dict)
            } else {
                // Default to SAM for .txt files
                Some(FileFormat::Sam)
            }
        }
        _ => None,
    }
}

/// Detect format from file content analysis
fn detect_format_from_content(content: &str) -> Result<FileFormat, FormatError> {
    let content_trimmed = content.trim();

    // Check for empty content
    if content_trimmed.is_empty() {
        return Err(FormatError::UnknownFormat);
    }

    // Check for binary content (non-UTF8 or control characters)
    if content
        .chars()
        .any(|c| c.is_control() && c != '\n' && c != '\r' && c != '\t')
    {
        return Err(FormatError::UnsupportedBinary);
    }

    let lines: Vec<&str> = content_trimmed.lines().take(20).collect(); // Sample first 20 lines

    // Picard dictionary: starts with @HD and has @SQ lines (check BEFORE Sam)
    if lines.iter().any(|line| line.starts_with("@HD\t"))
        && lines.iter().any(|line| line.starts_with("@SQ\t"))
    {
        return Ok(FileFormat::Dict);
    }

    // SAM header format: starts with @SQ lines
    if lines.iter().any(|line| line.starts_with("@SQ\t")) {
        return Ok(FileFormat::Sam);
    }

    // VCF format: starts with ## comments and has ##contig lines
    if lines
        .iter()
        .any(|line| line.starts_with("##fileformat=VCF"))
        || (lines.iter().any(|line| line.starts_with("##"))
            && lines.iter().any(|line| line.starts_with("##contig=")))
    {
        return Ok(FileFormat::Vcf);
    }

    // NCBI assembly report: has specific column headers (check all lines for header)
    if lines.iter().any(|line| {
        line.contains("Sequence-Name")
            && line.contains("Sequence-Role")
            && line.contains("Assigned-Molecule")
    }) {
        return Ok(FileFormat::NcbiReport);
    }

    // TSV format: tab-separated with consistent column count
    if lines.len() > 1 {
        let first_line_cols = lines[0].split('\t').count();
        if first_line_cols > 2
            && lines
                .iter()
                .take(5)
                .all(|line| line.split('\t').count() == first_line_cols)
        {
            // Check if it looks like sequence data (has length/size columns)
            if lines[0].to_lowercase().contains("length")
                || lines[0].to_lowercase().contains("size")
                || lines[0].to_lowercase().contains("sequence")
            {
                return Ok(FileFormat::Tsv);
            }
        }
    }

    // CSV format: comma-separated
    if lines.len() > 1 {
        let first_line_cols = lines[0].split(',').count();
        if first_line_cols > 2
            && lines
                .iter()
                .take(5)
                .all(|line| line.split(',').count() == first_line_cols)
            && (lines[0].to_lowercase().contains("length")
                || lines[0].to_lowercase().contains("size")
                || lines[0].to_lowercase().contains("sequence"))
        {
            return Ok(FileFormat::Tsv);
        }
    }

    // FAI format: exactly 5 tab-separated columns (name, length, offset, line_bases, line_width)
    // All non-empty, non-comment lines should have 5 columns with numeric values in columns 2-5
    if !lines.is_empty() {
        let fai_lines: Vec<&&str> = lines
            .iter()
            .filter(|line| !line.is_empty() && !line.starts_with('#'))
            .collect();

        if !fai_lines.is_empty()
            && fai_lines.iter().all(|line| {
                let fields: Vec<&str> = line.split('\t').collect();
                if fields.len() != 5 {
                    return false;
                }
                // All fields after the first should be numeric
                fields[1..].iter().all(|f| f.parse::<u64>().is_ok())
            })
        {
            return Ok(FileFormat::Fai);
        }
    }

    // If content looks like plain text with sequence-like data, assume SAM header
    if lines.iter().any(|line| {
        line.contains("chr")
            || line.contains("scaffold")
            || line.contains("contig")
            || line.to_lowercase().contains("sequence")
            || line.to_lowercase().contains("length")
    }) {
        return Ok(FileFormat::Sam);
    }

    Err(FormatError::UnknownFormat)
}

/// Validate that content matches the expected format
#[allow(clippy::trivially_copy_pass_by_ref)] // Clearer API with reference
fn validate_format_content(content: &str, format: &FileFormat) -> bool {
    match format {
        FileFormat::Sam => {
            content.contains("@SQ") || content.contains("SN:") || content.contains("LN:")
        }
        FileFormat::Dict => content.contains("@HD") && content.contains("@SQ"),
        FileFormat::Vcf => {
            content.contains("##")
                && (content.contains("##contig=") || content.contains("##fileformat=VCF"))
        }
        FileFormat::NcbiReport => {
            content.contains("Sequence-Name") || content.contains("Sequence-Role")
        }
        FileFormat::Tsv => {
            content.contains('\t')
                && (content.to_lowercase().contains("length")
                    || content.to_lowercase().contains("sequence"))
        }
        FileFormat::Fai => {
            // FAI format has 5 tab-separated columns
            let lines: Vec<&str> = content.lines().take(5).collect();
            lines.iter().any(|line| {
                let fields: Vec<&str> = line.split('\t').collect();
                fields.len() == 5 && fields[1..].iter().all(|f| f.parse::<u64>().is_ok())
            })
        }
        FileFormat::Bam | FileFormat::Cram | FileFormat::Fasta => {
            // Binary formats should not be validated against text content
            false
        }
        FileFormat::Auto => true, // Auto-detect always passes validation
    }
}

/// Parse content with the specified format
///
/// # Errors
///
/// Returns `ParseError::ParseFailed` if the content cannot be parsed with the
/// specified format.
pub fn parse_with_format(content: &str, format: FileFormat) -> Result<QueryHeader, ParseError> {
    match format {
        FileFormat::Sam => {
            crate::parsing::sam::parse_header_text(content).map_err(|e| ParseError::ParseFailed {
                format: FileFormat::Sam,
                message: e.to_string(),
            })
        }
        FileFormat::Dict => {
            crate::parsing::dict::parse_dict_text(content).map_err(|e| ParseError::ParseFailed {
                format: FileFormat::Dict,
                message: e.to_string(),
            })
        }
        FileFormat::Vcf => crate::parsing::vcf::parse_vcf_header_text(content).map_err(|e| {
            ParseError::ParseFailed {
                format: FileFormat::Vcf,
                message: e.to_string(),
            }
        }),
        FileFormat::NcbiReport => {
            // NCBI report parser returns Vec<NcbiContigEntry>, we need to convert
            match crate::parsing::ncbi_report::parse_ncbi_report_text(content) {
                Ok(entries) => {
                    let contigs = entries.into_iter().map(|entry| entry.to_contig()).collect();
                    Ok(crate::core::header::QueryHeader::new(contigs))
                }
                Err(e) => Err(ParseError::ParseFailed {
                    format: FileFormat::NcbiReport,
                    message: e.to_string(),
                }),
            }
        }
        FileFormat::Tsv => {
            // TSV parser requires delimiter - try tab first, then comma
            match crate::parsing::tsv::parse_tsv_text(content, '\t') {
                Ok(query) => Ok(query),
                Err(_) => crate::parsing::tsv::parse_tsv_text(content, ',').map_err(|e| {
                    ParseError::ParseFailed {
                        format: FileFormat::Tsv,
                        message: format!("Failed to parse as TSV or CSV: {e}"),
                    }
                }),
            }
        }
        FileFormat::Fai => {
            crate::parsing::fai::parse_fai_text(content).map_err(|e| ParseError::ParseFailed {
                format: FileFormat::Fai,
                message: e.to_string(),
            })
        }
        FileFormat::Bam => Err(ParseError::ParseFailed {
            format: FileFormat::Bam,
            message: "BAM files must be parsed as binary, not text".to_string(),
        }),
        FileFormat::Cram => Err(ParseError::ParseFailed {
            format: FileFormat::Cram,
            message: "CRAM files must be parsed as binary, not text".to_string(),
        }),
        FileFormat::Fasta => Err(ParseError::ParseFailed {
            format: FileFormat::Fasta,
            message: "FASTA files must be parsed as binary, not text".to_string(),
        }),
        FileFormat::Auto => {
            // For auto-detection, detect format first then parse
            let detected_format =
                detect_format_from_content(content).map_err(|e| ParseError::ParseFailed {
                    format: FileFormat::Auto,
                    message: format!("Auto-detection failed: {e}"),
                })?;

            parse_with_format(content, detected_format)
        }
    }
}

/// Parse binary file content (for BAM/CRAM files)
///
/// Note: This creates a secure temporary file since the underlying parsers expect file paths
///
/// # Errors
///
/// Returns `ParseError::Io` if the temporary file cannot be created or written,
/// or `ParseError::ParseFailed` if parsing fails.
pub fn parse_binary_file(
    file_content: &[u8],
    format: FileFormat,
) -> Result<QueryHeader, ParseError> {
    use std::io::Write;
    use tempfile::NamedTempFile;

    match format {
        FileFormat::Bam | FileFormat::Cram => {
            // Create a secure temporary file with cryptographically random name
            let file_extension = match format {
                FileFormat::Bam => ".bam",
                FileFormat::Cram => ".cram",
                _ => ".bin",
            };

            let mut temp_file =
                NamedTempFile::with_suffix(file_extension).map_err(ParseError::Io)?;

            // Write bytes to secure temporary file
            temp_file.write_all(file_content).map_err(ParseError::Io)?;

            // Parse the file using the secure temp file path
            let result = crate::parsing::sam::parse_file(temp_file.path());

            // File automatically deleted when NamedTempFile drops
            result.map_err(|_e| ParseError::ParseFailed {
                format,
                message: "Binary file parsing failed".to_string(), // Sanitized error message
            })
        }
        FileFormat::Fasta => {
            // Determine appropriate extension based on content (check for gzip magic bytes)
            let is_gzipped =
                file_content.len() >= 2 && file_content[0] == 0x1f && file_content[1] == 0x8b;
            let file_extension = if is_gzipped { ".fa.gz" } else { ".fa" };

            let mut temp_file =
                NamedTempFile::with_suffix(file_extension).map_err(ParseError::Io)?;

            // Write bytes to secure temporary file
            temp_file.write_all(file_content).map_err(ParseError::Io)?;

            // Parse the FASTA file using the secure temp file path
            let result = crate::parsing::fasta::parse_fasta_file(temp_file.path());

            // File automatically deleted when NamedTempFile drops
            result.map_err(|_e| ParseError::ParseFailed {
                format,
                message: "Binary file parsing failed".to_string(), // Sanitized error message
            })
        }
        _ => Err(ParseError::ParseFailed {
            format,
            message: "Format is not a binary file format".to_string(),
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_filename_detection() {
        assert_eq!(
            detect_format_from_filename("test.sam"),
            Some(FileFormat::Sam)
        );
        assert_eq!(
            detect_format_from_filename("test.bam"),
            Some(FileFormat::Bam)
        );
        assert_eq!(
            detect_format_from_filename("test.dict"),
            Some(FileFormat::Dict)
        );
        assert_eq!(
            detect_format_from_filename("test.vcf"),
            Some(FileFormat::Vcf)
        );
        assert_eq!(
            detect_format_from_filename("test.vcf.gz"),
            Some(FileFormat::Vcf)
        );
        assert_eq!(
            detect_format_from_filename("assembly_report.txt"),
            Some(FileFormat::NcbiReport)
        );
        assert_eq!(
            detect_format_from_filename("reference.fai"),
            Some(FileFormat::Fai)
        );
        assert_eq!(
            detect_format_from_filename("reference.fa"),
            Some(FileFormat::Fasta)
        );
        assert_eq!(
            detect_format_from_filename("reference.fasta"),
            Some(FileFormat::Fasta)
        );
        assert_eq!(
            detect_format_from_filename("reference.fa.gz"),
            Some(FileFormat::Fasta)
        );
        assert_eq!(
            detect_format_from_filename("reference.fasta.gz"),
            Some(FileFormat::Fasta)
        );
        assert_eq!(detect_format_from_filename("unknown.xyz"), None);
    }

    #[test]
    fn test_sam_header_detection() {
        let content = "@SQ\tSN:chr1\tLN:248956422\tM5:6aef897c3d6ff0c78aff06ac189178dd\n";
        assert_eq!(detect_format_from_content(content), Ok(FileFormat::Sam));
    }

    #[test]
    fn test_dict_detection() {
        let content = "@HD\tVN:1.0\tSO:coordinate\n@SQ\tSN:chr1\tLN:248956422\tM5:abc123\n";
        assert_eq!(detect_format_from_content(content), Ok(FileFormat::Dict));
    }

    #[test]
    fn test_vcf_detection() {
        let content = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=248956422>\n";
        assert_eq!(detect_format_from_content(content), Ok(FileFormat::Vcf));
    }

    #[test]
    fn test_ncbi_report_detection() {
        let content =
            "# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\n";
        assert_eq!(
            detect_format_from_content(content),
            Ok(FileFormat::NcbiReport)
        );
    }

    #[test]
    fn test_fai_detection() {
        let content = "chr1\t248956422\t112\t70\t71\nchr2\t242193529\t253404903\t70\t71\n";
        assert_eq!(detect_format_from_content(content), Ok(FileFormat::Fai));
    }

    #[test]
    fn test_fai_validation() {
        assert!(validate_format_content(
            "chr1\t248956422\t112\t70\t71",
            &FileFormat::Fai
        ));
        assert!(!validate_format_content(
            "chr1\t248956422\t112",
            &FileFormat::Fai
        ));
    }

    #[test]
    fn test_format_validation() {
        assert!(validate_format_content(
            "@SQ\tSN:chr1\tLN:123",
            &FileFormat::Sam
        ));
        assert!(!validate_format_content("random text", &FileFormat::Sam));

        assert!(validate_format_content(
            "##contig=<ID=chr1>",
            &FileFormat::Vcf
        ));
        assert!(!validate_format_content("@SQ\tSN:chr1", &FileFormat::Vcf));
    }

    #[test]
    fn test_combined_detection() {
        let content = "@SQ\tSN:chr1\tLN:248956422\n";
        assert_eq!(
            detect_format(content, Some("test.sam")),
            Ok(FileFormat::Sam)
        );
        assert_eq!(
            detect_format(content, Some("test.dict")),
            Ok(FileFormat::Sam)
        ); // Content overrides filename
        assert_eq!(detect_format(content, None), Ok(FileFormat::Sam));
    }
}
