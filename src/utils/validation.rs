//! Centralized validation and helper functions.

use crate::web::format_detection::FileFormat;
use std::collections::HashSet;

/// Maximum number of contigs allowed in a single file (DOS protection)
pub const MAX_CONTIGS: usize = 100_000;

/// Security-related constants for input validation
pub const MAX_FILENAME_LENGTH: usize = 255;
pub const MIN_FILE_CONTENT_SIZE: usize = 1;

/// Validate that a string is a valid MD5 checksum (32 hex characters).
///
/// # Examples
///
/// ```
/// use ref_solver::utils::validation::is_valid_md5;
///
/// assert!(is_valid_md5("6aef897c3d6ff0c78aff06ac189178dd"));
/// assert!(!is_valid_md5("not-an-md5"));
/// assert!(!is_valid_md5("6aef897c3d6ff0c78aff06ac189178d")); // 31 chars
/// ```
#[must_use]
pub fn is_valid_md5(s: &str) -> bool {
    s.len() == 32 && s.chars().all(|c| c.is_ascii_hexdigit())
}

/// Normalize an MD5 string to lowercase.
/// Returns None if the input is not a valid MD5.
#[must_use]
pub fn normalize_md5(s: &str) -> Option<String> {
    if is_valid_md5(s) {
        Some(s.to_lowercase())
    } else {
        None
    }
}

/// Compute a signature hash from a set of MD5 checksums.
///
/// The signature is computed by:
/// 1. Sorting the MD5s alphabetically
/// 2. Joining them with commas
/// 3. Computing MD5 of the concatenated string
///
/// This provides a deterministic identifier for a set of contigs.
#[must_use]
#[allow(clippy::implicit_hasher)] // Default hasher is fine for this use case
pub fn compute_signature(md5s: &HashSet<String>) -> String {
    if md5s.is_empty() {
        return String::new();
    }

    let mut sorted: Vec<&str> = md5s.iter().map(std::string::String::as_str).collect();
    sorted.sort_unstable();
    let concatenated = sorted.join(",");
    let digest = md5::compute(concatenated.as_bytes());
    format!("{digest:x}")
}

/// Check if adding another contig would exceed the maximum allowed.
///
/// Call this with the current count BEFORE adding a new contig.
/// Returns an error message if adding would exceed the limit, None if safe to add.
///
/// # Example
/// ```ignore
/// if check_contig_limit(contigs.len()).is_some() {
///     return Err(...);
/// }
/// contigs.push(new_contig); // Safe to add
/// ```
#[must_use]
pub fn check_contig_limit(count: usize) -> Option<String> {
    if count >= MAX_CONTIGS {
        Some(format!(
            "Too many contigs: adding another would exceed maximum of {MAX_CONTIGS}"
        ))
    } else {
        None
    }
}

/// Security validation error types
#[derive(Debug, thiserror::Error)]
pub enum ValidationError {
    #[error("Filename too long: exceeds {MAX_FILENAME_LENGTH} characters")]
    FilenameTooLong,
    #[error("Invalid filename: contains path traversal or invalid characters")]
    InvalidFilename,
    #[error("Empty filename provided")]
    EmptyFilename,
    #[error("File content appears malformed or invalid")]
    InvalidFileContent,
    #[error("File format validation failed")]
    FormatValidationFailed,
}

/// Secure filename validation to prevent directory traversal and other attacks
///
/// Validates and sanitizes filenames by:
/// - Checking length limits
/// - Preventing directory traversal (../, ..\\)
/// - Removing potentially dangerous characters
/// - Ensuring filename is not empty after sanitization
///
/// # Errors
///
/// Returns `ValidationError::EmptyFilename` if the filename is empty,
/// `ValidationError::FilenameTooLong` if it exceeds the limit, or
/// `ValidationError::InvalidFilename` if it contains invalid characters.
pub fn validate_filename(filename: &str) -> Result<String, ValidationError> {
    // Check if filename is empty
    if filename.trim().is_empty() {
        return Err(ValidationError::EmptyFilename);
    }

    // Check length limit
    if filename.len() > MAX_FILENAME_LENGTH {
        return Err(ValidationError::FilenameTooLong);
    }

    // Prevent directory traversal attacks
    if filename.contains("..") || filename.contains('/') || filename.contains('\\') {
        return Err(ValidationError::InvalidFilename);
    }

    // Check for null bytes and other dangerous characters
    if filename.contains('\0') || filename.chars().any(|c| ('\x01'..='\x1F').contains(&c)) {
        return Err(ValidationError::InvalidFilename);
    }

    // Sanitize filename by keeping only safe characters
    let sanitized = filename
        .chars()
        .filter(|c| c.is_ascii_alphanumeric() || *c == '.' || *c == '-' || *c == '_' || *c == ' ')
        .collect::<String>();

    // Ensure sanitized filename is not empty
    if sanitized.trim().is_empty() {
        return Err(ValidationError::InvalidFilename);
    }

    // Prevent hidden files (starting with .) unless it's a known extension
    if sanitized.starts_with('.') && !has_known_extension(&sanitized) {
        return Err(ValidationError::InvalidFilename);
    }

    Ok(sanitized)
}

/// Check if filename has a known safe extension
fn has_known_extension(filename: &str) -> bool {
    let safe_extensions = [
        ".sam",
        ".bam",
        ".cram",
        ".dict",
        ".vcf",
        ".txt",
        ".tsv",
        ".csv",
        ".gz",
        ".assembly_report.txt",
    ];

    safe_extensions
        .iter()
        .any(|ext| filename.to_lowercase().ends_with(ext))
}

/// Validate file content using magic numbers for known binary formats
///
/// Performs format validation by checking magic numbers (file signatures)
/// to prevent format confusion attacks and ensure file integrity
#[must_use]
pub fn validate_file_format(content: &[u8], expected_format: FileFormat) -> bool {
    if content.is_empty() {
        return false;
    }

    match expected_format {
        FileFormat::Bam => {
            // BAM files start with "BAM\x01"
            content.len() >= 4 && content.starts_with(b"BAM\x01")
        }
        FileFormat::Cram => {
            // CRAM files start with "CRAM"
            content.len() >= 4 && content.starts_with(b"CRAM")
        }
        FileFormat::Vcf => {
            // VCF files should start with "##fileformat=VCF"
            let content_str = std::str::from_utf8(content).unwrap_or("");
            content_str.starts_with("##fileformat=VCF")
        }
        FileFormat::Sam => {
            // SAM files are text-based, check for header indicators
            let content_str = std::str::from_utf8(content).unwrap_or("");
            content_str.contains("@SQ")
                || content_str.contains("@HD")
                || content_str.contains("SN:")
                || content_str.contains("LN:")
        }
        FileFormat::Dict => {
            // Picard dictionary files have @HD and @SQ headers
            let content_str = std::str::from_utf8(content).unwrap_or("");
            content_str.contains("@HD") && content_str.contains("@SQ")
        }
        FileFormat::NcbiReport => {
            // NCBI assembly reports have specific column headers
            let content_str = std::str::from_utf8(content).unwrap_or("");
            content_str.contains("Sequence-Name") || content_str.contains("Sequence-Role")
        }
        FileFormat::Tsv => {
            // TSV files should have tab-separated content
            let content_str = std::str::from_utf8(content).unwrap_or("");
            content_str.contains('\t')
                && (content_str.to_lowercase().contains("length")
                    || content_str.to_lowercase().contains("sequence")
                    || content_str.to_lowercase().contains("size"))
        }
        FileFormat::Auto => {
            // Auto-detection always passes initial validation
            true
        }
    }
}

/// Validate that file content is not malicious or malformed
///
/// Basic security checks for file content integrity:
/// - Minimum size requirements
/// - Binary content detection for text formats
/// - Basic malformation checks
///
/// # Errors
///
/// Returns `ValidationError::InvalidFileContent` if the content is too small,
/// contains unexpected binary data for text formats, or fails UTF-8 validation.
pub fn validate_file_content(content: &[u8], expected_text: bool) -> Result<(), ValidationError> {
    // Check minimum content size
    if content.len() < MIN_FILE_CONTENT_SIZE {
        return Err(ValidationError::InvalidFileContent);
    }

    // If we expect text content, validate it's not binary
    if expected_text {
        // Check for excessive non-printable characters
        let non_printable_count = content
            .iter()
            .filter(|&&b| b < 9 || (b > 13 && b < 32) || b > 126)
            .count();

        // Allow up to 5% non-printable characters for text files
        if content.len() > 100 && non_printable_count > content.len() / 20 {
            return Err(ValidationError::InvalidFileContent);
        }

        // Basic UTF-8 validation for text content
        if std::str::from_utf8(content).is_err() {
            return Err(ValidationError::InvalidFileContent);
        }
    }

    Ok(())
}

/// Comprehensive input validation combining filename and content checks
///
/// Performs complete security validation for file uploads:
/// - Filename sanitization and security checks
/// - File format validation via magic numbers
/// - Content integrity validation
///
/// # Errors
///
/// Returns a `ValidationError` if filename validation fails, the file format
/// doesn't match the expected format, or content validation fails.
pub fn validate_upload(
    filename: Option<&str>,
    content: &[u8],
    expected_format: FileFormat,
) -> Result<Option<String>, ValidationError> {
    // Validate filename if provided
    let validated_filename = if let Some(name) = filename {
        Some(validate_filename(name)?)
    } else {
        None
    };

    // Validate content integrity
    let is_text_format = matches!(
        expected_format,
        FileFormat::Sam
            | FileFormat::Dict
            | FileFormat::Vcf
            | FileFormat::NcbiReport
            | FileFormat::Tsv
            | FileFormat::Auto
    );

    validate_file_content(content, is_text_format)?;

    // Validate file format - even for auto-detection, check for obvious mismatches
    if expected_format == FileFormat::Auto {
        // For auto-detection, at least verify it's not a malformed binary file
        // Check if it looks like a known binary format but is malformed
        if content.len() >= 4 {
            let starts_with_bam = content.starts_with(b"BAM");
            let starts_with_cram = content.starts_with(b"CRAM");

            // If it looks like it should be BAM/CRAM but isn't valid, reject it
            if starts_with_bam && !validate_file_format(content, FileFormat::Bam) {
                return Err(ValidationError::FormatValidationFailed);
            }
            if starts_with_cram && !validate_file_format(content, FileFormat::Cram) {
                return Err(ValidationError::FormatValidationFailed);
            }
        }
    } else if !validate_file_format(content, expected_format) {
        return Err(ValidationError::FormatValidationFailed);
    }

    Ok(validated_filename)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_valid_md5() {
        assert!(is_valid_md5("6aef897c3d6ff0c78aff06ac189178dd"));
        assert!(is_valid_md5("AABBCCDD11223344556677889900AABB")); // uppercase ok
        assert!(!is_valid_md5("not-an-md5"));
        assert!(!is_valid_md5("6aef897c3d6ff0c78aff06ac189178d")); // 31 chars
        assert!(!is_valid_md5("6aef897c3d6ff0c78aff06ac189178ddd")); // 33 chars
        assert!(!is_valid_md5("")); // empty
        assert!(!is_valid_md5("6aef897c3d6ff0c78aff06ac189178dg")); // invalid char
    }

    #[test]
    fn test_normalize_md5() {
        assert_eq!(
            normalize_md5("6AEF897C3D6FF0C78AFF06AC189178DD"),
            Some("6aef897c3d6ff0c78aff06ac189178dd".to_string())
        );
        assert_eq!(normalize_md5("invalid"), None);
    }

    #[test]
    fn test_compute_signature() {
        let mut md5s = HashSet::new();
        md5s.insert("aaaa".repeat(8)); // fake MD5
        md5s.insert("bbbb".repeat(8));

        let sig = compute_signature(&md5s);
        assert_eq!(sig.len(), 32);

        // Same input should give same output
        let sig2 = compute_signature(&md5s);
        assert_eq!(sig, sig2);

        // Empty set gives empty string
        let empty: HashSet<String> = HashSet::new();
        assert_eq!(compute_signature(&empty), "");
    }

    #[test]
    fn test_check_contig_limit() {
        assert!(check_contig_limit(100).is_none());
        assert!(check_contig_limit(MAX_CONTIGS - 1).is_none());
        assert!(check_contig_limit(MAX_CONTIGS).is_some());
        assert!(check_contig_limit(MAX_CONTIGS + 1).is_some());
    }

    // Security validation tests
    #[test]
    fn test_validate_filename_safe() {
        assert!(validate_filename("test.sam").is_ok());
        assert!(validate_filename("my-file.bam").is_ok());
        assert!(validate_filename("data_file.txt").is_ok());
        assert!(validate_filename("sample 123.vcf").is_ok());
    }

    #[test]
    fn test_validate_filename_dangerous() {
        // Directory traversal attempts
        assert!(validate_filename("../etc/passwd").is_err());
        assert!(validate_filename("..\\windows\\system32").is_err());
        assert!(validate_filename("test/../../secret").is_err());

        // Null bytes and control characters
        assert!(validate_filename("test\0.txt").is_err());
        assert!(validate_filename("test\x01.txt").is_err());

        // Too long filename
        let long_name = "a".repeat(300);
        assert!(validate_filename(&long_name).is_err());

        // Empty or whitespace-only
        assert!(validate_filename("").is_err());
        assert!(validate_filename("   ").is_err());

        // Hidden files without known extensions
        assert!(validate_filename(".hidden").is_err());
    }

    #[test]
    fn test_validate_filename_sanitization() {
        // Should remove dangerous characters but keep safe ones
        let result = validate_filename("test@#$%file.txt").unwrap();
        assert_eq!(result, "testfile.txt");

        // Should preserve safe characters
        let result = validate_filename("my-file_123.sam").unwrap();
        assert_eq!(result, "my-file_123.sam");
    }

    #[test]
    fn test_validate_file_format_bam() {
        let bam_content = b"BAM\x01test_content";
        assert!(validate_file_format(bam_content, FileFormat::Bam));

        let invalid_bam = b"NOTBAM\x01";
        assert!(!validate_file_format(invalid_bam, FileFormat::Bam));
    }

    #[test]
    fn test_validate_file_format_cram() {
        let cram_content = b"CRAMtest_content";
        assert!(validate_file_format(cram_content, FileFormat::Cram));

        let invalid_cram = b"NOTCRAM";
        assert!(!validate_file_format(invalid_cram, FileFormat::Cram));
    }

    #[test]
    fn test_validate_file_format_vcf() {
        let vcf_content = b"##fileformat=VCFv4.2\n##contig=<ID=chr1>";
        assert!(validate_file_format(vcf_content, FileFormat::Vcf));

        let invalid_vcf = b"@SQ\tSN:chr1\tLN:123";
        assert!(!validate_file_format(invalid_vcf, FileFormat::Vcf));
    }

    #[test]
    fn test_validate_file_format_sam() {
        let sam_content = b"@SQ\tSN:chr1\tLN:123456";
        assert!(validate_file_format(sam_content, FileFormat::Sam));

        let sam_content2 = b"@HD\tVN:1.0\tSO:coordinate";
        assert!(validate_file_format(sam_content2, FileFormat::Sam));
    }

    #[test]
    fn test_validate_file_content_text() {
        let valid_text = b"@SQ\tSN:chr1\tLN:123456\n@SQ\tSN:chr2\tLN:654321";
        assert!(validate_file_content(valid_text, true).is_ok());

        // Too much binary data for text format
        let binary_data = vec![0u8; 1000];
        assert!(validate_file_content(&binary_data, true).is_err());

        // Empty content
        assert!(validate_file_content(b"", true).is_err());
    }

    #[test]
    fn test_validate_file_content_binary() {
        let binary_data = vec![0xABu8; 100];
        assert!(validate_file_content(&binary_data, false).is_ok());

        // Empty content still invalid for binary
        assert!(validate_file_content(b"", false).is_err());
    }

    #[test]
    fn test_validate_upload_complete() {
        let sam_content = b"@SQ\tSN:chr1\tLN:123456";

        // Valid upload with filename
        let result = validate_upload(Some("test.sam"), sam_content, FileFormat::Sam);
        assert!(result.is_ok());
        assert_eq!(result.unwrap().unwrap(), "test.sam");

        // Valid upload without filename
        let result = validate_upload(None, sam_content, FileFormat::Sam);
        assert!(result.is_ok());
        assert!(result.unwrap().is_none());

        // Invalid filename
        let result = validate_upload(Some("../etc/passwd"), sam_content, FileFormat::Sam);
        assert!(result.is_err());

        // Format mismatch
        let bam_content = b"BAM\x01test";
        let result = validate_upload(Some("test.sam"), bam_content, FileFormat::Sam);
        assert!(result.is_err());
    }

    #[test]
    fn test_has_known_extension() {
        assert!(has_known_extension(".sam"));
        assert!(has_known_extension(".bam"));
        assert!(has_known_extension(".vcf.gz"));
        assert!(has_known_extension("test.assembly_report.txt"));

        assert!(!has_known_extension(".exe"));
        assert!(!has_known_extension(".hidden"));
        assert!(!has_known_extension(".config"));
    }
}
