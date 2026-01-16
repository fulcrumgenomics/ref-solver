use std::path::Path;

use crate::core::header::QueryHeader;
use crate::parsing::sam::ParseError;

/// Parse a Picard sequence dictionary (.dict) file
///
/// # Errors
///
/// Returns `ParseError::Io` if the file cannot be read, or other parse errors
/// if the content is invalid.
pub fn parse_dict_file(path: &Path) -> Result<QueryHeader, ParseError> {
    let content = std::fs::read_to_string(path)?;
    parse_dict_text(&content)
}

/// Parse dictionary from text
///
/// # Errors
///
/// Returns `ParseError::InvalidFormat` if the text is not valid dictionary format,
/// or `ParseError::TooManyContigs` if the number of contigs exceeds the maximum.
pub fn parse_dict_text(text: &str) -> Result<QueryHeader, ParseError> {
    // .dict files are essentially SAM headers with only @HD and @SQ lines
    crate::parsing::sam::parse_header_text(text)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_dict_text() {
        let dict = r"@HD	VN:1.6
@SQ	SN:chr1	LN:248_956_422	M5:6aef897c3d6ff0c78aff06ac189178dd	UR:file:///reference/hg38.fa
@SQ	SN:chr2	LN:242_193_529	M5:f98db672eb0993dcfdabafe2a882905c	UR:file:///reference/hg38.fa
";

        let query = parse_dict_text(dict).unwrap();
        assert_eq!(query.contigs.len(), 2);
        assert_eq!(query.contigs[0].name, "chr1");
        assert_eq!(
            query.contigs[0].uri,
            Some("file:///reference/hg38.fa".to_string())
        );
    }
}
