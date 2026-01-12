//! Parsers for extracting sequence dictionaries from various file formats.
//!
//! This module provides parsers for:
//!
//! - **SAM/BAM/CRAM files**: Extract `@SQ` lines from alignment file headers
//! - **Picard .dict files**: Parse sequence dictionary files
//! - **FASTA index (.fai) files**: Parse FASTA index files
//! - **NCBI assembly reports**: Parse NCBI assembly reports with multiple naming conventions
//! - **VCF headers**: Extract `##contig` lines from VCF files
//! - **TSV/CSV files**: Parse tabular contig definitions
//!
//! ## Example
//!
//! ```rust,no_run
//! use ref_solver::parsing::sam::{parse_file, parse_header_text};
//! use std::path::Path;
//!
//! // Parse from a BAM file
//! let query = parse_file(Path::new("sample.bam")).unwrap();
//!
//! // Or parse from raw header text
//! let header = "@SQ\tSN:chr1\tLN:248956422\tM5:6aef897c3d6ff0c78aff06ac189178dd\n";
//! let query = parse_header_text(header).unwrap();
//! ```
//!
//! ## Supported Tags
//!
//! From SAM `@SQ` lines, the following tags are extracted:
//!
//! | Tag | Description | Required |
//! |-----|-------------|----------|
//! | SN  | Sequence name | Yes |
//! | LN  | Sequence length | Yes |
//! | M5  | MD5 checksum | No |
//! | AS  | Assembly identifier | No |
//! | UR  | URI for sequence | No |
//! | SP  | Species | No |
//! | AN  | Alternate names (aliases) | No |

pub mod dict;
pub mod fai;
pub mod fasta;
pub mod ncbi_report;
pub mod sam;
pub mod tsv;
pub mod vcf;
