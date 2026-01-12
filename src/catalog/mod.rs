//! Reference genome catalog storage and indexing.
//!
//! The catalog contains definitions of known human reference genomes with their
//! contigs, MD5 checksums, and metadata. An embedded catalog is compiled into
//! the binary, but custom catalogs can also be loaded from JSON files.
//!
//! ## Embedded Catalog
//!
//! The default catalog includes common reference genomes:
//!
//! - **GRCh38/hg38**: UCSC, NCBI, Broad analysis set, DRAGEN, hs38, hs38DH
//! - **GRCh37/hg19**: UCSC, NCBI, Broad b37, hs37d5
//! - **T2T-CHM13**: Complete telomere-to-telomere assembly
//!
//! ## Example
//!
//! ```rust,no_run
//! use ref_solver::ReferenceCatalog;
//! use ref_solver::core::types::ReferenceId;
//!
//! // Load embedded catalog
//! let catalog = ReferenceCatalog::load_embedded().unwrap();
//!
//! // List all references
//! for reference in &catalog.references {
//!     println!("{}", reference.id);
//! }
//!
//! // Get a specific reference
//! let hg38 = catalog.get(&ReferenceId::new("hg38_ucsc"));
//! ```
//!
//! ## Custom Catalogs
//!
//! Custom catalogs can be created by exporting and modifying the embedded catalog:
//!
//! ```rust,no_run
//! use ref_solver::ReferenceCatalog;
//! use std::path::Path;
//!
//! // Export to JSON
//! let catalog = ReferenceCatalog::load_embedded().unwrap();
//! let json = catalog.to_json().unwrap();
//!
//! // Load from custom file
//! let custom = ReferenceCatalog::load_from_file(Path::new("my_catalog.json")).unwrap();
//! ```

pub mod builder;
pub mod hierarchical;
pub mod index;
pub mod store;
