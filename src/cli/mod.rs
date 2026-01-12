//! Command-line interface for ref-solver.
//!
//! This module implements the CLI using clap. Available commands:
//!
//! - **identify**: Identify the reference genome from a BAM/SAM/CRAM file
//! - **compare**: Compare two headers or a header against a known reference
//! - **catalog**: List, show, or export references from the catalog
//! - **serve**: Start the interactive web interface
//!
//! ## Usage
//!
//! ```text
//! # Identify reference from a BAM file
//! ref-solver identify sample.bam
//!
//! # Pipe from samtools
//! samtools view -H sample.bam | ref-solver identify -
//!
//! # JSON output for scripting
//! ref-solver identify sample.bam --format json
//!
//! # Compare against a known reference
//! ref-solver compare sample.bam hg38_ucsc --reference
//!
//! # Start web UI
//! ref-solver serve --port 8080 --open
//! ```

use clap::{Parser, Subcommand};

pub mod catalog;
pub mod compare;
pub mod identify;

#[derive(Parser)]
#[command(name = "ref-solver")]
#[command(author = "Fulcrum Genomics")]
#[command(version)]
#[command(about = "Identify and match reference genomes from BAM/SAM headers")]
#[command(
    long_about = "ref-solver helps you identify which reference genome was used to align a BAM/SAM/CRAM file.\n\nIt matches the sequence dictionary from your file against a catalog of known human reference genomes and provides:\n- Exact matches when possible\n- Detailed diagnostics when differences exist\n- Actionable suggestions for fixing mismatches (renaming, reordering)"
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,

    /// Enable verbose output
    #[arg(short, long, global = true)]
    pub verbose: bool,

    /// Output format
    #[arg(short, long, global = true, default_value = "text")]
    pub format: OutputFormat,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Identify the reference genome used in a BAM/SAM file
    Identify(identify::IdentifyArgs),

    /// Compare two headers or references
    Compare(compare::CompareArgs),

    /// Manage the reference catalog
    Catalog(catalog::CatalogArgs),

    /// Start the web server
    Serve(ServeArgs),
}

#[derive(clap::Args)]
pub struct ServeArgs {
    /// Port to listen on
    #[arg(short, long, default_value = "8080")]
    pub port: u16,

    /// Address to bind to
    #[arg(short, long, default_value = "127.0.0.1")]
    pub address: String,

    /// Open browser automatically
    #[arg(long)]
    pub open: bool,
}

#[derive(Clone, Copy, Debug, clap::ValueEnum)]
pub enum OutputFormat {
    Text,
    Json,
    Tsv,
}
