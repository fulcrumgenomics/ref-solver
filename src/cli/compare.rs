use std::path::{Path, PathBuf};

use clap::Args;

use crate::catalog::store::ReferenceCatalog;
use crate::cli::OutputFormat;
use crate::core::header::QueryHeader;
use crate::core::reference::KnownReference;
use crate::core::types::{Assembly, ReferenceSource};
use crate::matching::scoring::MatchScore;
use crate::parsing;

#[derive(Args)]
pub struct CompareArgs {
    /// First input file (BAM, SAM, CRAM, .dict, or TSV)
    #[arg(required = true)]
    pub input_a: PathBuf,

    /// Second input file, or reference ID from catalog
    #[arg(required = true)]
    pub input_b: String,

    /// Treat second argument as a reference ID from the catalog
    #[arg(long)]
    pub reference: bool,

    /// Path to custom catalog file (only used with --reference)
    #[arg(long)]
    pub catalog: Option<PathBuf>,
}

pub fn run(args: CompareArgs, format: OutputFormat, verbose: bool) -> anyhow::Result<()> {
    // Parse first input
    let query_a = parse_input(&args.input_a)?;

    if verbose {
        eprintln!(
            "Input A: {} contigs ({}% have MD5)",
            query_a.contigs.len(),
            (query_a.md5_coverage() * 100.0) as u32
        );
    }

    // Parse second input
    let query_b = if args.reference {
        // Look up in catalog
        let catalog = if let Some(path) = &args.catalog {
            ReferenceCatalog::load_from_file(path)?
        } else {
            ReferenceCatalog::load_embedded()?
        };

        let ref_id = crate::core::types::ReferenceId::new(&args.input_b);
        let reference = catalog
            .get(&ref_id)
            .ok_or_else(|| anyhow::anyhow!("Reference '{}' not found in catalog", args.input_b))?;

        // Convert reference to query header for comparison
        QueryHeader::new(reference.contigs.clone())
    } else {
        // Parse as file
        let path = PathBuf::from(&args.input_b);
        parse_input(&path)?
    };

    if verbose {
        eprintln!(
            "Input B: {} contigs ({}% have MD5)",
            query_b.contigs.len(),
            (query_b.md5_coverage() * 100.0) as u32
        );
    }

    // Compare - use builder pattern to ensure indexes are properly computed
    let reference_b = KnownReference::new(
        "input_b",
        &args.input_b,
        Assembly::Other("unknown".to_string()),
        ReferenceSource::Custom("input".to_string()),
    )
    .with_contigs(query_b.contigs.clone());

    let score = MatchScore::calculate(&query_a, &reference_b);

    // Output results
    match format {
        OutputFormat::Text => print_text_comparison(&args, &query_a, &query_b, &score),
        OutputFormat::Json => print_json_comparison(&args, &query_a, &query_b, &score)?,
        OutputFormat::Tsv => print_tsv_comparison(&score),
    }

    Ok(())
}

fn parse_input(path: &Path) -> anyhow::Result<QueryHeader> {
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .map(|s| s.to_lowercase());

    match ext.as_deref() {
        Some("bam" | "sam" | "cram") => Ok(parsing::sam::parse_file(path)?),
        Some("dict") => Ok(parsing::dict::parse_dict_file(path)?),
        Some("tsv") => Ok(parsing::tsv::parse_tsv_file(path, '\t')?),
        Some("csv") => Ok(parsing::tsv::parse_tsv_file(path, ',')?),
        _ => Ok(parsing::sam::parse_file(path)?),
    }
}

fn print_text_comparison(
    args: &CompareArgs,
    query_a: &QueryHeader,
    query_b: &QueryHeader,
    score: &MatchScore,
) {
    println!("Comparison Results");
    println!("{}", "=".repeat(60));

    println!("\nInput A: {}", args.input_a.display());
    println!("  Contigs: {}", query_a.contigs.len());
    println!("  MD5 coverage: {:.0}%", query_a.md5_coverage() * 100.0);
    println!("  Naming convention: {:?}", query_a.naming_convention);

    println!("\nInput B: {}", args.input_b);
    println!("  Contigs: {}", query_b.contigs.len());
    println!("  MD5 coverage: {:.0}%", query_b.md5_coverage() * 100.0);
    println!("  Naming convention: {:?}", query_b.naming_convention);

    println!("\nSimilarity Scores:");
    println!("  MD5 Jaccard: {:.2}%", score.md5_jaccard * 100.0);
    println!(
        "  Name+Length Jaccard: {:.2}%",
        score.name_length_jaccard * 100.0
    );
    println!(
        "  MD5 Query Coverage: {:.2}%",
        score.md5_query_coverage * 100.0
    );
    println!(
        "  Name+Length Query Coverage: {:.2}%",
        score.name_length_query_coverage * 100.0
    );
    println!("  Order Preserved: {}", score.order_preserved);
    println!("  Order Score: {:.2}%", score.order_score * 100.0);
    println!("  Composite Score: {:.2}%", score.composite * 100.0);
    println!("  Confidence: {:?}", score.confidence);
}

fn print_json_comparison(
    args: &CompareArgs,
    query_a: &QueryHeader,
    query_b: &QueryHeader,
    score: &MatchScore,
) -> anyhow::Result<()> {
    let output = serde_json::json!({
        "input_a": {
            "path": args.input_a.display().to_string(),
            "contig_count": query_a.contigs.len(),
            "md5_coverage": query_a.md5_coverage(),
            "naming_convention": format!("{:?}", query_a.naming_convention),
        },
        "input_b": {
            "path": args.input_b,
            "contig_count": query_b.contigs.len(),
            "md5_coverage": query_b.md5_coverage(),
            "naming_convention": format!("{:?}", query_b.naming_convention),
        },
        "score": {
            "md5_jaccard": score.md5_jaccard,
            "name_length_jaccard": score.name_length_jaccard,
            "md5_query_coverage": score.md5_query_coverage,
            "name_length_query_coverage": score.name_length_query_coverage,
            "order_preserved": score.order_preserved,
            "order_score": score.order_score,
            "composite": score.composite,
            "confidence": format!("{:?}", score.confidence),
        }
    });

    println!("{}", serde_json::to_string_pretty(&output)?);
    Ok(())
}

fn print_tsv_comparison(score: &MatchScore) {
    println!(
        "md5_jaccard\tname_length_jaccard\tmd5_query_coverage\torder_preserved\tcomposite\tconfidence"
    );
    println!(
        "{:.4}\t{:.4}\t{:.4}\t{}\t{:.4}\t{:?}",
        score.md5_jaccard,
        score.name_length_jaccard,
        score.md5_query_coverage,
        score.order_preserved,
        score.composite,
        score.confidence,
    );
}
