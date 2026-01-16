//! Score command - compare two files directly using the scoring algorithm.
//!
//! This command compares a query file against a reference file without using
//! the reference catalog. Useful for comparing arbitrary files.

use std::path::{Path, PathBuf};

use clap::Args;

use crate::cli::OutputFormat;
use crate::core::header::QueryHeader;
use crate::core::reference::KnownReference;
use crate::core::types::{Assembly, ReferenceSource};
use crate::matching::engine::ScoringWeights;
use crate::matching::scoring::MatchScore;
use crate::parsing;

/// Arguments for the score command
#[derive(Args)]
pub struct ScoreArgs {
    /// Query file (the file you want to score)
    /// Supports: BAM, SAM, CRAM, FASTA, FAI, VCF, .dict, TSV, CSV
    #[arg(required = true)]
    pub query: PathBuf,

    /// Reference file (the file to compare against)
    /// Supports: BAM, SAM, CRAM, FASTA, FAI, VCF, .dict, TSV, CSV
    #[arg(required = true)]
    pub reference: PathBuf,

    /// Also compute the reverse comparison (reference as query, query as reference).
    /// By default, scoring is asymmetric: it measures how well the query matches
    /// the reference. With --symmetric, both directions are computed.
    #[arg(long)]
    pub symmetric: bool,

    // === Scoring weight options ===
    /// Weight for contig match score (0-100, default 70)
    /// How well query contigs match reference contigs
    #[arg(long, default_value = "70", value_parser = clap::value_parser!(u32).range(0..=100))]
    pub weight_match: u32,

    /// Weight for coverage score (0-100, default 20)
    /// What fraction of reference contigs are covered by query
    #[arg(long, default_value = "20", value_parser = clap::value_parser!(u32).range(0..=100))]
    pub weight_coverage: u32,

    /// Weight for order score (0-100, default 10)
    /// Whether contigs appear in the same order
    #[arg(long, default_value = "10", value_parser = clap::value_parser!(u32).range(0..=100))]
    pub weight_order: u32,
}

/// Result of scoring in one direction
struct ScoreResult {
    query_path: PathBuf,
    reference_path: PathBuf,
    query_header: QueryHeader,
    reference_header: QueryHeader,
    score: MatchScore,
}

/// Execute the score command
///
/// # Errors
///
/// Returns an error if inputs cannot be parsed or comparison fails.
#[allow(clippy::needless_pass_by_value)]
pub fn run(args: ScoreArgs, format: OutputFormat, verbose: bool) -> anyhow::Result<()> {
    // Build scoring weights from command line args
    let weights = ScoringWeights {
        contig_match: f64::from(args.weight_match) / 100.0,
        coverage: f64::from(args.weight_coverage) / 100.0,
        order: f64::from(args.weight_order) / 100.0,
        conflict_penalty: 0.1,
    };

    if verbose {
        eprintln!(
            "Scoring weights: {:.0}% match, {:.0}% coverage, {:.0}% order",
            weights.contig_match * 100.0,
            weights.coverage * 100.0,
            weights.order * 100.0,
        );
    }

    // Parse query file
    let query_header = parse_input(&args.query)?;
    if verbose {
        eprintln!(
            "Query: {} contigs ({:.0}% have MD5)",
            query_header.contigs.len(),
            query_header.md5_coverage() * 100.0,
        );
    }

    // Parse reference file
    let reference_header = parse_input(&args.reference)?;
    if verbose {
        eprintln!(
            "Reference: {} contigs ({:.0}% have MD5)",
            reference_header.contigs.len(),
            reference_header.md5_coverage() * 100.0,
        );
    }

    // Compute forward score (query vs reference)
    let forward_result = compute_score(
        args.query.clone(),
        args.reference.clone(),
        &query_header,
        &reference_header,
        &weights,
    );

    // Compute reverse score if symmetric
    let reverse_result = if args.symmetric {
        Some(compute_score(
            args.reference.clone(),
            args.query.clone(),
            &reference_header,
            &query_header,
            &weights,
        ))
    } else {
        None
    };

    // Output results
    match format {
        OutputFormat::Text => {
            print_text_result(&forward_result, &weights, "");
            if let Some(ref reverse) = reverse_result {
                println!("\n{}", "─".repeat(60));
                print_text_result(reverse, &weights, " (reverse)");
            }
        }
        OutputFormat::Json => {
            print_json_results(&forward_result, reverse_result.as_ref(), &weights)?;
        }
        OutputFormat::Tsv => {
            print_tsv_results(&forward_result, reverse_result.as_ref(), &weights);
        }
    }

    Ok(())
}

fn compute_score(
    query_path: PathBuf,
    reference_path: PathBuf,
    query_header: &QueryHeader,
    reference_header: &QueryHeader,
    weights: &ScoringWeights,
) -> ScoreResult {
    // Convert reference header to KnownReference for scoring
    let reference = KnownReference::new(
        "reference",
        reference_path.display().to_string().as_str(),
        Assembly::Other("unknown".to_string()),
        ReferenceSource::Custom("file".to_string()),
    )
    .with_contigs(reference_header.contigs.clone());

    let score = MatchScore::calculate_with_weights(query_header, &reference, weights);

    ScoreResult {
        query_path,
        reference_path,
        query_header: query_header.clone(),
        reference_header: reference_header.clone(),
        score,
    }
}

fn parse_input(path: &Path) -> anyhow::Result<QueryHeader> {
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .map(str::to_lowercase);

    match ext.as_deref() {
        Some("dict") => Ok(parsing::dict::parse_dict_file(path)?),
        Some("fai") => Ok(parsing::fai::parse_fai_file(path)?),
        Some("fa" | "fasta" | "fna") => Ok(parsing::fasta::parse_fasta_file(path)?),
        Some("vcf" | "vcf.gz") => Ok(parsing::vcf::parse_vcf_file(path)?),
        Some("tsv") => Ok(parsing::tsv::parse_tsv_file(path, '\t')?),
        Some("csv") => Ok(parsing::tsv::parse_tsv_file(path, ',')?),
        // Default to SAM/BAM/CRAM parsing
        _ => Ok(parsing::sam::parse_file(path)?),
    }
}

fn print_text_result(result: &ScoreResult, weights: &ScoringWeights, suffix: &str) {
    let norm = weights.normalized();

    println!("\nScoring{}: {} vs {}", suffix, result.query_path.display(), result.reference_path.display());

    // Score breakdown
    println!(
        "\n   Score: {:.1}% = {:.0}%×match + {:.0}%×coverage + {:.0}%×order",
        result.score.composite * 100.0,
        result.score.contig_match_score * 100.0,
        result.score.coverage_score * 100.0,
        result.score.order_score * 100.0,
    );
    println!(
        "          (weights: {:.0}% match, {:.0}% coverage, {:.0}% order)",
        norm.contig_match * 100.0,
        norm.coverage * 100.0,
        norm.order * 100.0,
    );

    // Query contigs
    let total_query = result.query_header.contigs.len();
    let exact = result.score.exact_matches;
    let name_len = result.score.name_length_matches;
    let conflicts = result.score.md5_conflicts;
    let unmatched = result.score.unmatched;

    println!(
        "\n   Query contigs: {total_query} total → {exact} exact, {name_len} name+length, {conflicts} conflicts, {unmatched} unmatched"
    );

    // Reference coverage
    let total_ref = result.reference_header.contigs.len();
    let matched_ref = exact + name_len;
    let uncovered_ref = total_ref.saturating_sub(matched_ref);
    println!(
        "   Reference contigs: {total_ref} total, {matched_ref} matched, {uncovered_ref} not in query"
    );

    // Order
    if !result.score.order_preserved {
        println!("   Order: DIFFERENT from reference");
    }

    // Confidence
    println!("   Confidence: {:?}", result.score.confidence);
}

fn print_json_results(
    forward: &ScoreResult,
    reverse: Option<&ScoreResult>,
    weights: &ScoringWeights,
) -> anyhow::Result<()> {
    let norm = weights.normalized();

    let make_result_json = |result: &ScoreResult| {
        let ref_total = result.reference_header.contigs.len();
        let ref_matched = result.score.exact_matches + result.score.name_length_matches;
        let ref_uncovered = ref_total.saturating_sub(ref_matched);

        serde_json::json!({
            "query": {
                "file": result.query_path.display().to_string(),
                "contigs": result.query_header.contigs.len(),
            },
            "reference": {
                "file": result.reference_path.display().to_string(),
                "contigs": result.reference_header.contigs.len(),
            },
            "score": {
                "composite": result.score.composite,
                "confidence": format!("{:?}", result.score.confidence),
                "contig_match_score": result.score.contig_match_score,
                "coverage_score": result.score.coverage_score,
                "order_score": result.score.order_score,
                "weights": {
                    "match": norm.contig_match,
                    "coverage": norm.coverage,
                    "order": norm.order,
                },
            },
            "query_contigs": {
                "exact_matches": result.score.exact_matches,
                "name_length_matches": result.score.name_length_matches,
                "md5_conflicts": result.score.md5_conflicts,
                "unmatched": result.score.unmatched,
            },
            "reference_coverage": {
                "total": ref_total,
                "matched": ref_matched,
                "not_in_query": ref_uncovered,
            },
            "order_preserved": result.score.order_preserved,
        })
    };

    let output = if let Some(rev) = reverse {
        serde_json::json!({
            "forward": make_result_json(forward),
            "reverse": make_result_json(rev),
        })
    } else {
        make_result_json(forward)
    };

    println!("{}", serde_json::to_string_pretty(&output)?);
    Ok(())
}

fn print_tsv_results(forward: &ScoreResult, reverse: Option<&ScoreResult>, weights: &ScoringWeights) {
    let norm = weights.normalized();

    // Header
    println!(
        "direction\tquery\treference\tscore\tmatch_score\tcoverage_score\torder_score\tweight_match\tweight_coverage\tweight_order\tconfidence\texact\tname_length\tconflicts\tunmatched\tref_total\tref_matched\tref_uncovered"
    );

    let print_row = |direction: &str, result: &ScoreResult| {
        let ref_total = result.reference_header.contigs.len();
        let ref_matched = result.score.exact_matches + result.score.name_length_matches;
        let ref_uncovered = ref_total.saturating_sub(ref_matched);

        println!(
            "{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.2}\t{:.2}\t{:.2}\t{:?}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            direction,
            result.query_path.display(),
            result.reference_path.display(),
            result.score.composite,
            result.score.contig_match_score,
            result.score.coverage_score,
            result.score.order_score,
            norm.contig_match,
            norm.coverage,
            norm.order,
            result.score.confidence,
            result.score.exact_matches,
            result.score.name_length_matches,
            result.score.md5_conflicts,
            result.score.unmatched,
            ref_total,
            ref_matched,
            ref_uncovered,
        );
    };

    print_row("forward", forward);
    if let Some(rev) = reverse {
        print_row("reverse", rev);
    }
}
