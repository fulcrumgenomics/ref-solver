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

    println!(
        "\nScoring{}: {} vs {}",
        suffix,
        result.query_path.display(),
        result.reference_path.display()
    );

    // Score breakdown
    println!(
        "\n   Score: {:.1}% = {:.0}%×match + {:.0}%×coverage + {:.0}%×order",
        result.score.composite * 100.0,
        result.score.match_quality * 100.0,
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
                "match_quality": result.score.match_quality,
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

fn print_tsv_results(
    forward: &ScoreResult,
    reverse: Option<&ScoreResult>,
    weights: &ScoringWeights,
) {
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
            result.score.match_quality,
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::contig::Contig;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_temp_dict_file(contigs: &[(&str, u64, Option<&str>)]) -> NamedTempFile {
        let mut file = NamedTempFile::with_suffix(".dict").unwrap();
        writeln!(file, "@HD\tVN:1.0\tSO:unsorted").unwrap();
        for (name, len, md5) in contigs {
            let md5_field = md5.map(|m| format!("\tM5:{m}")).unwrap_or_default();
            writeln!(file, "@SQ\tSN:{name}\tLN:{len}{md5_field}").unwrap();
        }
        file.flush().unwrap();
        file
    }

    #[test]
    fn test_parse_dict_input() {
        // MD5 must be exactly 32 hex characters
        let valid_md5 = "6aef897c3d6ff0c78aff06ac189178dd";
        let file = create_temp_dict_file(&[("chr1", 1000, Some(valid_md5)), ("chr2", 2000, None)]);

        let header = parse_input(file.path()).unwrap();
        assert_eq!(header.contigs.len(), 2);
        assert_eq!(header.contigs[0].name, "chr1");
        assert_eq!(header.contigs[0].length, 1000);
        assert_eq!(header.contigs[0].md5.as_deref(), Some(valid_md5));
        assert_eq!(header.contigs[1].name, "chr2");
        assert_eq!(header.contigs[1].length, 2000);
        assert!(header.contigs[1].md5.is_none());
    }

    #[test]
    fn test_compute_score_perfect_match() {
        let query_header =
            QueryHeader::new(vec![Contig::new("chr1", 1000), Contig::new("chr2", 2000)]);
        let reference_header =
            QueryHeader::new(vec![Contig::new("chr1", 1000), Contig::new("chr2", 2000)]);

        let weights = ScoringWeights::default();
        let result = compute_score(
            PathBuf::from("query.dict"),
            PathBuf::from("reference.dict"),
            &query_header,
            &reference_header,
            &weights,
        );

        // Perfect match: all contigs match by name+length
        assert_eq!(result.score.name_length_matches, 2);
        assert_eq!(result.score.unmatched, 0);
        assert!(
            result.score.composite > 0.9,
            "Perfect match should score > 90%"
        );
    }

    #[test]
    fn test_compute_score_partial_match() {
        let query_header = QueryHeader::new(vec![
            Contig::new("chr1", 1000),
            Contig::new("chr2", 2000),
            Contig::new("chr3", 3000),
        ]);
        let reference_header =
            QueryHeader::new(vec![Contig::new("chr1", 1000), Contig::new("chr2", 2000)]);

        let weights = ScoringWeights::default();
        let result = compute_score(
            PathBuf::from("query.dict"),
            PathBuf::from("reference.dict"),
            &query_header,
            &reference_header,
            &weights,
        );

        // 2 of 3 query contigs match
        assert_eq!(result.score.name_length_matches, 2);
        assert_eq!(result.score.unmatched, 1);
        assert!(
            result.score.match_quality < 1.0,
            "Partial match should have match_quality < 1.0"
        );
    }

    #[test]
    fn test_compute_score_asymmetric() {
        // Query has fewer contigs than reference
        let query_header = QueryHeader::new(vec![Contig::new("chr1", 1000)]);
        let reference_header = QueryHeader::new(vec![
            Contig::new("chr1", 1000),
            Contig::new("chr2", 2000),
            Contig::new("chr3", 3000),
        ]);

        let weights = ScoringWeights::default();

        // Forward: query → reference
        let forward = compute_score(
            PathBuf::from("query.dict"),
            PathBuf::from("reference.dict"),
            &query_header,
            &reference_header,
            &weights,
        );

        // Reverse: reference → query
        let reverse = compute_score(
            PathBuf::from("reference.dict"),
            PathBuf::from("query.dict"),
            &reference_header,
            &query_header,
            &weights,
        );

        // Forward: all query contigs match (100% contig match), but only 1/3 ref covered
        assert_eq!(forward.score.name_length_matches, 1);
        assert!(
            (forward.score.match_quality - 1.0).abs() < 0.001,
            "All query contigs match"
        );
        assert!(
            forward.score.coverage_score < 0.5,
            "Coverage should be 1/3 = 0.33"
        );

        // Reverse: only 1/3 query contigs match, but reference is fully covered
        assert_eq!(reverse.score.unmatched, 2);
        assert!(
            reverse.score.match_quality < 0.5,
            "Only 1/3 query contigs match"
        );
        assert!(
            (reverse.score.coverage_score - 1.0).abs() < 0.001,
            "Reference fully covered"
        );
    }

    #[test]
    fn test_custom_weights() {
        let query_header =
            QueryHeader::new(vec![Contig::new("chr1", 1000), Contig::new("chr2", 2000)]);
        let reference_header = QueryHeader::new(vec![
            Contig::new("chr1", 1000),
            Contig::new("chr2", 2000),
            Contig::new("chr3", 3000),
        ]);

        // Emphasize coverage over match
        let high_coverage_weights = ScoringWeights {
            contig_match: 0.2,
            coverage: 0.7,
            order: 0.1,
            conflict_penalty: 0.1,
        };

        // Emphasize match over coverage
        let high_match_weights = ScoringWeights {
            contig_match: 0.8,
            coverage: 0.1,
            order: 0.1,
            conflict_penalty: 0.1,
        };

        let result_high_cov = compute_score(
            PathBuf::from("q.dict"),
            PathBuf::from("r.dict"),
            &query_header,
            &reference_header,
            &high_coverage_weights,
        );

        let result_high_match = compute_score(
            PathBuf::from("q.dict"),
            PathBuf::from("r.dict"),
            &query_header,
            &reference_header,
            &high_match_weights,
        );

        // With emphasis on coverage (2/3), score should be lower
        // With emphasis on match (100%), score should be higher
        assert!(
            result_high_match.score.composite > result_high_cov.score.composite,
            "High match weight should yield higher score when matches are 100% but coverage is 66%"
        );
    }

    #[test]
    fn test_score_with_md5_match() {
        let query_header = QueryHeader::new(vec![
            Contig::new("chr1", 1000).with_md5("abc123"),
            Contig::new("chr2", 2000).with_md5("def456"),
        ]);
        let reference_header = QueryHeader::new(vec![
            Contig::new("chr1", 1000).with_md5("abc123"),
            Contig::new("chr2", 2000).with_md5("def456"),
        ]);

        let weights = ScoringWeights::default();
        let result = compute_score(
            PathBuf::from("query.dict"),
            PathBuf::from("reference.dict"),
            &query_header,
            &reference_header,
            &weights,
        );

        // With matching MD5s, these should be exact matches
        assert_eq!(result.score.exact_matches, 2);
        assert_eq!(result.score.name_length_matches, 0);
        assert!(
            result.score.composite > 0.95,
            "Exact MD5 match should score very high"
        );
    }

    #[test]
    fn test_score_with_md5_conflict() {
        let query_header = QueryHeader::new(vec![
            Contig::new("chr1", 1000).with_md5("abc123"),
            Contig::new("chr2", 2000).with_md5("def456"),
        ]);
        let reference_header = QueryHeader::new(vec![
            Contig::new("chr1", 1000).with_md5("DIFFERENT1"),
            Contig::new("chr2", 2000).with_md5("DIFFERENT2"),
        ]);

        let weights = ScoringWeights::default();
        let result = compute_score(
            PathBuf::from("query.dict"),
            PathBuf::from("reference.dict"),
            &query_header,
            &reference_header,
            &weights,
        );

        // MD5 conflicts should be penalized
        assert_eq!(result.score.md5_conflicts, 2);
        assert_eq!(result.score.exact_matches, 0);
        assert!(
            result.score.composite < 0.3,
            "MD5 conflicts should result in low score, got {:.1}%",
            result.score.composite * 100.0
        );
    }
}
