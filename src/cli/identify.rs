use std::path::{Path, PathBuf};

use clap::Args;

use crate::catalog::hierarchical::HierarchicalCatalog;
use crate::catalog::store::ReferenceCatalog;
use crate::cli::OutputFormat;
use crate::core::header::QueryHeader;
use crate::core::types::Confidence;
use crate::matching::engine::{MatchResult, MatchingConfig, MatchingEngine, ScoringWeights};
use crate::matching::hierarchical_engine::{HierarchicalMatchResult, HierarchicalMatchingEngine};
use crate::matching::Suggestion;
use crate::parsing;

/// How to handle references that have contigs missing from their FASTA
/// (e.g., CHM13 where MT is in assembly report but uses standard rCRS mitochondria)
#[derive(Clone, Copy, Debug, Default, clap::ValueEnum)]
pub enum MissingContigHandling {
    /// Show warnings when query has contigs that match a reference's missing contigs
    #[default]
    Warn,
    /// Treat missing contigs as errors that lower match confidence
    Strict,
    /// Don't mention missing contigs at all
    Silent,
}

#[derive(Args)]
pub struct IdentifyArgs {
    /// Input file (BAM, SAM, CRAM, FASTA, FAI, VCF, .dict, TSV, or CSV)
    /// Use '-' for stdin (expects header text)
    #[arg(required = true)]
    pub input: PathBuf,

    /// Input format (auto-detected by default)
    #[arg(long)]
    pub input_format: Option<InputFormat>,

    /// Number of matches to show
    #[arg(short = 'n', long, default_value = "5")]
    pub max_matches: usize,

    /// Only show exact or near-exact matches
    #[arg(long)]
    pub exact_only: bool,

    /// Path to custom catalog file
    #[arg(long)]
    pub catalog: Option<PathBuf>,

    /// Use hierarchical catalog format (required when --catalog points to a hierarchical catalog)
    #[arg(long)]
    pub hierarchical: bool,

    /// How to handle references with contigs missing from FASTA
    /// (e.g., CHM13 MT which uses standard rCRS mitochondria)
    #[arg(long, value_enum, default_value = "warn")]
    pub missing_contig_handling: MissingContigHandling,

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

#[derive(Clone, Copy, Debug, clap::ValueEnum)]
pub enum InputFormat {
    Sam,
    Bam,
    Cram,
    Dict,
    Fai,
    Fasta,
    Vcf,
    Tsv,
    Csv,
}

/// Execute identify subcommand
///
/// # Errors
///
/// Returns an error if the input cannot be parsed or identification fails.
#[allow(clippy::needless_pass_by_value)] // CLI entry point, values from clap
pub fn run(args: IdentifyArgs, format: OutputFormat, verbose: bool) -> anyhow::Result<()> {
    // Parse input first (needed for both catalog types)
    let query = parse_input(&args)?;

    if verbose {
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)] // Percentage 0-100
        let md5_pct = (query.md5_coverage() * 100.0) as u32;
        eprintln!(
            "Parsed {} contigs from input ({md5_pct}% have MD5)",
            query.contigs.len(),
        );
    }

    // Use hierarchical or flat catalog based on flag
    if args.hierarchical {
        run_hierarchical(&args, &query, format, verbose)
    } else {
        run_flat(&args, &query, format, verbose)
    }
}

fn run_flat(
    args: &IdentifyArgs,
    query: &QueryHeader,
    format: OutputFormat,
    verbose: bool,
) -> anyhow::Result<()> {
    // Load flat catalog
    let catalog = if let Some(path) = &args.catalog {
        ReferenceCatalog::load_from_file(path)?
    } else {
        ReferenceCatalog::load_embedded()?
    };

    if verbose {
        eprintln!("Loaded flat catalog with {} references", catalog.len());
    }

    if catalog.is_empty() {
        eprintln!("Warning: Catalog is empty, no references to match against.");
        return Ok(());
    }

    // Build scoring weights from command line args
    let scoring_weights = ScoringWeights {
        contig_match: f64::from(args.weight_match) / 100.0,
        coverage: f64::from(args.weight_coverage) / 100.0,
        order: f64::from(args.weight_order) / 100.0,
        conflict_penalty: 0.1, // Default: 10% credit for MD5 conflicts
    };

    if verbose {
        eprintln!(
            "Scoring weights: {:.0}% match, {:.0}% coverage, {:.0}% order",
            scoring_weights.contig_match * 100.0,
            scoring_weights.coverage * 100.0,
            scoring_weights.order * 100.0,
        );
    }

    // Find matches with custom config
    let config = MatchingConfig {
        min_score: 0.1,
        scoring_weights: scoring_weights.clone(),
    };
    let engine = MatchingEngine::new(&catalog, config);
    let matches = engine.find_matches(query, args.max_matches);

    if matches.is_empty() {
        eprintln!("No matching references found.");
        return Ok(());
    }

    // Output results
    match format {
        OutputFormat::Text => {
            print_text_results(
                &matches,
                query,
                verbose,
                args.missing_contig_handling,
                &scoring_weights,
            );
        }
        OutputFormat::Json => {
            print_json_results(&matches, args.missing_contig_handling, &scoring_weights)?;
        }
        OutputFormat::Tsv => print_tsv_results(&matches, &scoring_weights),
    }

    Ok(())
}

fn run_hierarchical(
    args: &IdentifyArgs,
    query: &QueryHeader,
    format: OutputFormat,
    verbose: bool,
) -> anyhow::Result<()> {
    // Load hierarchical catalog
    let catalog_path = args
        .catalog
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("--catalog is required when using --hierarchical"))?;

    let catalog = HierarchicalCatalog::load(catalog_path)?;

    if verbose {
        eprintln!(
            "Loaded hierarchical catalog v{} with {} assemblies",
            catalog.version,
            catalog.assemblies.len()
        );
    }

    // Find matches
    let engine = HierarchicalMatchingEngine::new(&catalog);
    let matches = engine.find_matches(query, args.max_matches);

    if matches.is_empty() {
        eprintln!("No matching references found.");
        return Ok(());
    }

    // Output results
    match format {
        OutputFormat::Text => print_hierarchical_text_results(&matches, query, verbose),
        OutputFormat::Json => print_hierarchical_json_results(&matches)?,
        OutputFormat::Tsv => print_hierarchical_tsv_results(&matches),
    }

    Ok(())
}

fn parse_input(args: &IdentifyArgs) -> anyhow::Result<QueryHeader> {
    use std::io::{self, Read};

    // Handle stdin
    if args.input.to_string_lossy() == "-" {
        let mut buffer = String::new();
        io::stdin().read_to_string(&mut buffer)?;
        return Ok(parsing::sam::parse_header_text(&buffer)?);
    }

    // Auto-detect or use specified format
    let format = args
        .input_format
        .unwrap_or_else(|| detect_format(&args.input));

    match format {
        InputFormat::Sam | InputFormat::Bam | InputFormat::Cram => {
            Ok(parsing::sam::parse_file(&args.input)?)
        }
        InputFormat::Dict => Ok(parsing::dict::parse_dict_file(&args.input)?),
        InputFormat::Fai => Ok(parsing::fai::parse_fai_file(&args.input)?),
        InputFormat::Fasta => Ok(parsing::fasta::parse_fasta_file(&args.input)?),
        InputFormat::Vcf => Ok(parsing::vcf::parse_vcf_file(&args.input)?),
        InputFormat::Tsv => Ok(parsing::tsv::parse_tsv_file(&args.input, '\t')?),
        InputFormat::Csv => Ok(parsing::tsv::parse_tsv_file(&args.input, ',')?),
    }
}

/// Detect input format from file extension
fn detect_format(path: &Path) -> InputFormat {
    let path_str = path.to_string_lossy().to_lowercase();

    // Check for FASTA files (including gzipped)
    if parsing::fasta::is_fasta_file(path) {
        return InputFormat::Fasta;
    }

    // Check for gzipped VCF
    if path_str.ends_with(".vcf.gz") || path_str.ends_with(".vcf.bgz") {
        return InputFormat::Vcf;
    }

    // Get the extension for simple cases
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .map(str::to_lowercase);

    match ext.as_deref() {
        Some("bam") => InputFormat::Bam,
        Some("cram") => InputFormat::Cram,
        Some("dict") => InputFormat::Dict,
        Some("fai") => InputFormat::Fai,
        Some("vcf") => InputFormat::Vcf,
        Some("tsv") => InputFormat::Tsv,
        Some("csv") => InputFormat::Csv,
        _ => InputFormat::Sam, // Default to SAM for unknown extensions
    }
}

#[allow(clippy::too_many_lines)] // TODO: Refactor into smaller functions
fn print_text_results(
    matches: &[MatchResult],
    query: &QueryHeader,
    verbose: bool,
    missing_handling: MissingContigHandling,
    weights: &ScoringWeights,
) {
    for (i, result) in matches.iter().enumerate() {
        if i > 0 {
            println!("\n{}", "─".repeat(60));
        }

        // Header
        let confidence_str = match result.score.confidence {
            Confidence::Exact => "EXACT",
            Confidence::High => "HIGH",
            Confidence::Medium => "MEDIUM",
            Confidence::Low => "LOW",
        };

        println!(
            "\n#{} {} ({})",
            i + 1,
            result.reference.display_name,
            confidence_str
        );
        println!("   ID: {}", result.reference.id);
        println!("   Assembly: {}", result.reference.assembly);
        println!("   Source: {}", result.reference.source);
        println!("   Match Type: {:?}", result.diagnosis.match_type);

        // Score breakdown: show component scores and final composite
        // Normalize weights for display
        let norm = weights.normalized();
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

        // Check for contigs missing from FASTA
        if !result.reference.contigs_missing_from_fasta.is_empty() {
            // Find query contigs that match the missing contigs (by name, case-insensitive)
            let missing_set: std::collections::HashSet<String> = result
                .reference
                .contigs_missing_from_fasta
                .iter()
                .map(|s| s.to_lowercase())
                .collect();

            let query_has_missing: Vec<&str> = query
                .contigs
                .iter()
                .filter(|c| missing_set.contains(&c.name.to_lowercase()))
                .map(|c| c.name.as_str())
                .collect();

            match missing_handling {
                MissingContigHandling::Silent => {}
                MissingContigHandling::Warn => {
                    if !query_has_missing.is_empty() {
                        println!(
                            "\n   Warning: Query has contig(s) not in reference FASTA: {}",
                            query_has_missing.join(", ")
                        );
                        println!(
                            "   Note: {} uses external sequence(s) for: {}",
                            result.reference.display_name,
                            result.reference.contigs_missing_from_fasta.join(", ")
                        );
                    }
                }
                MissingContigHandling::Strict => {
                    if !query_has_missing.is_empty() {
                        println!(
                            "\n   ERROR: Query has contig(s) not in reference FASTA: {}",
                            query_has_missing.join(", ")
                        );
                        println!(
                            "   The reference {} does not include: {}",
                            result.reference.display_name,
                            result.reference.contigs_missing_from_fasta.join(", ")
                        );
                    }
                }
            }
        }

        // Match details - query contigs
        let total_query = query.contigs.len();
        let exact = result.score.exact_matches;
        let name_len = result.score.name_length_matches;
        let conflicts = result.score.md5_conflicts;
        let unmatched = result.score.unmatched;

        println!(
            "\n   Query contigs: {total_query} total → {exact} exact, {name_len} name+length, {conflicts} conflicts, {unmatched} unmatched"
        );

        // Reference coverage info
        let total_ref = result.reference.contigs.len();
        let matched_ref = exact + name_len; // Good matches that cover reference
        let uncovered_ref = total_ref.saturating_sub(matched_ref);
        println!(
            "   Reference contigs: {total_ref} total, {matched_ref} matched, {uncovered_ref} not in query"
        );

        if result.diagnosis.reordered {
            println!("   Order: DIFFERENT from reference");
        }

        // Conflicts
        if !result.diagnosis.conflicts.is_empty() {
            println!("\n   Conflicts:");
            for conflict in &result.diagnosis.conflicts {
                println!("   - {}", conflict.description);
            }
        }

        // Suggestions
        if !result.diagnosis.suggestions.is_empty() {
            println!("\n   Suggestions:");
            for suggestion in &result.diagnosis.suggestions {
                match suggestion {
                    Suggestion::RenameContigs { command_hint, .. } => {
                        println!("   - Rename contigs:");
                        for line in command_hint.lines() {
                            println!("     {line}");
                        }
                    }
                    Suggestion::ReorderContigs { command_hint } => {
                        println!("   - Reorder contigs:");
                        for line in command_hint.lines() {
                            println!("     {line}");
                        }
                    }
                    Suggestion::ReplaceContig {
                        contig_name,
                        reason,
                        ..
                    } => {
                        println!("   - Replace {contig_name}: {reason}");
                    }
                    Suggestion::UseAsIs { warnings } => {
                        if warnings.is_empty() {
                            println!("   - Safe to use as-is");
                        } else {
                            println!("   - Safe to use with warnings:");
                            for w in warnings {
                                println!("     - {w}");
                            }
                        }
                    }
                    Suggestion::Realign {
                        reason,
                        suggested_reference,
                    } => {
                        println!("   - Realignment needed: {reason}");
                        println!("     Suggested reference: {suggested_reference}");
                    }
                }
            }
        }

        // Download URL
        if let Some(url) = &result.reference.download_url {
            println!("\n   Download: {url}");
        }

        // Verbose details
        if verbose && !result.diagnosis.renamed_matches.is_empty() {
            println!("\n   Rename mappings:");
            for r in &result.diagnosis.renamed_matches {
                println!("     {} -> {}", r.query_name, r.reference_name);
            }
        }
    }

    println!();
}

fn print_json_results(
    matches: &[MatchResult],
    missing_handling: MissingContigHandling,
    weights: &ScoringWeights,
) -> anyhow::Result<()> {
    let norm = weights.normalized();
    // Create serializable output
    let output: Vec<serde_json::Value> = matches
        .iter()
        .map(|m| {
            // Calculate reference coverage
            let ref_total = m.reference.contigs.len();
            let ref_matched = m.score.exact_matches + m.score.name_length_matches;
            let ref_uncovered = ref_total.saturating_sub(ref_matched);

            let mut json = serde_json::json!({
                "reference": {
                    "id": m.reference.id.0,
                    "display_name": m.reference.display_name,
                    "assembly": format!("{}", m.reference.assembly),
                    "source": format!("{}", m.reference.source),
                    "download_url": m.reference.download_url,
                    "total_contigs": ref_total,
                },
                "score": {
                    "composite": m.score.composite,
                    "confidence": format!("{:?}", m.score.confidence),
                    // Component scores (these make up the composite)
                    "match_quality": m.score.match_quality,
                    "coverage_score": m.score.coverage_score,
                    "order_score": m.score.order_score,
                    // Weights used
                    "weights": {
                        "match": norm.contig_match,
                        "coverage": norm.coverage,
                        "order": norm.order,
                    },
                },
                "query_contigs": {
                    "exact_matches": m.score.exact_matches,
                    "name_length_matches": m.score.name_length_matches,
                    "md5_conflicts": m.score.md5_conflicts,
                    "unmatched": m.score.unmatched,
                },
                "reference_coverage": {
                    "total": ref_total,
                    "matched": ref_matched,
                    "not_in_query": ref_uncovered,
                },
                "match_type": format!("{:?}", m.diagnosis.match_type),
                "reordered": m.diagnosis.reordered,
            });

            // Add missing contig info unless silent
            if !matches!(missing_handling, MissingContigHandling::Silent)
                && !m.reference.contigs_missing_from_fasta.is_empty()
            {
                json["reference"]["contigs_missing_from_fasta"] =
                    serde_json::json!(&m.reference.contigs_missing_from_fasta);
            }

            json
        })
        .collect();

    println!("{}", serde_json::to_string_pretty(&output)?);
    Ok(())
}

fn print_tsv_results(matches: &[MatchResult], weights: &ScoringWeights) {
    let norm = weights.normalized();
    // Header with all fields
    println!(
        "rank\tid\tdisplay_name\tassembly\tsource\tmatch_type\tscore\tmatch_score\tcoverage_score\torder_score\tweight_match\tweight_coverage\tweight_order\tconfidence\texact\tname_length\tconflicts\tunmatched\tref_total\tref_matched\tref_uncovered"
    );
    for (i, m) in matches.iter().enumerate() {
        let ref_total = m.reference.contigs.len();
        let ref_matched = m.score.exact_matches + m.score.name_length_matches;
        let ref_uncovered = ref_total.saturating_sub(ref_matched);

        println!(
            "{}\t{}\t{}\t{}\t{}\t{:?}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.2}\t{:.2}\t{:.2}\t{:?}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            i + 1,
            m.reference.id,
            m.reference.display_name,
            m.reference.assembly,
            m.reference.source,
            m.diagnosis.match_type,
            m.score.composite,
            m.score.match_quality,
            m.score.coverage_score,
            m.score.order_score,
            norm.contig_match,
            norm.coverage,
            norm.order,
            m.score.confidence,
            m.score.exact_matches,
            m.score.name_length_matches,
            m.score.md5_conflicts,
            m.score.unmatched,
            ref_total,
            ref_matched,
            ref_uncovered,
        );
    }
}

// ============================================================================
// Hierarchical catalog output functions
// ============================================================================

fn print_hierarchical_text_results(
    matches: &[HierarchicalMatchResult],
    query: &QueryHeader,
    verbose: bool,
) {
    for (i, result) in matches.iter().enumerate() {
        if i > 0 {
            println!("\n{}", "─".repeat(60));
        }

        // Header with match type
        let match_str = format!("{:?}", result.match_type).to_uppercase();
        println!("\n#{} {} ({})", i + 1, result.display_name, match_str);

        // Distribution info
        println!("   Distribution ID: {}", result.distribution_id);

        // Assembly info (if available)
        if !result.assembly_id.is_empty() {
            println!(
                "   Assembly: {} ({})",
                result.assembly_name, result.assembly_id
            );
            if !result.version_string.is_empty() {
                println!(
                    "   Version: {} ({})",
                    result.version_string, result.version_id
                );
            }
        }

        // Match score
        println!("   Score: {:.1}%", result.match_percentage());

        // Contig summary
        println!("\n   Contig Summary:");
        println!("   - Your file: {} contigs", result.total_query_contigs);
        println!(
            "   - This distribution: {} contigs",
            result.total_distribution_contigs
        );
        println!("   - Matched: {} contigs", result.matched_contigs);

        if result.extra_in_query > 0 {
            println!("   - Extra in your file: {}", result.extra_in_query);
        }
        if result.missing_from_query > 0 {
            println!("   - Missing from your file: {}", result.missing_from_query);
        }

        // Presence breakdown (only show if there's assembly linkage)
        let counts = &result.presence_counts;
        if counts.in_both > 0 || counts.fasta_only > 0 || counts.report_only > 0 {
            println!("\n   Presence Breakdown:");
            if counts.in_both > 0 {
                println!("   - In both (FASTA + report): {} contigs", counts.in_both);
            }
            if counts.fasta_only > 0 {
                println!("   - FASTA-only (decoy/HLA): {} contigs", counts.fasta_only);
            }
            if counts.report_only > 0 {
                println!(
                    "   - Report-only (not in FASTA): {} contigs",
                    counts.report_only
                );
            }
        }

        // Verbose details
        if verbose {
            println!("\n   Query contigs: {}", query.contigs.len());
            let md5_count = query.contigs.iter().filter(|c| c.md5.is_some()).count();
            println!("   Query contigs with MD5: {md5_count}");
        }
    }

    println!();
}

fn print_hierarchical_json_results(matches: &[HierarchicalMatchResult]) -> anyhow::Result<()> {
    let output: Vec<serde_json::Value> = matches
        .iter()
        .map(|m| {
            serde_json::json!({
                "distribution": {
                    "id": m.distribution_id,
                    "display_name": m.display_name,
                },
                "assembly": {
                    "id": m.assembly_id,
                    "name": m.assembly_name,
                    "version_id": m.version_id,
                    "version": m.version_string,
                },
                "match_type": format!("{:?}", m.match_type),
                "score": m.score,
                "matched_contigs": m.matched_contigs,
                "total_query_contigs": m.total_query_contigs,
                "total_distribution_contigs": m.total_distribution_contigs,
                "extra_in_query": m.extra_in_query,
                "missing_from_query": m.missing_from_query,
                "presence_counts": {
                    "in_both": m.presence_counts.in_both,
                    "fasta_only": m.presence_counts.fasta_only,
                    "report_only": m.presence_counts.report_only,
                },
            })
        })
        .collect();

    println!("{}", serde_json::to_string_pretty(&output)?);
    Ok(())
}

fn print_hierarchical_tsv_results(matches: &[HierarchicalMatchResult]) {
    println!("rank\tdistribution_id\tdisplay_name\tassembly_id\tversion_id\tmatch_type\tscore\tmatched\tquery_total\tdist_total\tin_both\tfasta_only");
    for (i, m) in matches.iter().enumerate() {
        println!(
            "{}\t{}\t{}\t{}\t{}\t{:?}\t{:.4}\t{}\t{}\t{}\t{}\t{}",
            i + 1,
            m.distribution_id,
            m.display_name,
            m.assembly_id,
            m.version_id,
            m.match_type,
            m.score,
            m.matched_contigs,
            m.total_query_contigs,
            m.total_distribution_contigs,
            m.presence_counts.in_both,
            m.presence_counts.fasta_only,
        );
    }
}
