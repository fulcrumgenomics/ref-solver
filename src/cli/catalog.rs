use std::path::PathBuf;

use clap::{Args, Subcommand};

use crate::catalog::builder::{InputFormat, ReferenceBuilder};
use crate::catalog::hierarchical::HierarchicalCatalog;
use crate::catalog::store::ReferenceCatalog;
use crate::cli::OutputFormat;
use crate::core::types::{Assembly, ReferenceSource};

/// Helper function to convert usize count to f64 with explicit precision loss allowance
#[inline]
fn count_to_f64(count: usize) -> f64 {
    #[allow(clippy::cast_precision_loss)]
    {
        count as f64
    }
}

#[derive(Args)]
pub struct CatalogArgs {
    #[command(subcommand)]
    pub command: CatalogCommands,
}

#[derive(Subcommand)]
pub enum CatalogCommands {
    /// List all references in the catalog
    List {
        /// Path to custom catalog file
        #[arg(long)]
        catalog: Option<PathBuf>,

        /// Filter by assembly (e.g., "GRCh38")
        #[arg(long)]
        assembly: Option<String>,

        /// Filter by source (e.g., "UCSC")
        #[arg(long)]
        source: Option<String>,
    },

    /// Show details of a specific reference
    Show {
        /// Reference ID
        #[arg(required = true)]
        id: String,

        /// Path to custom catalog file
        #[arg(long)]
        catalog: Option<PathBuf>,

        /// Show all contigs
        #[arg(long)]
        all_contigs: bool,
    },

    /// Export the catalog to a file
    Export {
        /// Output file path
        #[arg(required = true)]
        output: PathBuf,

        /// Path to custom catalog file to export (defaults to embedded)
        #[arg(long)]
        catalog: Option<PathBuf>,
    },

    /// List hierarchical catalog contents (assemblies, versions, distributions)
    ListHierarchical {
        /// Path to hierarchical catalog file
        #[arg(required = true)]
        catalog: PathBuf,
    },

    /// Build a hierarchical catalog entry (FastaDistribution)
    BuildHierarchical {
        /// Distribution ID (e.g., "hg38_custom")
        #[arg(long, required = true)]
        id: String,

        /// Display name (e.g., "hg38 Custom Build")
        #[arg(long, required = true)]
        name: String,

        /// Input file(s) - can be specified multiple times
        #[arg(short, long = "input", required = true, num_args = 1..)]
        inputs: Vec<PathBuf>,

        /// Assembly ID to attach to (e.g., "grch38")
        #[arg(long)]
        assembly_id: Option<String>,

        /// Version ID to attach to (e.g., "grch38_p14")
        #[arg(long)]
        version_id: Option<String>,

        /// Source organization (ucsc, ncbi, broad, ensembl, 1kg, dragen, gdc, or custom)
        #[arg(long)]
        source: Option<String>,

        /// Reference FASTA download URL
        #[arg(long)]
        download_url: Option<String>,

        /// Tags (comma-separated)
        #[arg(long)]
        tags: Option<String>,

        /// Output file (creates new hierarchical catalog or standalone distribution JSON)
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Append to existing hierarchical catalog
        #[arg(long)]
        append_to: Option<PathBuf>,

        /// Overwrite if distribution ID already exists
        #[arg(long)]
        force: bool,

        /// Require MD5 checksums for all contigs
        #[arg(long)]
        require_md5: bool,

        /// Infer base assembly by matching MD5s against an existing catalog
        /// If no path given, uses the embedded catalog (or --append-to catalog)
        #[arg(long)]
        infer_assembly: Option<Option<PathBuf>>,
    },

    /// Build a new reference entry from input files
    Build {
        /// Unique reference ID (e.g., "grch38_custom")
        #[arg(long, required = true)]
        id: String,

        /// Display name (e.g., "GRCh38 Custom Build")
        #[arg(long, required = true)]
        name: String,

        /// Input file(s) - can be specified multiple times
        /// Supported formats: .dict, .fai, .sam, .bam, .cram, .vcf, _assembly_report.txt
        #[arg(short, long = "input", required = true, num_args = 1..)]
        inputs: Vec<PathBuf>,

        /// Assembly version (grch37, grch38, or custom name)
        #[arg(long)]
        assembly: Option<String>,

        /// Source organization (ucsc, ncbi, broad, ensembl, illumina, 1kg, or custom)
        #[arg(long)]
        source: Option<String>,

        /// Description text
        #[arg(long)]
        description: Option<String>,

        /// Reference FASTA download URL
        #[arg(long)]
        download_url: Option<String>,

        /// NCBI assembly report URL
        #[arg(long)]
        assembly_report_url: Option<String>,

        /// Comma-separated tags
        #[arg(long)]
        tags: Option<String>,

        /// Output file (JSON). If not specified, prints to stdout
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Append to existing catalog file
        #[arg(long)]
        append_to: Option<PathBuf>,

        /// Force overwrite if ID already exists in catalog
        #[arg(long)]
        force: bool,

        /// Force input format instead of auto-detection
        #[arg(long, value_enum)]
        input_format: Option<InputFormatArg>,

        /// Error if any contig lacks MD5 checksum
        #[arg(long)]
        require_md5: bool,
    },
}

/// Input format argument for CLI
#[derive(Clone, Copy, Debug, clap::ValueEnum)]
pub enum InputFormatArg {
    Dict,
    Fai,
    Fasta,
    NcbiReport,
    Sam,
    Bam,
    Cram,
    Vcf,
    Tsv,
}

impl From<InputFormatArg> for InputFormat {
    fn from(arg: InputFormatArg) -> Self {
        match arg {
            InputFormatArg::Dict => InputFormat::Dict,
            InputFormatArg::Fai => InputFormat::Fai,
            InputFormatArg::Fasta => InputFormat::Fasta,
            InputFormatArg::NcbiReport => InputFormat::NcbiReport,
            InputFormatArg::Sam => InputFormat::Sam,
            InputFormatArg::Bam => InputFormat::Bam,
            InputFormatArg::Cram => InputFormat::Cram,
            InputFormatArg::Vcf => InputFormat::Vcf,
            InputFormatArg::Tsv => InputFormat::Tsv,
        }
    }
}

pub fn run(args: CatalogArgs, format: OutputFormat, verbose: bool) -> anyhow::Result<()> {
    match args.command {
        CatalogCommands::List {
            catalog,
            assembly,
            source,
        } => run_list(
            catalog,
            assembly.as_deref(),
            source.as_deref(),
            format,
            verbose,
        ),
        CatalogCommands::Show {
            id,
            catalog,
            all_contigs,
        } => run_show(id, catalog, all_contigs, format),
        CatalogCommands::Export { output, catalog } => run_export(output, catalog),
        CatalogCommands::ListHierarchical { catalog } => {
            run_list_hierarchical(catalog, format, verbose)
        }
        CatalogCommands::BuildHierarchical {
            id,
            name,
            inputs,
            assembly_id,
            version_id,
            source,
            download_url,
            tags,
            output,
            append_to,
            force,
            require_md5,
            infer_assembly,
        } => run_build_hierarchical(
            id,
            name,
            inputs,
            assembly_id,
            version_id,
            source,
            download_url,
            tags,
            output,
            append_to,
            force,
            require_md5,
            infer_assembly,
            format,
            verbose,
        ),
        CatalogCommands::Build {
            id,
            name,
            inputs,
            assembly,
            source,
            description,
            download_url,
            assembly_report_url,
            tags,
            output,
            append_to,
            force,
            input_format,
            require_md5,
        } => run_build(
            id,
            name,
            inputs,
            assembly,
            source,
            description,
            download_url,
            assembly_report_url,
            tags,
            output,
            append_to,
            force,
            input_format,
            require_md5,
            format,
            verbose,
        ),
    }
}

fn run_list(
    catalog_path: Option<PathBuf>,
    assembly_filter: Option<&str>,
    source_filter: Option<&str>,
    format: OutputFormat,
    verbose: bool,
) -> anyhow::Result<()> {
    let catalog = if let Some(path) = catalog_path {
        ReferenceCatalog::load_from_file(&path)?
    } else {
        ReferenceCatalog::load_embedded()?
    };

    if verbose {
        eprintln!("Loaded catalog with {} references", catalog.len());
    }

    // Filter references
    let filtered: Vec<_> = catalog
        .references
        .iter()
        .filter(|r| {
            if let Some(assembly) = &assembly_filter {
                let ref_assembly = format!("{}", r.assembly).to_lowercase();
                if !ref_assembly.contains(&assembly.to_lowercase()) {
                    return false;
                }
            }
            if let Some(source) = &source_filter {
                let ref_source = format!("{}", r.source).to_lowercase();
                if !ref_source.contains(&source.to_lowercase()) {
                    return false;
                }
            }
            true
        })
        .collect();

    match format {
        OutputFormat::Text => {
            // Calculate column widths dynamically
            let id_width = filtered
                .iter()
                .map(|r| r.id.0.len())
                .max()
                .unwrap_or(2)
                .max(2);
            let name_width = filtered
                .iter()
                .map(|r| r.display_name.len().min(35))
                .max()
                .unwrap_or(4)
                .max(4);
            let assembly_width = filtered
                .iter()
                .map(|r| format!("{}", r.assembly).len())
                .max()
                .unwrap_or(8)
                .max(8);
            let source_width = filtered
                .iter()
                .map(|r| format!("{}", r.source).len())
                .max()
                .unwrap_or(6)
                .max(6);

            let total_width = id_width + name_width + assembly_width + source_width + 8 + 8;

            println!("Reference Catalog ({} references)\n", filtered.len());
            println!(
                "{:<id_w$} {:<name_w$} {:<asm_w$} {:<src_w$} {:>8}",
                "ID",
                "Name",
                "Assembly",
                "Source",
                "Contigs",
                id_w = id_width,
                name_w = name_width,
                asm_w = assembly_width,
                src_w = source_width
            );
            println!("{}", "-".repeat(total_width));

            for r in &filtered {
                println!(
                    "{:<id_w$} {:<name_w$} {:<asm_w$} {:<src_w$} {:>8}",
                    r.id.0,
                    truncate(&r.display_name, name_width),
                    format!("{}", r.assembly),
                    format!("{}", r.source),
                    r.contigs.len(),
                    id_w = id_width,
                    name_w = name_width,
                    asm_w = assembly_width,
                    src_w = source_width
                );
                if verbose {
                    let md5_count = r.contigs.iter().filter(|c| c.md5.is_some()).count();
                    let md5_pct = if r.contigs.is_empty() {
                        0.0
                    } else {
                        100.0 * count_to_f64(md5_count) / count_to_f64(r.contigs.len())
                    };
                    if let Some(url) = &r.download_url {
                        println!(
                            "  └─ MD5: {}/{} ({:.0}%)  URL: {}",
                            md5_count,
                            r.contigs.len(),
                            md5_pct,
                            url
                        );
                    } else {
                        println!(
                            "  └─ MD5: {}/{} ({:.0}%)",
                            md5_count,
                            r.contigs.len(),
                            md5_pct
                        );
                    }
                }
            }
        }
        OutputFormat::Json => {
            let output: Vec<serde_json::Value> = filtered
                .iter()
                .map(|r| {
                    let md5_count = r.contigs.iter().filter(|c| c.md5.is_some()).count();
                    let role_counts = r.role_counts();
                    let mut json = serde_json::json!({
                        "id": r.id.0,
                        "display_name": r.display_name,
                        "assembly": format!("{}", r.assembly),
                        "source": format!("{}", r.source),
                        "contig_count": r.contigs.len(),
                        "md5_count": md5_count,
                        "has_decoy": r.has_decoy(),
                        "has_alt": r.has_alt(),
                        "fasta_url": r.download_url,
                        "assembly_report_url": r.assembly_report_url,
                        "role_counts": {
                            "assembled_molecule": role_counts.assembled_molecule,
                            "alt_scaffold": role_counts.alt_scaffold,
                            "fix_patch": role_counts.fix_patch,
                            "novel_patch": role_counts.novel_patch,
                            "unlocalized_scaffold": role_counts.unlocalized_scaffold,
                            "unplaced_scaffold": role_counts.unplaced_scaffold,
                            "unknown": role_counts.unknown,
                        },
                        "tags": r.tags,
                    });
                    // Add contigs_missing_from_fasta only if non-empty
                    if !r.contigs_missing_from_fasta.is_empty() {
                        json["contigs_missing_from_fasta"] =
                            serde_json::json!(&r.contigs_missing_from_fasta);
                    }
                    json
                })
                .collect();
            println!("{}", serde_json::to_string_pretty(&output)?);
        }
        OutputFormat::Tsv => {
            println!("id\tdisplay_name\tassembly\tsource\tcontig_count\tmd5_count\thas_decoy\thas_alt\tdownload_url");
            for r in &filtered {
                let md5_count = r.contigs.iter().filter(|c| c.md5.is_some()).count();
                println!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    r.id.0,
                    r.display_name,
                    r.assembly,
                    r.source,
                    r.contigs.len(),
                    md5_count,
                    r.has_decoy(),
                    r.has_alt(),
                    r.download_url.as_deref().unwrap_or("")
                );
            }
        }
    }

    Ok(())
}

fn run_show(
    id: String,
    catalog_path: Option<PathBuf>,
    all_contigs: bool,
    format: OutputFormat,
) -> anyhow::Result<()> {
    let catalog = if let Some(path) = catalog_path {
        ReferenceCatalog::load_from_file(&path)?
    } else {
        ReferenceCatalog::load_embedded()?
    };

    let ref_id = crate::core::types::ReferenceId::new(&id);
    let reference = catalog
        .get(&ref_id)
        .ok_or_else(|| anyhow::anyhow!("Reference '{}' not found", id))?;

    match format {
        OutputFormat::Text => {
            println!("Reference: {}\n", reference.display_name);
            println!("ID:       {}", reference.id);
            println!("Assembly: {}", reference.assembly);
            println!("Source:   {}", reference.source);
            println!("Naming:   {:?}", reference.naming_convention);
            println!("Contigs:  {}", reference.contigs.len());
            println!("Has Decoy: {}", reference.has_decoy());
            println!("Has ALT:   {}", reference.has_alt());

            if let Some(desc) = &reference.description {
                println!("\nDescription: {}", desc);
            }

            if let Some(url) = &reference.download_url {
                println!("\nDownload URL: {}", url);
            }

            if !reference.tags.is_empty() {
                println!("\nTags: {}", reference.tags.join(", "));
            }

            let contigs_to_show = if all_contigs {
                &reference.contigs[..]
            } else {
                &reference.contigs[..reference.contigs.len().min(25)]
            };

            println!("\nContigs:");
            println!("{:<25} {:>15} MD5", "Name", "Length");
            println!("{}", "-".repeat(80));
            for contig in contigs_to_show {
                println!(
                    "{:<25} {:>15} {}",
                    contig.name,
                    contig.length,
                    contig.md5.as_deref().unwrap_or("-")
                );
            }

            if !all_contigs && reference.contigs.len() > 25 {
                println!(
                    "\n... and {} more contigs (use --all-contigs to show all)",
                    reference.contigs.len() - 25
                );
            }
        }
        OutputFormat::Json => {
            println!("{}", serde_json::to_string_pretty(&reference)?);
        }
        OutputFormat::Tsv => {
            println!("name\tlength\tmd5");
            for contig in &reference.contigs {
                println!(
                    "{}\t{}\t{}",
                    contig.name,
                    contig.length,
                    contig.md5.as_deref().unwrap_or("")
                );
            }
        }
    }

    Ok(())
}

fn run_export(output: PathBuf, catalog_path: Option<PathBuf>) -> anyhow::Result<()> {
    let catalog = if let Some(path) = catalog_path {
        ReferenceCatalog::load_from_file(&path)?
    } else {
        ReferenceCatalog::load_embedded()?
    };

    let json = catalog.to_json()?;
    std::fs::write(&output, json)?;

    println!(
        "Exported {} references to {}",
        catalog.len(),
        output.display()
    );

    Ok(())
}

fn run_list_hierarchical(
    catalog_path: PathBuf,
    format: OutputFormat,
    verbose: bool,
) -> anyhow::Result<()> {
    let catalog = HierarchicalCatalog::load(&catalog_path)?;

    if verbose {
        eprintln!(
            "Loaded hierarchical catalog v{} with {} assemblies",
            catalog.version,
            catalog.assemblies.len()
        );
    }

    match format {
        OutputFormat::Text => {
            println!("Hierarchical Reference Catalog (v{})\n", catalog.version);

            // Count totals
            let mut total_versions = 0;
            let mut total_distributions = 0;
            let mut total_contigs = 0;

            for assembly in &catalog.assemblies {
                total_versions += assembly.versions.len();
                for version in &assembly.versions {
                    total_distributions += version.fasta_distributions.len();
                    for dist in &version.fasta_distributions {
                        total_contigs += dist.contigs.len();
                    }
                }
            }

            println!(
                "Summary: {} assemblies, {} versions, {} distributions, {} total contigs\n",
                catalog.assemblies.len(),
                total_versions,
                total_distributions,
                total_contigs
            );

            // List hierarchy
            for assembly in &catalog.assemblies {
                println!("{} ({})", assembly.name, assembly.id);
                println!("  Organism: {}", assembly.organism);

                for version in &assembly.versions {
                    println!("\n  Version: {} ({})", version.version, version.id);

                    // Show report source
                    match &version.source {
                        crate::core::assembly::ReportSource::Ncbi { accession, .. } => {
                            println!("    Source: NCBI ({})", accession);
                        }
                        crate::core::assembly::ReportSource::DerivedFromFasta {
                            base_assembly,
                            ..
                        } => {
                            if let Some(base) = base_assembly {
                                println!("    Source: Derived from FASTA (base: {})", base);
                            } else {
                                println!("    Source: Derived from FASTA");
                            }
                        }
                        crate::core::assembly::ReportSource::Manual { .. } => {
                            println!("    Source: Manual");
                        }
                    }

                    if !version.report_contigs.is_empty() {
                        println!("    Report contigs: {}", version.report_contigs.len());
                    }

                    println!("    Distributions:");
                    for dist in &version.fasta_distributions {
                        let md5_count = dist.contigs.iter().filter(|c| !c.md5.is_empty()).count();
                        println!(
                            "      - {} ({}): {} contigs, {} with MD5",
                            dist.display_name,
                            dist.id,
                            dist.contigs.len(),
                            md5_count
                        );

                        if verbose {
                            // Show presence breakdown
                            let counts = dist.presence_counts();
                            if counts.in_both > 0 || counts.fasta_only > 0 {
                                println!(
                                    "        Presence: {} in-both, {} fasta-only",
                                    counts.in_both, counts.fasta_only
                                );
                            }

                            if let Some(url) = &dist.download_url {
                                println!("        URL: {}", url);
                            }
                        }
                    }
                }
                println!();
            }

            // Standalone distributions
            if !catalog.standalone_distributions.is_empty() {
                println!("Standalone Distributions:");
                for dist in &catalog.standalone_distributions {
                    let md5_count = dist.contigs.iter().filter(|c| !c.md5.is_empty()).count();
                    println!(
                        "  - {} ({}): {} contigs, {} with MD5",
                        dist.display_name,
                        dist.id,
                        dist.contigs.len(),
                        md5_count
                    );
                }
            }
        }
        OutputFormat::Json => {
            println!("{}", serde_json::to_string_pretty(&catalog)?);
        }
        OutputFormat::Tsv => {
            println!(
                "assembly_id\tversion_id\tdistribution_id\tdisplay_name\tcontig_count\tmd5_count"
            );
            for assembly in &catalog.assemblies {
                for version in &assembly.versions {
                    for dist in &version.fasta_distributions {
                        let md5_count = dist.contigs.iter().filter(|c| !c.md5.is_empty()).count();
                        println!(
                            "{}\t{}\t{}\t{}\t{}\t{}",
                            assembly.id,
                            version.id,
                            dist.id,
                            dist.display_name,
                            dist.contigs.len(),
                            md5_count
                        );
                    }
                }
            }
            // Standalone
            for dist in &catalog.standalone_distributions {
                let md5_count = dist.contigs.iter().filter(|c| !c.md5.is_empty()).count();
                println!(
                    "\t\t{}\t{}\t{}\t{}",
                    dist.id,
                    dist.display_name,
                    dist.contigs.len(),
                    md5_count
                );
            }
        }
    }

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn run_build_hierarchical(
    id: String,
    name: String,
    inputs: Vec<PathBuf>,
    assembly_id: Option<String>,
    version_id: Option<String>,
    source: Option<String>,
    download_url: Option<String>,
    tags: Option<String>,
    output: Option<PathBuf>,
    append_to: Option<PathBuf>,
    force: bool,
    require_md5: bool,
    infer_assembly: Option<Option<PathBuf>>,
    format: OutputFormat,
    verbose: bool,
) -> anyhow::Result<()> {
    use crate::catalog::builder::DistributionBuilder;

    // Parse source
    let ref_source = source
        .map(|s| parse_reference_source(&s))
        .unwrap_or(ReferenceSource::Custom("custom".to_string()));

    // Parse tags
    let tags: Vec<String> = tags
        .map(|s| s.split(',').map(|t| t.trim().to_string()).collect())
        .unwrap_or_default();

    // Create builder
    let mut builder = DistributionBuilder::new(&id)
        .with_display_name(&name)
        .with_source(ref_source);

    if let Some(url) = download_url {
        builder = builder.with_download_url(url);
    }
    if !tags.is_empty() {
        builder = builder.with_tags(tags);
    }

    // Process input files
    for input_path in &inputs {
        if !input_path.exists() {
            anyhow::bail!("Input file not found: {}", input_path.display());
        }

        if verbose {
            eprintln!("Processing: {}", input_path.display());
        }

        builder.add_input(input_path)?;
    }

    // Build the distribution
    let dist = builder.build()?;

    // Check MD5 requirement
    if require_md5 {
        let missing_md5: Vec<_> = dist
            .contigs
            .iter()
            .filter(|c| c.md5.is_empty())
            .map(|c| c.name.as_str())
            .collect();

        if !missing_md5.is_empty() {
            anyhow::bail!(
                "MD5 required but {} contig(s) lack MD5: {}",
                missing_md5.len(),
                missing_md5.join(", ")
            );
        }
    }

    // Summary
    let md5_count = dist.contigs.iter().filter(|c| !c.md5.is_empty()).count();
    if verbose {
        eprintln!(
            "Built distribution '{}' with {} contigs ({} with MD5)",
            id,
            dist.contigs.len(),
            md5_count
        );
    }

    // Inference of base assembly
    let (inferred_assembly_id, inferred_version_id) = if infer_assembly.is_some() {
        // Load catalog for inference
        let infer_catalog = match &infer_assembly {
            Some(Some(path)) => Some(HierarchicalCatalog::load(path)?),
            Some(None) => {
                // Try to use append_to catalog first, otherwise try embedded (which doesn't exist in hierarchical format)
                if let Some(ref append_path) = append_to {
                    Some(HierarchicalCatalog::load(append_path)?)
                } else {
                    if verbose {
                        eprintln!("Warning: No catalog specified for inference. Use --infer-assembly=<path> or --append-to");
                    }
                    None
                }
            }
            None => None,
        };

        if let Some(ref catalog) = infer_catalog {
            match catalog.infer_base_assembly_default(&dist.contigs) {
                Some(inferred) => {
                    if verbose {
                        eprintln!(
                            "Inferred base assembly: {} {} ({:.1}% match, {}/{} contigs)",
                            inferred.assembly_name,
                            inferred.version_string,
                            inferred.match_rate * 100.0,
                            inferred.matched_contigs,
                            inferred.total_input_contigs
                        );
                    }
                    (
                        assembly_id.clone().or(Some(inferred.assembly_id)),
                        version_id.clone().or(Some(inferred.version_id)),
                    )
                }
                None => {
                    if verbose {
                        eprintln!("Could not infer base assembly (no match above 90% threshold)");
                    }
                    (assembly_id.clone(), version_id.clone())
                }
            }
        } else {
            (assembly_id.clone(), version_id.clone())
        }
    } else {
        (assembly_id.clone(), version_id.clone())
    };

    // Output handling
    if let Some(append_path) = append_to {
        // Append to existing catalog
        let mut catalog = HierarchicalCatalog::load(&append_path)?;

        // Check if we need to add to a specific assembly/version (using inferred if available)
        if let (Some(asm_id), Some(ver_id)) = (&inferred_assembly_id, &inferred_version_id) {
            let mut found = false;
            for assembly in &mut catalog.assemblies {
                if assembly.id == *asm_id {
                    for version in &mut assembly.versions {
                        if version.id == *ver_id {
                            // Check for existing distribution
                            if !force && version.fasta_distributions.iter().any(|d| d.id == id) {
                                anyhow::bail!(
                                    "Distribution '{}' already exists in version '{}'. Use --force to overwrite.",
                                    id,
                                    ver_id
                                );
                            }

                            // Remove existing if force
                            version.fasta_distributions.retain(|d| d.id != id);
                            version.fasta_distributions.push(dist.clone());
                            found = true;
                            break;
                        }
                    }
                }
            }
            if !found {
                anyhow::bail!(
                    "Assembly '{}' with version '{}' not found in catalog",
                    asm_id,
                    ver_id
                );
            }
        } else {
            // Add as standalone distribution
            if !force && catalog.standalone_distributions.iter().any(|d| d.id == id) {
                anyhow::bail!(
                    "Standalone distribution '{}' already exists. Use --force to overwrite.",
                    id
                );
            }
            catalog.standalone_distributions.retain(|d| d.id != id);
            catalog.standalone_distributions.push(dist.clone());
        }

        catalog.save(&append_path)?;
        eprintln!("Added distribution '{}' to {}", id, append_path.display());
    } else if let Some(out_path) = output {
        // Create new output
        if out_path.exists() && !force {
            anyhow::bail!(
                "Output file '{}' exists. Use --force to overwrite.",
                out_path.display()
            );
        }

        // Output as standalone distribution JSON or wrap in catalog
        match format {
            OutputFormat::Json => {
                // Just output the distribution as JSON
                let json = serde_json::to_string_pretty(&dist)?;
                std::fs::write(&out_path, json)?;
                eprintln!("Wrote distribution to {}", out_path.display());
            }
            _ => {
                // Create a catalog with just this distribution
                let catalog = HierarchicalCatalog::new().with_standalone_distribution(dist);
                catalog.save(&out_path)?;
                eprintln!("Wrote hierarchical catalog to {}", out_path.display());
            }
        }
    } else {
        // Print to stdout
        match format {
            OutputFormat::Json => {
                println!("{}", serde_json::to_string_pretty(&dist)?);
            }
            OutputFormat::Text => {
                print_distribution_summary(&dist);
            }
            OutputFormat::Tsv => {
                println!("name\tlength\tmd5\treport_contig_id");
                for c in &dist.contigs {
                    println!(
                        "{}\t{}\t{}\t{}",
                        c.name,
                        c.length,
                        c.md5,
                        c.report_contig_id
                            .map(|i| i.to_string())
                            .unwrap_or_default()
                    );
                }
            }
        }
    }

    Ok(())
}

fn print_distribution_summary(dist: &crate::core::assembly::FastaDistribution) {
    println!("Distribution: {} ({})", dist.display_name, dist.id);
    println!("Source: {:?}", dist.source);
    if let Some(url) = &dist.download_url {
        println!("Download URL: {}", url);
    }
    if !dist.tags.is_empty() {
        println!("Tags: {}", dist.tags.join(", "));
    }
    println!("Contigs: {}", dist.contigs.len());

    let md5_count = dist.contigs.iter().filter(|c| !c.md5.is_empty()).count();
    println!("With MD5: {}", md5_count);

    let linked = dist
        .contigs
        .iter()
        .filter(|c| c.report_contig_id.is_some())
        .count();
    println!("Linked to report: {}", linked);

    // Show presence counts
    let counts = dist.presence_counts();
    if counts.in_both > 0 || counts.fasta_only > 0 {
        println!(
            "Presence: {} in-both, {} fasta-only",
            counts.in_both, counts.fasta_only
        );
    }
}

fn parse_reference_source(s: &str) -> ReferenceSource {
    match s.to_lowercase().as_str() {
        "ucsc" => ReferenceSource::Ucsc,
        "ncbi" => ReferenceSource::Ncbi,
        "broad" => ReferenceSource::Broad,
        "ensembl" => ReferenceSource::Ensembl,
        "1kg" | "1000genomes" => ReferenceSource::OneThousandGenomes,
        "dragen" | "illumina" => ReferenceSource::Illumina,
        _ => ReferenceSource::Custom(s.to_string()),
    }
}

fn truncate(s: &str, max_len: usize) -> String {
    if s.len() <= max_len {
        s.to_string()
    } else {
        format!("{}...", &s[..max_len - 3])
    }
}

#[allow(clippy::too_many_arguments)]
fn run_build(
    id: String,
    name: String,
    inputs: Vec<PathBuf>,
    assembly: Option<String>,
    source: Option<String>,
    description: Option<String>,
    download_url: Option<String>,
    assembly_report_url: Option<String>,
    tags: Option<String>,
    output: Option<PathBuf>,
    append_to: Option<PathBuf>,
    force: bool,
    input_format: Option<InputFormatArg>,
    require_md5: bool,
    format: OutputFormat,
    verbose: bool,
) -> anyhow::Result<()> {
    // Parse assembly
    let assembly = assembly.map(|s| parse_assembly(&s));

    // Parse source
    let source = source.map(|s| parse_source(&s));

    // Parse tags
    let tags: Vec<String> = tags
        .map(|s| s.split(',').map(|t| t.trim().to_string()).collect())
        .unwrap_or_default();

    // Create builder
    let mut builder = ReferenceBuilder::new(&id, &name);

    if let Some(assembly) = assembly {
        builder = builder.assembly(assembly);
    }
    if let Some(source) = source {
        builder = builder.source(source);
    }
    if let Some(desc) = description {
        builder = builder.description(desc);
    }
    if let Some(url) = download_url {
        builder = builder.download_url(url);
    }
    if let Some(url) = assembly_report_url {
        builder = builder.assembly_report_url(url);
    }
    if !tags.is_empty() {
        builder = builder.tags(tags);
    }

    // Process input files
    for input_path in &inputs {
        if !input_path.exists() {
            anyhow::bail!("Input file not found: {}", input_path.display());
        }

        if verbose {
            eprintln!("Processing: {}", input_path.display());
        }

        if let Some(fmt) = input_format {
            builder.add_input_with_format(input_path, fmt.into())?;
        } else {
            builder.add_input(input_path)?;
        }
    }

    // Get summary before building
    let summary = builder.summary();

    // Check for conflicts
    if !summary.conflicts.is_empty() {
        eprintln!("Build failed due to conflicts:");
        for conflict in &summary.conflicts {
            eprintln!("  - {}", conflict);
        }
        anyhow::bail!(
            "Build failed: {} conflict(s) detected",
            summary.conflicts.len()
        );
    }

    // Check MD5 requirement
    if require_md5 && summary.with_md5 < summary.total_contigs {
        anyhow::bail!(
            "Build failed: --require-md5 specified but only {}/{} contigs have MD5",
            summary.with_md5,
            summary.total_contigs
        );
    }

    // Build the reference
    let reference = builder.build()?;

    // Print summary
    if verbose || matches!(format, OutputFormat::Text) {
        eprintln!("{}", summary);
    }

    // Handle output
    if let Some(catalog_path) = append_to {
        // Append to existing catalog
        let mut catalog = if catalog_path.exists() {
            ReferenceCatalog::load_from_file(&catalog_path)?
        } else {
            ReferenceCatalog::new()
        };

        // Check if ID already exists
        let ref_id = crate::core::types::ReferenceId::new(&id);
        if catalog.get(&ref_id).is_some() {
            if force {
                eprintln!("Warning: Overwriting existing reference '{}'", id);
                // Remove old reference by rebuilding catalog without it
                let refs: Vec<_> = catalog
                    .references
                    .into_iter()
                    .filter(|r| r.id != ref_id)
                    .collect();
                catalog = ReferenceCatalog::new();
                for r in refs {
                    catalog.add_reference(r);
                }
            } else {
                anyhow::bail!(
                    "Reference '{}' already exists in catalog. Use --force to overwrite.",
                    id
                );
            }
        }

        catalog.add_reference(reference);
        let json = catalog.to_json()?;
        std::fs::write(&catalog_path, json)?;

        println!(
            "Added reference '{}' to {} ({} total references)",
            id,
            catalog_path.display(),
            catalog.len()
        );
    } else if let Some(output_path) = output {
        // Write single reference to file
        let json = serde_json::to_string_pretty(&reference)?;
        std::fs::write(&output_path, &json)?;
        println!("Wrote reference '{}' to {}", id, output_path.display());
    } else {
        // Print to stdout
        match format {
            OutputFormat::Json => {
                println!("{}", serde_json::to_string_pretty(&reference)?);
            }
            OutputFormat::Text | OutputFormat::Tsv => {
                // Print summary info
                println!("Reference: {}", reference.display_name);
                println!("ID:        {}", reference.id);
                println!("Assembly:  {}", reference.assembly);
                println!("Source:    {}", reference.source);
                println!("Contigs:   {}", reference.contigs.len());
                println!();
                println!("Use --output <file> to save as JSON");
            }
        }
    }

    Ok(())
}

fn parse_assembly(s: &str) -> Assembly {
    let lower = s.to_lowercase();
    match lower.as_str() {
        "grch37" | "hg19" | "b37" => Assembly::Grch37,
        "grch38" | "hg38" => Assembly::Grch38,
        _ => Assembly::Other(s.to_string()),
    }
}

fn parse_source(s: &str) -> ReferenceSource {
    let lower = s.to_lowercase();
    match lower.as_str() {
        "ucsc" => ReferenceSource::Ucsc,
        "ncbi" | "grc" => ReferenceSource::Ncbi,
        "broad" => ReferenceSource::Broad,
        "ensembl" => ReferenceSource::Ensembl,
        "illumina" | "dragen" => ReferenceSource::Illumina,
        "1kg" | "1000genomes" => ReferenceSource::OneThousandGenomes,
        _ => ReferenceSource::Custom(s.to_string()),
    }
}
