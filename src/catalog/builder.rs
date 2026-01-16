//! Reference builder for creating catalog entries from multiple input sources.
//!
//! The `ReferenceBuilder` collates metadata from multiple input files (dict, FAI,
//! NCBI assembly reports, VCF headers, etc.) and produces a `KnownReference` with
//! all contigs and aliases properly merged.

use std::collections::{HashMap, HashSet};
use std::path::Path;
use thiserror::Error;

use crate::core::assembly::{ContigMergeError, FastaContig, FastaDistribution};
use crate::core::contig::{Contig, SequenceRole};
use crate::core::reference::KnownReference;
use crate::core::types::{Assembly, ReferenceSource};

#[derive(Error, Debug)]
pub enum BuilderError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Parse error: {0}")]
    Parse(String),

    #[error("Conflict: {0}")]
    Conflict(String),

    #[error("Missing required field: {0}")]
    MissingField(String),

    #[error("Merge error: {0}")]
    Merge(#[from] ContigMergeError),
}

/// Input format for auto-detection
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputFormat {
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

impl InputFormat {
    /// Detect format from file extension
    #[must_use]
    pub fn from_path(path: &Path) -> Option<Self> {
        let name = path.file_name()?.to_str()?;
        let name_lower = name.to_lowercase();

        // Check for NCBI assembly report pattern first
        // name_lower is already lowercase, so case-sensitive check is fine
        #[allow(clippy::case_sensitive_file_extension_comparisons)]
        if name_lower.contains("_assembly_report") && name_lower.ends_with(".txt") {
            return Some(Self::NcbiReport);
        }

        // Check for gzipped FASTA files (.fa.gz, .fasta.gz, .fna.gz)
        if name_lower.ends_with(".fa.gz")
            || name_lower.ends_with(".fasta.gz")
            || name_lower.ends_with(".fna.gz")
            || name_lower.ends_with(".fa.bgz")
            || name_lower.ends_with(".fasta.bgz")
            || name_lower.ends_with(".fna.bgz")
        {
            return Some(Self::Fasta);
        }

        let ext = path.extension()?.to_str()?.to_lowercase();
        match ext.as_str() {
            "dict" => Some(Self::Dict),
            "fai" => Some(Self::Fai),
            "fa" | "fasta" | "fna" => Some(Self::Fasta),
            "sam" => Some(Self::Sam),
            "bam" => Some(Self::Bam),
            "cram" => Some(Self::Cram),
            "vcf" => Some(Self::Vcf),
            "tsv" | "txt" => Some(Self::Tsv),
            "gz" => {
                // Check for .vcf.gz (stem is lowercased)
                let stem = path.file_stem()?.to_str()?.to_lowercase();
                #[allow(clippy::case_sensitive_file_extension_comparisons)]
                if stem.ends_with(".vcf") {
                    Some(Self::Vcf)
                } else {
                    None
                }
            }
            _ => None,
        }
    }
}

/// Metadata for a single contig, collected from multiple sources
#[derive(Debug, Clone)]
pub struct ContigMetadata {
    /// Primary name (exact, from first source)
    pub primary_name: String,

    /// Length (must be consistent across sources)
    pub length: Option<u64>,

    /// MD5 checksum
    pub md5: Option<String>,

    /// Explicit aliases (from AN tag or NCBI report)
    pub aliases: HashSet<String>,

    /// Assembly tag (AS)
    pub assembly: Option<String>,

    /// URI (UR)
    pub uri: Option<String>,

    /// Species (SP)
    pub species: Option<String>,

    /// Sequence role from NCBI assembly report
    pub sequence_role: SequenceRole,

    /// Sources that contributed to this contig
    pub sources: Vec<String>,
}

impl ContigMetadata {
    fn new(name: String) -> Self {
        Self {
            primary_name: name,
            length: None,
            md5: None,
            aliases: HashSet::new(),
            assembly: None,
            uri: None,
            species: None,
            sequence_role: SequenceRole::Unknown,
            sources: Vec::new(),
        }
    }

    /// Convert to a Contig
    fn to_contig(&self) -> Option<Contig> {
        let length = self.length?;
        let mut contig = Contig::new(&self.primary_name, length);
        contig.md5.clone_from(&self.md5);
        contig.assembly.clone_from(&self.assembly);
        contig.uri.clone_from(&self.uri);
        contig.species.clone_from(&self.species);
        contig.aliases = self.aliases.iter().cloned().collect();
        contig.sequence_role = self.sequence_role;
        Some(contig)
    }
}

/// Record of a processed input file
#[derive(Debug, Clone)]
pub struct InputRecord {
    pub path: String,
    pub format: InputFormat,
    pub contigs_found: usize,
    pub contigs_merged: usize,
    pub aliases_added: usize,
}

/// Builder that collates metadata from multiple input sources
pub struct ReferenceBuilder {
    id: String,
    display_name: String,
    assembly: Option<Assembly>,
    source: Option<ReferenceSource>,
    description: Option<String>,
    download_url: Option<String>,
    assembly_report_url: Option<String>,
    tags: Vec<String>,

    /// Contigs keyed by EXACT primary name
    contigs: HashMap<String, ContigMetadata>,

    /// Preserve original contig order (first seen)
    contig_order: Vec<String>,

    /// Reverse lookup: alias -> primary name
    alias_to_primary: HashMap<String, String>,

    /// Records of processed inputs
    inputs_processed: Vec<InputRecord>,

    /// Conflicts detected
    conflicts: Vec<String>,

    /// Warnings
    warnings: Vec<String>,
}

impl ReferenceBuilder {
    /// Create a new builder with required fields
    pub fn new(id: impl Into<String>, display_name: impl Into<String>) -> Self {
        Self {
            id: id.into(),
            display_name: display_name.into(),
            assembly: None,
            source: None,
            description: None,
            download_url: None,
            assembly_report_url: None,
            tags: Vec::new(),
            contigs: HashMap::new(),
            contig_order: Vec::new(),
            alias_to_primary: HashMap::new(),
            inputs_processed: Vec::new(),
            conflicts: Vec::new(),
            warnings: Vec::new(),
        }
    }

    #[must_use]
    pub fn assembly(mut self, assembly: Assembly) -> Self {
        self.assembly = Some(assembly);
        self
    }

    #[must_use]
    pub fn source(mut self, source: ReferenceSource) -> Self {
        self.source = Some(source);
        self
    }

    #[must_use]
    pub fn description(mut self, description: impl Into<String>) -> Self {
        self.description = Some(description.into());
        self
    }

    #[must_use]
    pub fn download_url(mut self, url: impl Into<String>) -> Self {
        self.download_url = Some(url.into());
        self
    }

    #[must_use]
    pub fn assembly_report_url(mut self, url: impl Into<String>) -> Self {
        self.assembly_report_url = Some(url.into());
        self
    }

    #[must_use]
    pub fn tags(mut self, tags: Vec<String>) -> Self {
        self.tags = tags;
        self
    }

    /// Add an input file (auto-detect format)
    ///
    /// # Errors
    ///
    /// Returns `BuilderError::Parse` if format cannot be detected, or other
    /// errors from parsing the specific format.
    pub fn add_input(&mut self, path: &Path) -> Result<(), BuilderError> {
        let format = InputFormat::from_path(path).ok_or_else(|| {
            BuilderError::Parse(format!("Cannot detect format for: {}", path.display()))
        })?;
        self.add_input_with_format(path, format)
    }

    /// Add an input file with explicit format
    ///
    /// # Errors
    ///
    /// Returns `BuilderError::Io` if the file cannot be read, `BuilderError::Parse`
    /// if parsing fails, or `BuilderError::Conflict` if contig data conflicts.
    pub fn add_input_with_format(
        &mut self,
        path: &Path,
        format: InputFormat,
    ) -> Result<(), BuilderError> {
        let path_str = path.display().to_string();

        match format {
            InputFormat::Dict | InputFormat::Sam => {
                self.add_dict_or_sam(path, &path_str, format)?;
            }
            InputFormat::Bam => {
                self.add_bam(path, &path_str)?;
            }
            InputFormat::Cram => {
                self.add_cram(path, &path_str)?;
            }
            InputFormat::Fai => {
                self.add_fai(path, &path_str)?;
            }
            InputFormat::Fasta => {
                self.add_fasta(path, &path_str)?;
            }
            InputFormat::NcbiReport => {
                self.add_ncbi_report(path, &path_str)?;
            }
            InputFormat::Vcf => {
                self.add_vcf(path, &path_str)?;
            }
            InputFormat::Tsv => {
                self.add_tsv(path, &path_str)?;
            }
        }

        Ok(())
    }

    fn add_dict_or_sam(
        &mut self,
        path: &Path,
        path_str: &str,
        format: InputFormat,
    ) -> Result<(), BuilderError> {
        let content = std::fs::read_to_string(path)?;
        let query = crate::parsing::sam::parse_header_text(&content)
            .map_err(|e| BuilderError::Parse(e.to_string()))?;

        let mut record = InputRecord {
            path: path_str.to_string(),
            format,
            contigs_found: query.contigs.len(),
            contigs_merged: 0,
            aliases_added: 0,
        };

        for contig in query.contigs {
            let (merged, aliases) = self.merge_contig(&contig, path_str)?;
            if merged {
                record.contigs_merged += 1;
            }
            record.aliases_added += aliases;
        }

        self.inputs_processed.push(record);
        Ok(())
    }

    fn add_bam(&mut self, path: &Path, path_str: &str) -> Result<(), BuilderError> {
        let query = crate::parsing::sam::parse_file(path)
            .map_err(|e| BuilderError::Parse(e.to_string()))?;

        let mut record = InputRecord {
            path: path_str.to_string(),
            format: InputFormat::Bam,
            contigs_found: query.contigs.len(),
            contigs_merged: 0,
            aliases_added: 0,
        };

        for contig in query.contigs {
            let (merged, aliases) = self.merge_contig(&contig, path_str)?;
            if merged {
                record.contigs_merged += 1;
            }
            record.aliases_added += aliases;
        }

        self.inputs_processed.push(record);
        Ok(())
    }

    fn add_cram(&mut self, path: &Path, path_str: &str) -> Result<(), BuilderError> {
        let query = crate::parsing::sam::parse_file(path)
            .map_err(|e| BuilderError::Parse(e.to_string()))?;

        let mut record = InputRecord {
            path: path_str.to_string(),
            format: InputFormat::Cram,
            contigs_found: query.contigs.len(),
            contigs_merged: 0,
            aliases_added: 0,
        };

        for contig in query.contigs {
            let (merged, aliases) = self.merge_contig(&contig, path_str)?;
            if merged {
                record.contigs_merged += 1;
            }
            record.aliases_added += aliases;
        }

        self.inputs_processed.push(record);
        Ok(())
    }

    fn add_fai(&mut self, path: &Path, path_str: &str) -> Result<(), BuilderError> {
        let content = std::fs::read_to_string(path)?;
        let query = crate::parsing::fai::parse_fai_text(&content)
            .map_err(|e| BuilderError::Parse(e.to_string()))?;

        let mut record = InputRecord {
            path: path_str.to_string(),
            format: InputFormat::Fai,
            contigs_found: query.contigs.len(),
            contigs_merged: 0,
            aliases_added: 0,
        };

        for contig in query.contigs {
            let (merged, aliases) = self.merge_contig(&contig, path_str)?;
            if merged {
                record.contigs_merged += 1;
            }
            record.aliases_added += aliases;
        }

        self.inputs_processed.push(record);
        Ok(())
    }

    fn add_ncbi_report(&mut self, path: &Path, path_str: &str) -> Result<(), BuilderError> {
        let content = std::fs::read_to_string(path)?;
        let entries = crate::parsing::ncbi_report::parse_ncbi_report_text(&content)
            .map_err(|e| BuilderError::Parse(e.to_string()))?;

        let mut record = InputRecord {
            path: path_str.to_string(),
            format: InputFormat::NcbiReport,
            contigs_found: entries.len(),
            contigs_merged: 0,
            aliases_added: 0,
        };

        for entry in entries {
            let contig = entry.to_contig();
            let (merged, aliases) = self.merge_contig(&contig, path_str)?;
            if merged {
                record.contigs_merged += 1;
            }
            record.aliases_added += aliases;
        }

        self.inputs_processed.push(record);
        Ok(())
    }

    fn add_vcf(&mut self, path: &Path, path_str: &str) -> Result<(), BuilderError> {
        let content = std::fs::read_to_string(path)?;
        let query = crate::parsing::vcf::parse_vcf_header_text(&content)
            .map_err(|e| BuilderError::Parse(e.to_string()))?;

        let mut record = InputRecord {
            path: path_str.to_string(),
            format: InputFormat::Vcf,
            contigs_found: query.contigs.len(),
            contigs_merged: 0,
            aliases_added: 0,
        };

        for contig in query.contigs {
            let (merged, aliases) = self.merge_contig(&contig, path_str)?;
            if merged {
                record.contigs_merged += 1;
            }
            record.aliases_added += aliases;
        }

        self.inputs_processed.push(record);
        Ok(())
    }

    fn add_tsv(&mut self, path: &Path, path_str: &str) -> Result<(), BuilderError> {
        let content = std::fs::read_to_string(path)?;
        let query = crate::parsing::tsv::parse_tsv_text(&content, '\t')
            .map_err(|e| BuilderError::Parse(e.to_string()))?;

        let mut record = InputRecord {
            path: path_str.to_string(),
            format: InputFormat::Tsv,
            contigs_found: query.contigs.len(),
            contigs_merged: 0,
            aliases_added: 0,
        };

        for contig in query.contigs {
            let (merged, aliases) = self.merge_contig(&contig, path_str)?;
            if merged {
                record.contigs_merged += 1;
            }
            record.aliases_added += aliases;
        }

        self.inputs_processed.push(record);
        Ok(())
    }

    fn add_fasta(&mut self, path: &Path, path_str: &str) -> Result<(), BuilderError> {
        let query = crate::parsing::fasta::parse_fasta_file_with_md5(path)
            .map_err(|e| BuilderError::Parse(e.to_string()))?;

        let mut record = InputRecord {
            path: path_str.to_string(),
            format: InputFormat::Fasta,
            contigs_found: query.contigs.len(),
            contigs_merged: 0,
            aliases_added: 0,
        };

        for contig in query.contigs {
            let (merged, aliases) = self.merge_contig(&contig, path_str)?;
            if merged {
                record.contigs_merged += 1;
            }
            record.aliases_added += aliases;
        }

        self.inputs_processed.push(record);
        Ok(())
    }

    /// Merge a contig into the builder.
    /// Returns (`was_merged_into_existing`, `num_aliases_added`)
    fn merge_contig(
        &mut self,
        contig: &Contig,
        source: &str,
    ) -> Result<(bool, usize), BuilderError> {
        // Find existing entry by exact name or alias
        let mut existing_primary = self.find_existing_primary(&contig.name);

        // Also check if any of the incoming contig's aliases match existing entries
        if existing_primary.is_none() {
            for alias in &contig.aliases {
                if let Some(primary) = self.find_existing_primary(alias) {
                    existing_primary = Some(primary);
                    break;
                }
            }
        }

        if let Some(primary) = existing_primary {
            // Merge into existing
            let metadata = self.contigs.get_mut(&primary).unwrap();
            let aliases_before = metadata.aliases.len();

            // Check length conflict
            if let Some(existing_len) = metadata.length {
                if existing_len != contig.length {
                    let msg = format!(
                        "Length conflict for '{}': {} vs {} (from {})",
                        contig.name, existing_len, contig.length, source
                    );
                    self.conflicts.push(msg.clone());
                    return Err(BuilderError::Conflict(msg));
                }
            } else {
                metadata.length = Some(contig.length);
            }

            // Check MD5 conflict
            if let (Some(existing_md5), Some(new_md5)) = (&metadata.md5, &contig.md5) {
                if existing_md5 != new_md5 {
                    let msg = format!(
                        "MD5 conflict for '{}': {} vs {} (from {})",
                        contig.name, existing_md5, new_md5, source
                    );
                    self.conflicts.push(msg.clone());
                    return Err(BuilderError::Conflict(msg));
                }
            } else if metadata.md5.is_none() && contig.md5.is_some() {
                metadata.md5.clone_from(&contig.md5);
            }

            // Merge aliases
            for alias in &contig.aliases {
                if !metadata.aliases.contains(alias) && alias != &metadata.primary_name {
                    // Check alias doesn't conflict with another contig
                    if let Some(other_primary) = self.alias_to_primary.get(alias) {
                        if other_primary != &primary {
                            self.warnings.push(format!(
                                "Alias '{alias}' already mapped to '{other_primary}', skipping for '{primary}'"
                            ));
                            continue;
                        }
                    }
                    metadata.aliases.insert(alias.clone());
                    self.alias_to_primary.insert(alias.clone(), primary.clone());
                }
            }

            // Also add the contig's name as an alias if different from primary
            if contig.name != primary && !metadata.aliases.contains(&contig.name) {
                metadata.aliases.insert(contig.name.clone());
                self.alias_to_primary
                    .insert(contig.name.clone(), primary.clone());
            }

            // Fill other fields
            if metadata.assembly.is_none() && contig.assembly.is_some() {
                metadata.assembly.clone_from(&contig.assembly);
            }
            if metadata.uri.is_none() && contig.uri.is_some() {
                metadata.uri.clone_from(&contig.uri);
            }
            if metadata.species.is_none() && contig.species.is_some() {
                metadata.species.clone_from(&contig.species);
            }
            // Update sequence role if not already set
            if matches!(metadata.sequence_role, SequenceRole::Unknown) {
                metadata.sequence_role = contig.sequence_role;
            }

            metadata.sources.push(source.to_string());

            let aliases_added = metadata.aliases.len() - aliases_before;
            Ok((true, aliases_added))
        } else {
            // Create new entry
            let mut metadata = ContigMetadata::new(contig.name.clone());
            metadata.length = Some(contig.length);
            metadata.md5.clone_from(&contig.md5);
            metadata.assembly.clone_from(&contig.assembly);
            metadata.uri.clone_from(&contig.uri);
            metadata.species.clone_from(&contig.species);
            metadata.sequence_role = contig.sequence_role;
            metadata.sources.push(source.to_string());

            // Add aliases
            let mut aliases_added = 0;
            for alias in &contig.aliases {
                if alias != &contig.name {
                    if let Some(other_primary) = self.alias_to_primary.get(alias) {
                        self.warnings.push(format!(
                            "Alias '{}' already mapped to '{}', skipping for '{}'",
                            alias, other_primary, contig.name
                        ));
                        continue;
                    }
                    metadata.aliases.insert(alias.clone());
                    self.alias_to_primary
                        .insert(alias.clone(), contig.name.clone());
                    aliases_added += 1;
                }
            }

            self.contig_order.push(contig.name.clone());
            self.contigs.insert(contig.name.clone(), metadata);

            Ok((false, aliases_added))
        }
    }

    /// Find the primary name for a contig name (checks exact match and aliases)
    fn find_existing_primary(&self, name: &str) -> Option<String> {
        // Check exact name match
        if self.contigs.contains_key(name) {
            return Some(name.to_string());
        }

        // Check alias match
        if let Some(primary) = self.alias_to_primary.get(name) {
            return Some(primary.clone());
        }

        None
    }

    /// Build the final `KnownReference`
    ///
    /// # Errors
    ///
    /// Returns `BuilderError::MissingField` if no contigs were added or required
    /// fields are missing.
    pub fn build(self) -> Result<KnownReference, BuilderError> {
        // Validate
        if self.contigs.is_empty() {
            return Err(BuilderError::MissingField("No contigs added".to_string()));
        }

        // Check for missing lengths
        let mut missing_length = Vec::new();
        for name in &self.contig_order {
            if let Some(metadata) = self.contigs.get(name) {
                if metadata.length.is_none() {
                    missing_length.push(name.clone());
                }
            }
        }
        if !missing_length.is_empty() {
            return Err(BuilderError::MissingField(format!(
                "Missing length for contigs: {missing_length:?}"
            )));
        }

        // Build contigs in order
        let contigs: Vec<Contig> = self
            .contig_order
            .iter()
            .filter_map(|name| self.contigs.get(name))
            .filter_map(ContigMetadata::to_contig)
            .collect();

        // Find contigs that are in assembly report but missing MD5 (i.e., not in FASTA)
        // These are typically assembled-molecule contigs that use external references (like MT in CHM13)
        let mut contigs_missing_from_fasta: Vec<String> = Vec::new();
        for (name, meta) in &self.contigs {
            if meta.md5.is_none() && matches!(meta.sequence_role, SequenceRole::AssembledMolecule) {
                contigs_missing_from_fasta.push(name.clone());
            }
        }
        contigs_missing_from_fasta.sort();

        // Determine assembly if not set
        let assembly = self
            .assembly
            .unwrap_or_else(|| detect_assembly_from_name(&self.display_name));

        // Determine source if not set
        let source = self
            .source
            .unwrap_or(ReferenceSource::Custom("Unknown".to_string()));

        // Determine naming convention
        let naming_convention = crate::core::contig::detect_naming_convention(&contigs);

        let mut reference = KnownReference {
            id: crate::core::types::ReferenceId::new(&self.id),
            display_name: self.display_name,
            assembly,
            source,
            naming_convention,
            download_url: self.download_url,
            assembly_report_url: self.assembly_report_url,
            contigs,
            description: self.description,
            tags: self.tags,
            contigs_missing_from_fasta,
            md5_set: HashSet::new(),
            name_length_set: HashSet::new(),
            signature: None,
        };

        reference.rebuild_indexes();
        Ok(reference)
    }

    /// Get summary of build
    #[must_use]
    pub fn summary(&self) -> BuildSummary {
        let total_contigs = self.contigs.len();
        let with_length = self.contigs.values().filter(|m| m.length.is_some()).count();
        let with_md5 = self.contigs.values().filter(|m| m.md5.is_some()).count();
        let with_aliases = self
            .contigs
            .values()
            .filter(|m| !m.aliases.is_empty())
            .count();

        // Collect assembled-molecule contigs missing MD5 (these are important)
        let mut missing_md5_assembled: Vec<String> = Vec::new();
        for (name, meta) in &self.contigs {
            if meta.md5.is_none() {
                // Check if it's an assembled-molecule (primary chromosome)
                if matches!(meta.sequence_role, SequenceRole::AssembledMolecule) {
                    missing_md5_assembled.push(name.clone());
                }
            }
        }

        // Count primary chromosomes
        let primary_count = self
            .contigs
            .keys()
            .filter(|name| {
                let contig = Contig::new(name.as_str(), 0);
                contig.is_primary_chromosome()
            })
            .count();

        // Count ALT contigs
        let alt_count = self
            .contigs
            .keys()
            .filter(|name| name.ends_with("_alt") || name.contains("_alt_"))
            .count();

        BuildSummary {
            id: self.id.clone(),
            display_name: self.display_name.clone(),
            assembly: self.assembly.clone(),
            source: self.source.clone(),
            inputs: self.inputs_processed.clone(),
            total_contigs,
            primary_chromosomes: primary_count,
            alt_contigs: alt_count,
            other_contigs: total_contigs.saturating_sub(primary_count + alt_count),
            with_length,
            with_md5,
            with_aliases,
            missing_md5_assembled,
            conflicts: self.conflicts.clone(),
            warnings: self.warnings.clone(),
        }
    }
}

/// Summary of the build process
#[derive(Debug, Clone)]
pub struct BuildSummary {
    pub id: String,
    pub display_name: String,
    pub assembly: Option<Assembly>,
    pub source: Option<ReferenceSource>,
    pub inputs: Vec<InputRecord>,
    pub total_contigs: usize,
    pub primary_chromosomes: usize,
    pub alt_contigs: usize,
    pub other_contigs: usize,
    pub with_length: usize,
    pub with_md5: usize,
    pub with_aliases: usize,
    pub missing_md5_assembled: Vec<String>,
    pub conflicts: Vec<String>,
    pub warnings: Vec<String>,
}

impl std::fmt::Display for BuildSummary {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Reference Builder Summary")?;
        writeln!(f, "=========================")?;
        writeln!(f, "ID:       {}", self.id)?;
        writeln!(f, "Name:     {}", self.display_name)?;
        if let Some(ref assembly) = self.assembly {
            writeln!(f, "Assembly: {assembly:?}")?;
        }
        if let Some(ref source) = self.source {
            writeln!(f, "Source:   {source:?}")?;
        }
        writeln!(f)?;

        writeln!(f, "Inputs:")?;
        for (i, input) in self.inputs.iter().enumerate() {
            writeln!(
                f,
                "  [{}] {} ({:?}) -> {} contigs, {} merged, {} aliases",
                i + 1,
                input.path,
                input.format,
                input.contigs_found,
                input.contigs_merged,
                input.aliases_added
            )?;
        }
        writeln!(f)?;

        writeln!(f, "Contigs: {} total", self.total_contigs)?;
        writeln!(f, "  - Primary chromosomes: {}", self.primary_chromosomes)?;
        writeln!(f, "  - ALT contigs: {}", self.alt_contigs)?;
        writeln!(f, "  - Other: {}", self.other_contigs)?;
        writeln!(f)?;

        writeln!(f, "Coverage:")?;
        let pct = |n: usize, total: usize| {
            if total == 0 {
                0
            } else {
                (n * 100) / total
            }
        };
        let check = |n: usize, total: usize| {
            if n == total {
                "+"
            } else {
                "o"
            }
        };
        writeln!(
            f,
            "  {} Length:  {}/{} ({}%)",
            check(self.with_length, self.total_contigs),
            self.with_length,
            self.total_contigs,
            pct(self.with_length, self.total_contigs)
        )?;
        writeln!(
            f,
            "  {} MD5:     {}/{} ({}%)",
            check(self.with_md5, self.total_contigs),
            self.with_md5,
            self.total_contigs,
            pct(self.with_md5, self.total_contigs)
        )?;
        writeln!(
            f,
            "  {} Aliases: {}/{} ({}%)",
            check(self.with_aliases, self.total_contigs),
            self.with_aliases,
            self.total_contigs,
            pct(self.with_aliases, self.total_contigs)
        )?;
        writeln!(f)?;

        writeln!(f, "Conflicts: {}", self.conflicts.len())?;
        for conflict in &self.conflicts {
            writeln!(f, "  - {conflict}")?;
        }

        let total_warnings = self.warnings.len() + self.missing_md5_assembled.len();
        writeln!(f, "Warnings: {total_warnings}")?;
        for warning in &self.warnings {
            writeln!(f, "  - {warning}")?;
        }
        if !self.missing_md5_assembled.is_empty() {
            writeln!(
                f,
                "  - Missing MD5 for assembled-molecule contigs: {}",
                self.missing_md5_assembled.join(", ")
            )?;
        }

        Ok(())
    }
}

/// Builder for creating `FastaDistribution` from multiple input files
///
/// This builder merges contig metadata from multiple sources (dict, VCF, SAM/BAM, etc.)
/// keyed by (name, length). It handles:
/// - MD5 merging (first non-empty wins, error on conflict)
/// - Alias merging (union of all)
/// - Sort order preservation (first seen order)
#[derive(Debug)]
pub struct DistributionBuilder {
    id: String,
    display_name: String,
    source: ReferenceSource,
    download_url: Option<String>,
    tags: Vec<String>,

    /// Contigs keyed by (name, length) for merging
    contigs: HashMap<(String, u64), FastaContig>,

    /// Track insertion order
    insertion_order: Vec<(String, u64)>,

    /// Source files processed
    source_files: Vec<String>,
}

impl Default for DistributionBuilder {
    fn default() -> Self {
        Self::new("")
    }
}

impl DistributionBuilder {
    /// Create a new builder with the given distribution ID
    pub fn new(id: impl Into<String>) -> Self {
        Self {
            id: id.into(),
            display_name: String::new(),
            source: ReferenceSource::Custom("custom".to_string()),
            download_url: None,
            tags: Vec::new(),
            contigs: HashMap::new(),
            insertion_order: Vec::new(),
            source_files: Vec::new(),
        }
    }

    /// Set the display name
    #[must_use]
    pub fn with_display_name(mut self, name: impl Into<String>) -> Self {
        self.display_name = name.into();
        self
    }

    /// Set the reference source
    #[must_use]
    pub fn with_source(mut self, source: ReferenceSource) -> Self {
        self.source = source;
        self
    }

    /// Set the download URL
    #[must_use]
    pub fn with_download_url(mut self, url: impl Into<String>) -> Self {
        self.download_url = Some(url.into());
        self
    }

    /// Set tags
    #[must_use]
    pub fn with_tags(mut self, tags: Vec<String>) -> Self {
        self.tags = tags;
        self
    }

    /// Add an input file (auto-detect format)
    ///
    /// # Errors
    ///
    /// Returns `BuilderError::Parse` if format cannot be detected, or other
    /// errors from parsing the specific format.
    pub fn add_input(&mut self, path: &Path) -> Result<&mut Self, BuilderError> {
        let format = InputFormat::from_path(path).ok_or_else(|| {
            BuilderError::Parse(format!("Cannot detect format for: {}", path.display()))
        })?;
        self.add_input_with_format(path, format)
    }

    /// Add an input file with explicit format
    ///
    /// # Errors
    ///
    /// Returns `BuilderError::Io` if the file cannot be read, `BuilderError::Parse`
    /// if parsing fails, or `BuilderError::Conflict` if contig data conflicts.
    pub fn add_input_with_format(
        &mut self,
        path: &Path,
        format: InputFormat,
    ) -> Result<&mut Self, BuilderError> {
        let path_str = path.display().to_string();
        self.source_files.push(path_str.clone());

        let contigs = self.parse_input(path, format)?;

        for contig in contigs {
            let key = (contig.name.clone(), contig.length);
            #[allow(clippy::cast_possible_truncation)] // Contig count limited by MAX_CONTIGS (50k)
            let fasta_contig = FastaContig {
                name: contig.name,
                length: contig.length,
                md5: contig.md5.unwrap_or_default(),
                sort_order: self.insertion_order.len() as u32,
                report_contig_id: None,
                aliases: contig.aliases,
            };

            if let Some(existing) = self.contigs.get_mut(&key) {
                existing.merge(&fasta_contig)?;
            } else {
                self.insertion_order.push(key.clone());
                self.contigs.insert(key, fasta_contig);
            }
        }

        Ok(self)
    }

    /// Parse contigs from an input file
    #[allow(clippy::unused_self)] // Kept as method for API consistency
    fn parse_input(&self, path: &Path, format: InputFormat) -> Result<Vec<Contig>, BuilderError> {
        match format {
            InputFormat::Dict | InputFormat::Sam => {
                let content = std::fs::read_to_string(path)?;
                let query = crate::parsing::sam::parse_header_text(&content)
                    .map_err(|e| BuilderError::Parse(e.to_string()))?;
                Ok(query.contigs)
            }
            InputFormat::Bam | InputFormat::Cram => {
                let query = crate::parsing::sam::parse_file(path)
                    .map_err(|e| BuilderError::Parse(e.to_string()))?;
                Ok(query.contigs)
            }
            InputFormat::Fai => {
                let content = std::fs::read_to_string(path)?;
                let query = crate::parsing::fai::parse_fai_text(&content)
                    .map_err(|e| BuilderError::Parse(e.to_string()))?;
                Ok(query.contigs)
            }
            InputFormat::Fasta => {
                // Parse FASTA with MD5 computation
                let query = crate::parsing::fasta::parse_fasta_file_with_md5(path)
                    .map_err(|e| BuilderError::Parse(e.to_string()))?;
                Ok(query.contigs)
            }
            InputFormat::NcbiReport => {
                let content = std::fs::read_to_string(path)?;
                let entries = crate::parsing::ncbi_report::parse_ncbi_report_text(&content)
                    .map_err(|e| BuilderError::Parse(e.to_string()))?;
                Ok(entries.into_iter().map(|e| e.to_contig()).collect())
            }
            InputFormat::Vcf => {
                let content = std::fs::read_to_string(path)?;
                let query = crate::parsing::vcf::parse_vcf_header_text(&content)
                    .map_err(|e| BuilderError::Parse(e.to_string()))?;
                Ok(query.contigs)
            }
            InputFormat::Tsv => {
                let content = std::fs::read_to_string(path)?;
                let query = crate::parsing::tsv::parse_tsv_text(&content, '\t')
                    .map_err(|e| BuilderError::Parse(e.to_string()))?;
                Ok(query.contigs)
            }
        }
    }

    /// Build the `FastaDistribution`
    ///
    /// # Errors
    ///
    /// Returns `BuilderError::MissingField` if no contigs were found.
    pub fn build(self) -> Result<FastaDistribution, BuilderError> {
        if self.contigs.is_empty() {
            return Err(BuilderError::MissingField("No contigs found".to_string()));
        }

        // Build contigs in insertion order
        let mut contigs: Vec<FastaContig> = Vec::with_capacity(self.insertion_order.len());
        for (i, key) in self.insertion_order.iter().enumerate() {
            if let Some(mut contig) = self.contigs.get(key).cloned() {
                #[allow(clippy::cast_possible_truncation)] // Contig count limited
                {
                    contig.sort_order = i as u32;
                }
                contigs.push(contig);
            }
        }

        Ok(FastaDistribution {
            id: self.id,
            display_name: self.display_name,
            source: self.source,
            download_url: self.download_url,
            tags: self.tags,
            contigs,
        })
    }
}

/// Detect assembly version from display name
#[must_use]
pub fn detect_assembly_from_name(display_name: &str) -> Assembly {
    let lower = display_name.to_lowercase();
    if lower.contains("chm13") || lower.contains("t2t") {
        Assembly::Other("CHM13".to_string())
    } else if lower.contains("38") || lower.contains("hg38") || lower.contains("hs38") {
        Assembly::Grch38
    } else if lower.contains("37")
        || lower.contains("19")
        || lower.contains("hg19")
        || lower.contains("hs37")
        || lower.contains("b37")
    {
        Assembly::Grch37
    } else {
        Assembly::Other(display_name.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_input_format_detection() {
        assert_eq!(
            InputFormat::from_path(Path::new("test.dict")),
            Some(InputFormat::Dict)
        );
        assert_eq!(
            InputFormat::from_path(Path::new("test.fai")),
            Some(InputFormat::Fai)
        );
        assert_eq!(
            InputFormat::from_path(Path::new("test.fa")),
            Some(InputFormat::Fasta)
        );
        assert_eq!(
            InputFormat::from_path(Path::new("test.fasta")),
            Some(InputFormat::Fasta)
        );
        assert_eq!(
            InputFormat::from_path(Path::new("test.fna")),
            Some(InputFormat::Fasta)
        );
        assert_eq!(
            InputFormat::from_path(Path::new("test.fa.gz")),
            Some(InputFormat::Fasta)
        );
        assert_eq!(
            InputFormat::from_path(Path::new("test.vcf")),
            Some(InputFormat::Vcf)
        );
        assert_eq!(
            InputFormat::from_path(Path::new("test.vcf.gz")),
            Some(InputFormat::Vcf)
        );
        assert_eq!(
            InputFormat::from_path(Path::new("GRCh38_assembly_report.txt")),
            Some(InputFormat::NcbiReport)
        );
    }

    #[test]
    fn test_detect_assembly() {
        assert!(matches!(
            detect_assembly_from_name("GRCh38 (Broad)"),
            Assembly::Grch38
        ));
        assert!(matches!(
            detect_assembly_from_name("hg38 UCSC"),
            Assembly::Grch38
        ));
        assert!(matches!(
            detect_assembly_from_name("GRCh37"),
            Assembly::Grch37
        ));
        assert!(matches!(
            detect_assembly_from_name("hg19"),
            Assembly::Grch37
        ));
        assert!(matches!(
            detect_assembly_from_name("T2T-CHM13"),
            Assembly::Other(_)
        ));
    }

    #[test]
    fn test_builder_single_input() {
        let mut builder = ReferenceBuilder::new("test_ref", "Test Reference")
            .assembly(Assembly::Grch38)
            .source(ReferenceSource::Custom("test".to_string()));

        // Manually add a contig
        let contig = Contig::new("chr1", 248_956_422);
        builder.merge_contig(&contig, "test").unwrap();

        let reference = builder.build().unwrap();
        assert_eq!(reference.id.0, "test_ref");
        assert_eq!(reference.contigs.len(), 1);
        assert_eq!(reference.contigs[0].name, "chr1");
    }

    #[test]
    fn test_builder_merge_with_aliases() {
        let mut builder = ReferenceBuilder::new("test_ref", "Test Reference");

        // Add first contig
        let contig1 = Contig::new("chr1", 248_956_422);
        builder.merge_contig(&contig1, "source1").unwrap();

        // Add same contig with different name (as alias)
        let mut contig2 = Contig::new("1", 248_956_422);
        contig2.aliases = vec!["chr1".to_string()];
        builder.merge_contig(&contig2, "source2").unwrap();

        // Should have merged into one entry
        let summary = builder.summary();
        assert_eq!(summary.total_contigs, 1);

        let reference = builder.build().unwrap();
        assert_eq!(reference.contigs.len(), 1);
        // The alias "1" should be in aliases
        assert!(reference.contigs[0].aliases.contains(&"1".to_string()));
    }

    #[test]
    fn test_builder_conflict_detection() {
        let mut builder = ReferenceBuilder::new("test_ref", "Test Reference");

        // Add first contig
        let contig1 = Contig::new("chr1", 248_956_422);
        builder.merge_contig(&contig1, "source1").unwrap();

        // Add conflicting contig (different length, same name)
        let contig2 = Contig::new("chr1", 100_000);
        let result = builder.merge_contig(&contig2, "source2");
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), BuilderError::Conflict(_)));
    }

    #[test]
    fn test_distribution_builder_single_dict() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut file = NamedTempFile::with_suffix(".dict").unwrap();
        writeln!(file, "@HD\tVN:1.6").unwrap();
        writeln!(
            file,
            "@SQ\tSN:chr1\tLN:1000\tM5:6aef897c3d6ff0c78aff06ac189178dd"
        )
        .unwrap();
        writeln!(
            file,
            "@SQ\tSN:chr2\tLN:2000\tM5:f98db672eb0993dcfdabafe2a882905c"
        )
        .unwrap();

        let mut builder = DistributionBuilder::new("test_ref");
        builder.add_input(file.path()).unwrap();
        let dist = builder.build().unwrap();

        assert_eq!(dist.contigs.len(), 2);
        assert_eq!(dist.contigs[0].name, "chr1");
        assert_eq!(dist.contigs[0].md5, "6aef897c3d6ff0c78aff06ac189178dd");
        assert_eq!(dist.contigs[1].name, "chr2");
        assert_eq!(dist.contigs[1].md5, "f98db672eb0993dcfdabafe2a882905c");
    }

    #[test]
    fn test_distribution_builder_merges_inputs() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        // Dict without MD5
        let mut dict = NamedTempFile::with_suffix(".dict").unwrap();
        writeln!(dict, "@HD\tVN:1.6").unwrap();
        writeln!(dict, "@SQ\tSN:chr1\tLN:1000").unwrap();

        // VCF with MD5
        let mut vcf = NamedTempFile::with_suffix(".vcf").unwrap();
        writeln!(vcf, "##fileformat=VCFv4.2").unwrap();
        writeln!(
            vcf,
            "##contig=<ID=chr1,length=1000,md5=6aef897c3d6ff0c78aff06ac189178dd>"
        )
        .unwrap();
        writeln!(vcf, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();

        let mut builder = DistributionBuilder::new("test_ref");
        builder.add_input(dict.path()).unwrap();
        builder.add_input(vcf.path()).unwrap();
        let distribution = builder.build().unwrap();

        assert_eq!(distribution.contigs.len(), 1);
        assert_eq!(
            distribution.contigs[0].md5,
            "6aef897c3d6ff0c78aff06ac189178dd"
        );
    }

    #[test]
    fn test_distribution_builder_md5_conflict() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut dict = NamedTempFile::with_suffix(".dict").unwrap();
        writeln!(dict, "@HD\tVN:1.6").unwrap();
        writeln!(
            dict,
            "@SQ\tSN:chr1\tLN:1000\tM5:6aef897c3d6ff0c78aff06ac189178dd"
        )
        .unwrap();

        let mut vcf = NamedTempFile::with_suffix(".vcf").unwrap();
        writeln!(vcf, "##fileformat=VCFv4.2").unwrap();
        writeln!(
            vcf,
            "##contig=<ID=chr1,length=1000,md5=f98db672eb0993dcfdabafe2a882905c>"
        )
        .unwrap();
        writeln!(vcf, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();

        let mut builder = DistributionBuilder::new("test_ref");
        builder.add_input(dict.path()).unwrap();
        let result = builder.add_input(vcf.path());

        assert!(matches!(result, Err(BuilderError::Merge(_))));
    }

    #[test]
    fn test_distribution_builder_preserves_order() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut file = NamedTempFile::with_suffix(".dict").unwrap();
        writeln!(file, "@HD\tVN:1.6").unwrap();
        writeln!(
            file,
            "@SQ\tSN:chrM\tLN:16569\tM5:d2ed829b8a1628d16cbeee88e88e39eb"
        )
        .unwrap();
        writeln!(
            file,
            "@SQ\tSN:chr1\tLN:248956422\tM5:6aef897c3d6ff0c78aff06ac189178dd"
        )
        .unwrap();
        writeln!(
            file,
            "@SQ\tSN:chr2\tLN:242193529\tM5:f98db672eb0993dcfdabafe2a882905c"
        )
        .unwrap();

        let mut builder = DistributionBuilder::new("test_ref");
        builder.add_input(file.path()).unwrap();
        let dist = builder.build().unwrap();

        assert_eq!(dist.contigs[0].name, "chrM");
        assert_eq!(dist.contigs[0].sort_order, 0);
        assert_eq!(dist.contigs[1].name, "chr1");
        assert_eq!(dist.contigs[1].sort_order, 1);
        assert_eq!(dist.contigs[2].name, "chr2");
        assert_eq!(dist.contigs[2].sort_order, 2);
    }
}
