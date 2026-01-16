use std::collections::{HashMap, HashSet};

use crate::core::contig::Contig;
use crate::core::header::QueryHeader;
use crate::core::reference::KnownReference;
use crate::core::types::MatchType;

/// Detailed diagnosis of differences between query and reference
#[derive(Debug, Clone)]
pub struct MatchDiagnosis {
    /// Type of match
    pub match_type: MatchType,

    /// Query contigs that match reference exactly
    pub exact_matches: Vec<ContigMatch>,

    /// Query contigs that match by MD5 but have different names
    pub renamed_matches: Vec<RenamedContig>,

    /// Query contigs that match by name+length but MD5 differs or is missing
    pub name_length_only_matches: Vec<ContigMatch>,

    /// Query contigs with no match in reference
    pub query_only: Vec<Contig>,

    /// Contigs that are in a different order
    pub reordered: bool,

    /// Potential source conflicts (e.g., mito from different build)
    pub conflicts: Vec<ContigConflict>,

    /// Suggested fixes
    pub suggestions: Vec<Suggestion>,
}

#[derive(Debug, Clone)]
pub struct ContigMatch;

#[derive(Debug, Clone)]
pub struct RenamedContig {
    pub query_name: String,
    pub reference_name: String,
}

#[derive(Debug, Clone)]
pub struct ContigConflict {
    pub query_contig: Contig,
    pub expected: Option<Contig>,
    pub conflict_type: ConflictType,
    pub description: String,
}

#[derive(Debug, Clone)]
pub enum ConflictType {
    /// Same name, different sequence (different MD5/length)
    SequenceMismatch,
    /// Mitochondrial from different build
    MitochondrialMismatch,
    /// Unknown contig not in reference
    UnknownContig,
}

#[derive(Debug, Clone)]
pub enum Suggestion {
    /// Rename contigs using fgbio or similar
    RenameContigs { command_hint: String },
    /// Reorder contigs to match reference
    ReorderContigs { command_hint: String },
    /// Replace a specific contig (e.g., mito)
    ReplaceContig {
        contig_name: String,
        reason: String,
        source: String,
    },
    /// Reference is a close match, safe to use
    UseAsIs { warnings: Vec<String> },
    /// Need to realign with different reference
    Realign {
        reason: String,
        suggested_reference: String,
    },
}

impl MatchDiagnosis {
    #[must_use]
    #[allow(clippy::too_many_lines)] // TODO: Refactor into smaller functions
    pub fn analyze(query: &QueryHeader, reference: &KnownReference) -> Self {
        let mut exact_matches = Vec::new();
        let mut renamed_matches = Vec::new();
        let mut name_length_only_matches = Vec::new();
        let mut query_only = Vec::new();
        let mut conflicts = Vec::new();

        // Build lookup maps for reference
        let ref_by_md5: HashMap<&str, &Contig> = reference
            .contigs
            .iter()
            .filter_map(|c| c.md5.as_ref().map(|m| (m.as_str(), c)))
            .collect();

        // Use exact names for matching (no normalization)
        let ref_by_name_length: HashMap<(String, u64), &Contig> = reference
            .contigs
            .iter()
            .map(|c| ((c.name.clone(), c.length), c))
            .collect();

        let mut matched_ref_md5s: HashSet<&str> = HashSet::new();
        let mut matched_ref_name_lengths: HashSet<(String, u64)> = HashSet::new();

        // Analyze each query contig using exact names
        for q_contig in &query.contigs {
            let q_key = (q_contig.name.clone(), q_contig.length);

            // Try MD5 match first
            if let Some(q_md5) = &q_contig.md5 {
                if let Some(r_contig) = ref_by_md5.get(q_md5.as_str()) {
                    matched_ref_md5s.insert(q_md5.as_str());

                    if q_contig.name == r_contig.name {
                        // Exact match
                        exact_matches.push(ContigMatch);
                    } else {
                        // Same sequence, different name
                        renamed_matches.push(RenamedContig {
                            query_name: q_contig.name.clone(),
                            reference_name: r_contig.name.clone(),
                        });
                    }
                    continue;
                }
            }

            // Try name+length match (direct name or via alias)
            let matched_ref = ref_by_name_length.get(&q_key).copied().or_else(|| {
                // Try matching via query aliases (reverse alias matching)
                q_contig.aliases.iter().find_map(|alias| {
                    let alias_key = (alias.clone(), q_contig.length);
                    ref_by_name_length.get(&alias_key).copied()
                })
            });

            if let Some(r_contig) = matched_ref {
                let matched_key = (r_contig.name.clone(), r_contig.length);
                matched_ref_name_lengths.insert(matched_key);

                // Check if this might be a conflict (same position but different sequence)
                if let (Some(q_md5), Some(r_md5)) = (&q_contig.md5, &r_contig.md5) {
                    if q_md5 != r_md5 {
                        // Different sequence!
                        let conflict_type = if q_contig.is_mitochondrial() {
                            ConflictType::MitochondrialMismatch
                        } else {
                            ConflictType::SequenceMismatch
                        };

                        conflicts.push(ContigConflict {
                            query_contig: q_contig.clone(),
                            expected: Some(r_contig.clone()),
                            conflict_type,
                            description: format!(
                                "Contig {} has same name/length but different MD5 (query: {}, ref: {})",
                                q_contig.name, q_md5, r_md5
                            ),
                        });
                        continue;
                    }
                }

                // If names differ but match via alias, count as renamed
                if q_contig.name == r_contig.name {
                    name_length_only_matches.push(ContigMatch);
                } else {
                    renamed_matches.push(RenamedContig {
                        query_name: q_contig.name.clone(),
                        reference_name: r_contig.name.clone(),
                    });
                }
                continue;
            }

            // No match found
            let conflict_type = if q_contig.is_mitochondrial() {
                // Special handling for mitochondrial
                ConflictType::MitochondrialMismatch
            } else {
                ConflictType::UnknownContig
            };

            if q_contig.is_primary_chromosome() || q_contig.is_mitochondrial() {
                conflicts.push(ContigConflict {
                    query_contig: q_contig.clone(),
                    expected: None,
                    conflict_type,
                    description: format!(
                        "No match found for {} (length: {})",
                        q_contig.name, q_contig.length
                    ),
                });
            } else {
                query_only.push(q_contig.clone());
            }
        }

        // Determine match type and reordering
        let reordered = !super::scoring::MatchScore::calculate(query, reference).order_preserved;

        let match_type = determine_match_type(
            &exact_matches,
            &renamed_matches,
            &name_length_only_matches,
            &query_only,
            &conflicts,
            reordered,
        );

        // Generate suggestions
        let suggestions = generate_suggestions(
            &match_type,
            &renamed_matches,
            &conflicts,
            reordered,
            reference,
        );

        Self {
            match_type,
            exact_matches,
            renamed_matches,
            name_length_only_matches,
            query_only,
            reordered,
            conflicts,
            suggestions,
        }
    }
}

fn determine_match_type(
    exact_matches: &[ContigMatch],
    renamed_matches: &[RenamedContig],
    name_length_only: &[ContigMatch],
    query_only: &[Contig],
    conflicts: &[ContigConflict],
    reordered: bool,
) -> MatchType {
    let total_matched = exact_matches.len() + renamed_matches.len() + name_length_only.len();

    // No matches at all
    if total_matched == 0 {
        return MatchType::NoMatch;
    }

    // All exact matches
    if !exact_matches.is_empty()
        && renamed_matches.is_empty()
        && conflicts.is_empty()
        && query_only.is_empty()
    {
        if reordered {
            return MatchType::Reordered;
        }
        return MatchType::Exact;
    }

    // All matches but need renaming
    if !renamed_matches.is_empty() && conflicts.is_empty() && query_only.is_empty() {
        if reordered {
            return MatchType::ReorderedAndRenamed;
        }
        return MatchType::Renamed;
    }

    // Has conflicts or unmatched - partial or mixed
    if !conflicts.is_empty() {
        // Check if conflicts suggest mixed sources
        let has_mito_conflict = conflicts
            .iter()
            .any(|c| matches!(c.conflict_type, ConflictType::MitochondrialMismatch));

        if has_mito_conflict && total_matched > conflicts.len() * 2 {
            return MatchType::Mixed;
        }
    }

    MatchType::Partial
}

fn generate_suggestions(
    match_type: &MatchType,
    renamed_matches: &[RenamedContig],
    conflicts: &[ContigConflict],
    reordered: bool,
    reference: &KnownReference,
) -> Vec<Suggestion> {
    let mut suggestions = Vec::new();

    // Renaming suggestion
    if !renamed_matches.is_empty() {
        suggestions.push(Suggestion::RenameContigs {
            command_hint: "fgbio UpdateSequenceDictionary --in input.bam --out output.bam \\\n  \
                 --sequence-dictionary reference.dict"
                .to_string(),
        });
    }

    // Reordering suggestion
    if reordered {
        suggestions.push(Suggestion::ReorderContigs {
            command_hint: "picard ReorderSam \\\n  \
                 I=input.bam \\\n  \
                 O=output.bam \\\n  \
                 REFERENCE=reference.fa"
                .to_string(),
        });
    }

    // Conflict-specific suggestions
    for conflict in conflicts {
        match &conflict.conflict_type {
            ConflictType::MitochondrialMismatch => {
                suggestions.push(Suggestion::ReplaceContig {
                    contig_name: conflict.query_contig.name.clone(),
                    reason: format!(
                        "Mitochondrial sequence differs: {} ({}bp) vs expected ({}bp)",
                        conflict.query_contig.name,
                        conflict.query_contig.length,
                        conflict.expected.as_ref().map_or(0, |c| c.length)
                    ),
                    source: reference.download_url.clone().unwrap_or_default(),
                });
            }
            ConflictType::SequenceMismatch => {
                suggestions.push(Suggestion::Realign {
                    reason: format!(
                        "Contig {} has different sequence than expected",
                        conflict.query_contig.name
                    ),
                    suggested_reference: reference.id.to_string(),
                });
            }
            ConflictType::UnknownContig => {
                // Usually fine - just extra contigs
            }
        }
    }

    // Safe to use suggestion
    if matches!(match_type, MatchType::Exact) && conflicts.is_empty() {
        suggestions.push(Suggestion::UseAsIs {
            warnings: Vec::new(),
        });
    }

    suggestions
}
