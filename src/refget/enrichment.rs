//! Async enrichment logic for querying refget servers about unknown contigs.

use std::sync::Arc;

use tokio::sync::Semaphore;

use crate::core::contig::Contig;
use crate::refget::{EnrichedContig, RefgetAlias, RefgetConfig, RefgetLookupResult};

/// Look up metadata for a set of contigs from a refget server.
///
/// Only contigs that have an MD5 digest are queried. Lookups run concurrently
/// up to `config.max_concurrent` at a time. Errors are captured per-contig
/// and never propagate — the caller always gets results back.
///
/// # Panics
///
/// Panics if the internal semaphore is closed, which should not happen
/// during normal operation.
pub async fn enrich_contigs(contigs: &[Contig], config: &RefgetConfig) -> Vec<EnrichedContig> {
    let semaphore = Arc::new(Semaphore::new(config.max_concurrent));

    // Build a reqwest client with the configured timeout
    let http_client = match reqwest::Client::builder().timeout(config.timeout).build() {
        Ok(c) => c,
        Err(e) => {
            tracing::warn!("Failed to create HTTP client for refget: {e}");
            return error_for_all(contigs, &format!("Failed to create HTTP client: {e}"));
        }
    };

    let refget_client =
        match refget_client::RefgetClient::with_client(http_client, &config.server_url) {
            Ok(c) => Arc::new(c),
            Err(e) => {
                tracing::warn!("Failed to create refget client: {e}");
                return error_for_all(contigs, &format!("Failed to create refget client: {e}"));
            }
        };

    let mut join_set = tokio::task::JoinSet::new();

    for (idx, contig) in contigs.iter().enumerate() {
        let md5 = match &contig.md5 {
            Some(md5) => md5.clone(),
            None => {
                // No MD5 — we can't query refget, skip this contig.
                // We'll handle these in the result collection below.
                continue;
            }
        };

        let sem = Arc::clone(&semaphore);
        let client = Arc::clone(&refget_client);

        join_set.spawn(async move {
            let _permit = sem.acquire().await.expect("semaphore closed unexpectedly");
            let result = lookup_single(&client, &md5).await;
            (idx, result)
        });
    }

    // Collect results, preserving original order
    let mut results: Vec<Option<RefgetLookupResult>> = vec![None; contigs.len()];

    while let Some(join_result) = join_set.join_next().await {
        match join_result {
            Ok((idx, lookup_result)) => {
                results[idx] = Some(lookup_result);
            }
            Err(e) => {
                tracing::warn!("Refget lookup task panicked: {e}");
            }
        }
    }

    contigs
        .iter()
        .enumerate()
        .map(|(idx, contig)| {
            let refget_metadata = results[idx].take().unwrap_or_else(|| {
                if contig.md5.is_none() {
                    RefgetLookupResult::Error {
                        message: "No MD5 digest available for lookup".to_string(),
                    }
                } else {
                    RefgetLookupResult::Error {
                        message: "Lookup task failed".to_string(),
                    }
                }
            });

            EnrichedContig {
                name: contig.name.clone(),
                md5: contig.md5.clone(),
                refget_metadata,
            }
        })
        .collect()
}

/// Create error results for all contigs when the client cannot be initialized.
fn error_for_all(contigs: &[Contig], message: &str) -> Vec<EnrichedContig> {
    contigs
        .iter()
        .map(|c| EnrichedContig {
            name: c.name.clone(),
            md5: c.md5.clone(),
            refget_metadata: RefgetLookupResult::Error {
                message: message.to_string(),
            },
        })
        .collect()
}

/// Query the refget server for a single MD5 digest.
async fn lookup_single(client: &refget_client::RefgetClient, md5: &str) -> RefgetLookupResult {
    match client.get_metadata(md5).await {
        Ok(Some(metadata)) => RefgetLookupResult::Found {
            aliases: metadata
                .aliases
                .into_iter()
                .map(|a| RefgetAlias {
                    naming_authority: a.naming_authority,
                    value: a.value,
                })
                .collect(),
            sha512t24u: metadata.sha512t24u,
            circular: metadata.circular,
        },
        Ok(None) => RefgetLookupResult::NotFound,
        Err(e) => {
            tracing::debug!("Refget lookup failed for {md5}: {e}");
            RefgetLookupResult::Error {
                message: e.to_string(),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::contig::Contig;

    #[tokio::test]
    async fn test_contigs_without_md5_are_skipped() {
        // Use a non-routable address so no real HTTP happens
        let config = RefgetConfig::new("http://192.0.2.1:1");
        let contigs = vec![Contig::new("chr_no_md5".to_string(), 1000)];

        let results = enrich_contigs(&contigs, &config).await;

        assert_eq!(results.len(), 1);
        assert_eq!(results[0].name, "chr_no_md5");
        assert!(results[0].md5.is_none());
        // Should be an error since no MD5 was available
        assert!(matches!(
            results[0].refget_metadata,
            RefgetLookupResult::Error { .. }
        ));
    }

    #[tokio::test]
    async fn test_invalid_server_produces_errors() {
        // Use a non-routable address with very short timeout
        let mut config = RefgetConfig::new("http://192.0.2.1:1");
        config.timeout = std::time::Duration::from_millis(100);

        let contigs = vec![
            Contig::new("chr1".to_string(), 1000).with_md5("6aef897c3d6ff0c78aff06ac189178dd")
        ];

        let results = enrich_contigs(&contigs, &config).await;

        assert_eq!(results.len(), 1);
        assert_eq!(results[0].name, "chr1");
        // Should be an error due to connection failure/timeout
        assert!(matches!(
            results[0].refget_metadata,
            RefgetLookupResult::Error { .. }
        ));
    }

    #[tokio::test]
    async fn test_mixed_contigs_with_and_without_md5() {
        let mut config = RefgetConfig::new("http://192.0.2.1:1");
        config.timeout = std::time::Duration::from_millis(100);

        let contigs = vec![
            Contig::new("chr1".to_string(), 1000).with_md5("6aef897c3d6ff0c78aff06ac189178dd"),
            Contig::new("chrUn".to_string(), 1000),
            Contig::new("chr2".to_string(), 1000).with_md5("b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2b2"),
        ];

        let results = enrich_contigs(&contigs, &config).await;

        assert_eq!(results.len(), 3);
        assert_eq!(results[0].name, "chr1");
        assert_eq!(results[1].name, "chrUn");
        assert_eq!(results[2].name, "chr2");

        // chrUn has no MD5 — should get error
        assert!(results[1].md5.is_none());
        assert!(matches!(
            results[1].refget_metadata,
            RefgetLookupResult::Error { .. }
        ));
    }

    #[tokio::test]
    async fn test_concurrency_limit_respected() {
        let mut config = RefgetConfig::new("http://192.0.2.1:1");
        config.timeout = std::time::Duration::from_millis(100);
        config.max_concurrent = 2; // Limit to 2

        // Create more contigs than the concurrency limit
        let contigs: Vec<Contig> = (0..5)
            .map(|i| Contig::new(format!("chr{i}"), 1000).with_md5(format!("{i:032x}")))
            .collect();

        let results = enrich_contigs(&contigs, &config).await;

        // All should complete (with errors due to unreachable server)
        assert_eq!(results.len(), 5);
        for result in &results {
            assert!(matches!(
                result.refget_metadata,
                RefgetLookupResult::Error { .. }
            ));
        }
    }
}
