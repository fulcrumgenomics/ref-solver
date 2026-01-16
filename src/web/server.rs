use axum::http::header;
use axum::{
    extract::{DefaultBodyLimit, Multipart, Query, State},
    http::{HeaderName, HeaderValue, StatusCode},
    response::{Html, IntoResponse, Json, Response},
    routing::{get, post},
    Router,
};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::net::SocketAddr;
use std::sync::Arc;
use std::time::Duration;
use tokio::net::TcpListener;
use tower::limit::ConcurrencyLimitLayer;
use tower::ServiceBuilder;
use tower_governor::{governor::GovernorConfigBuilder, GovernorLayer};
use tower_http::set_header::SetResponseHeaderLayer;
use tower_http::timeout::TimeoutLayer;

use crate::catalog::store::ReferenceCatalog;
use crate::cli::ServeArgs;
use crate::matching::engine::{MatchingConfig, MatchingEngine, ScoringWeights};
use crate::matching::Suggestion;
use crate::utils::validation::{validate_upload, ValidationError};
use crate::web::format_detection::{
    detect_format, parse_binary_file, parse_with_format, FileFormat,
};

/// Security configuration constants to prevent `DoS` attacks
pub const MAX_MULTIPART_FIELDS: usize = 10;
pub const MAX_FILE_FIELD_SIZE: usize = 16 * 1024 * 1024; // 16MB
pub const MAX_TEXT_FIELD_SIZE: usize = 1024 * 1024; // 1MB

/// Helper function to convert usize count to f64 with explicit precision loss allowance
#[inline]
fn count_to_f64(count: usize) -> f64 {
    #[allow(clippy::cast_precision_loss)]
    {
        count as f64
    }
}

/// Shared application state
pub struct AppState {
    pub catalog: ReferenceCatalog,
}

/// Input data extracted from multipart form
#[derive(Debug)]
struct InputData {
    /// Text content (if provided via textarea or text file)
    text_content: Option<String>,
    /// Binary file content (if provided)
    binary_content: Option<Vec<u8>>,
    /// Original filename
    filename: Option<String>,
    /// Detected or specified format
    format: Option<FileFormat>,
}

/// Enhanced error response
#[derive(Serialize)]
pub struct ErrorResponse {
    pub error: String,
    pub error_type: String,
    pub details: Option<String>,
}

#[derive(Serialize)]
struct ConfigurationInfo {
    score_threshold: f64,
    result_limit: usize,
    scoring_weights: ScoringWeights,
}

/// Query parameters for detailed mode
#[derive(Deserialize)]
struct DetailedQueryParams {
    /// Mode: "detailed" for detailed contig breakdown, omit for summary
    mode: Option<String>,
    /// Match index to get details for (0-based)
    match_id: Option<usize>,
    /// Page number for query contigs (0-based)
    query_page: Option<usize>,
    /// Page size for query contigs (default: 100, max: 500)
    query_page_size: Option<usize>,
    /// Page number for reference contigs (0-based)
    ref_page: Option<usize>,
    /// Page size for reference contigs (default: 100, max: 500)
    ref_page_size: Option<usize>,
}

/// Create a safe error response that prevents information disclosure
/// while logging detailed errors server-side for debugging
pub fn create_safe_error_response(
    error_type: &str,
    user_message: &str,
    internal_error: Option<&str>,
) -> ErrorResponse {
    // Log detailed error server-side for debugging (not exposed to client)
    if let Some(internal_msg) = internal_error {
        tracing::error!("Internal error ({}): {}", error_type, internal_msg);
    }

    ErrorResponse {
        error: user_message.to_string(),
        error_type: error_type.to_string(),
        details: None, // Never expose internal details to prevent information disclosure
    }
}

/// Run the web server
///
/// # Errors
///
/// Returns an error if the tokio runtime cannot be created or the server fails to start.
pub fn run(args: ServeArgs) -> anyhow::Result<()> {
    // Build tokio runtime
    let rt = tokio::runtime::Runtime::new()?;
    rt.block_on(async move { run_server(args).await })
}

/// Create the application router with all routes and middleware configured.
///
/// # Errors
///
/// Returns an error if the catalog cannot be loaded.
#[allow(clippy::missing_panics_doc)] // Panics only on invalid governor config (constants are valid)
pub fn create_router() -> anyhow::Result<Router> {
    // Load catalog
    let catalog = ReferenceCatalog::load_embedded()?;
    let state = Arc::new(AppState { catalog });

    // Configure IP-based rate limiting
    let governor_conf = GovernorConfigBuilder::default()
        .per_second(10) // 10 requests per second per IP
        .burst_size(50) // Allow bursts of 50 requests
        .finish()
        .unwrap();

    // Build router with comprehensive security layers
    let app = Router::new()
        .route("/", get(index_handler))
        .route("/api/identify", post(identify_handler))
        .route("/api/catalog", get(catalog_handler))
        // Static file routes
        .route("/static/css/styles.css", get(styles_css_handler))
        .route("/static/js/main.js", get(main_js_handler))
        .route("/static/js/utils/helpers.js", get(helpers_js_handler))
        .route(
            "/static/js/managers/ConfigurationManager.js",
            get(config_manager_js_handler),
        )
        .route(
            "/static/js/managers/TabManager.js",
            get(tab_manager_js_handler),
        )
        .route(
            "/static/js/managers/ResultsManager.js",
            get(results_manager_js_handler),
        )
        .route(
            "/static/js/managers/SplitViewManager.js",
            get(split_view_manager_js_handler),
        )
        .with_state(state)
        .layer(
            ServiceBuilder::new()
                // Security headers for browser protection
                .layer(SetResponseHeaderLayer::if_not_present(
                    HeaderName::from_static("x-content-type-options"),
                    HeaderValue::from_static("nosniff"),
                ))
                .layer(SetResponseHeaderLayer::if_not_present(
                    HeaderName::from_static("x-frame-options"),
                    HeaderValue::from_static("DENY"),
                ))
                .layer(SetResponseHeaderLayer::if_not_present(
                    HeaderName::from_static("x-xss-protection"),
                    HeaderValue::from_static("1; mode=block"),
                ))
                .layer(SetResponseHeaderLayer::if_not_present(
                    HeaderName::from_static("strict-transport-security"),
                    HeaderValue::from_static("max-age=31536000; includeSubDomains"),
                ))
                .layer(SetResponseHeaderLayer::if_not_present(
                    HeaderName::from_static("referrer-policy"),
                    HeaderValue::from_static("strict-origin-when-cross-origin"),
                ))
                // IP-based rate limiting to prevent abuse
                .layer(GovernorLayer {
                    config: Arc::new(governor_conf),
                })
                // Request timeout to prevent slow client attacks
                .layer(TimeoutLayer::with_status_code(
                    StatusCode::REQUEST_TIMEOUT,
                    Duration::from_secs(30),
                ))
                // Limit concurrent requests to prevent DOS
                .layer(ConcurrencyLimitLayer::new(100))
                // Limit request body size (accommodate largest file + multipart overhead)
                .layer(DefaultBodyLimit::max(20 * 1024 * 1024)), // 20MB limit
        );

    Ok(app)
}

async fn run_server(args: ServeArgs) -> anyhow::Result<()> {
    let app = create_router()?;

    let addr = format!("{}:{}", args.address, args.port);
    println!("Starting ref-solver web server at http://{addr}");

    if args.open {
        let _ = open::that(format!("http://{addr}"));
    }

    let listener = TcpListener::bind(&addr).await?;
    axum::serve(
        listener,
        app.into_make_service_with_connect_info::<SocketAddr>(),
    )
    .await?;

    Ok(())
}

/// Main page handler
async fn index_handler() -> Html<&'static str> {
    Html(include_str!("templates/index.html"))
}

/// Static CSS handler
async fn styles_css_handler() -> impl IntoResponse {
    (
        [(header::CONTENT_TYPE, "text/css; charset=utf-8")],
        include_str!("static/css/styles.css"),
    )
}

/// Static JS handlers for ES6 modules
async fn main_js_handler() -> impl IntoResponse {
    (
        [(
            header::CONTENT_TYPE,
            "application/javascript; charset=utf-8",
        )],
        include_str!("static/js/main.js"),
    )
}

async fn helpers_js_handler() -> impl IntoResponse {
    (
        [(
            header::CONTENT_TYPE,
            "application/javascript; charset=utf-8",
        )],
        include_str!("static/js/utils/helpers.js"),
    )
}

async fn config_manager_js_handler() -> impl IntoResponse {
    (
        [(
            header::CONTENT_TYPE,
            "application/javascript; charset=utf-8",
        )],
        include_str!("static/js/managers/ConfigurationManager.js"),
    )
}

async fn tab_manager_js_handler() -> impl IntoResponse {
    (
        [(
            header::CONTENT_TYPE,
            "application/javascript; charset=utf-8",
        )],
        include_str!("static/js/managers/TabManager.js"),
    )
}

async fn results_manager_js_handler() -> impl IntoResponse {
    (
        [(
            header::CONTENT_TYPE,
            "application/javascript; charset=utf-8",
        )],
        include_str!("static/js/managers/ResultsManager.js"),
    )
}

async fn split_view_manager_js_handler() -> impl IntoResponse {
    (
        [(
            header::CONTENT_TYPE,
            "application/javascript; charset=utf-8",
        )],
        include_str!("static/js/managers/SplitViewManager.js"),
    )
}

/// API endpoint for identifying references
#[allow(clippy::too_many_lines)] // TODO: Refactor into smaller functions
async fn identify_handler(
    State(state): State<Arc<AppState>>,
    Query(params): Query<DetailedQueryParams>,
    mut multipart: Multipart,
) -> impl IntoResponse {
    let start_time = std::time::Instant::now();

    // Extract input data and configuration from multipart form
    let (input_data, config) = match extract_request_data(&mut multipart).await {
        Ok(data) => data,
        Err(error_response) => return error_response,
    };

    // Parse input using intelligent format detection
    let query = match parse_input_data(&input_data) {
        Ok(query) => query,
        Err(error_response) => return *error_response,
    };

    // Create matching engine with configuration
    let matching_config = MatchingConfig {
        min_score: config.score_threshold,
        scoring_weights: config.scoring_weights.clone(),
    };

    let engine = MatchingEngine::new(&state.catalog, matching_config);
    let matches = engine.find_matches(&query, config.result_limit);

    // Check if detailed mode is requested
    if params.mode.as_deref() == Some("detailed") {
        return handle_detailed_response(&params, &matches, &query, start_time, &config).await;
    }

    // Build enhanced response
    let results: Vec<serde_json::Value> = matches
        .iter()
        .map(|m| {
            serde_json::json!({
                "reference": {
                    "id": m.reference.id.0,
                    "display_name": m.reference.display_name,
                    "assembly": format!("{}", m.reference.assembly),
                    "source": format!("{}", m.reference.source),
                    "download_url": m.reference.download_url,
                },
                "score": {
                    "composite": m.score.composite,
                    "confidence": format!("{:?}", m.score.confidence),
                    "detailed_scores": {
                        "md5_jaccard": m.score.md5_jaccard,
                        "name_length_jaccard": m.score.name_length_jaccard,
                        "md5_query_coverage": m.score.md5_query_coverage,
                        "order_score": m.score.order_score,
                    },
                },
                "match_type": format!("{:?}", m.diagnosis.match_type),
                "reordered": m.diagnosis.reordered,
                "exact_matches": m.diagnosis.exact_matches.len(),
                "renamed_matches": m.diagnosis.renamed_matches.len(),
                "conflicts": m.diagnosis.conflicts.len(),
                "query_only": m.diagnosis.query_only.len(),
                "diagnosis": {
                    "exact_matches": m.diagnosis.exact_matches.iter().map(|_| {
                        serde_json::json!({"type": "exact"})
                    }).collect::<Vec<_>>(),
                    "renamed_matches": m.diagnosis.renamed_matches.iter().map(|r| {
                        serde_json::json!({
                            "query_name": r.query_name,
                            "reference_name": r.reference_name
                        })
                    }).collect::<Vec<_>>(),
                    "conflicts": m.diagnosis.conflicts.iter().map(|c| {
                        serde_json::json!({
                            "query_contig": {
                                "name": c.query_contig.name,
                                "length": c.query_contig.length,
                                "md5": c.query_contig.md5
                            },
                            "conflict_type": format!("{:?}", c.conflict_type),
                            "description": c.description
                        })
                    }).collect::<Vec<_>>(),
                },
                "suggestions": m.diagnosis.suggestions.iter().map(|s| {
                    match s {
                        Suggestion::RenameContigs { command_hint, .. } => {
                            serde_json::json!({"type": "rename", "command": command_hint})
                        }
                        Suggestion::ReorderContigs { command_hint } => {
                            serde_json::json!({"type": "reorder", "command": command_hint})
                        }
                        Suggestion::ReplaceContig { contig_name, reason, source } => {
                            serde_json::json!({"type": "replace", "contig": contig_name, "reason": reason, "source": source})
                        }
                        Suggestion::UseAsIs { warnings } => {
                            serde_json::json!({"type": "use_as_is", "warnings": warnings})
                        }
                        Suggestion::Realign { reason, suggested_reference } => {
                            serde_json::json!({"type": "realign", "reason": reason, "reference": suggested_reference})
                        }
                    }
                }).collect::<Vec<_>>(),
            })
        })
        .collect();

    #[allow(clippy::cast_possible_truncation)] // Processing time won't exceed u64
    let processing_time = start_time.elapsed().as_millis() as u64;

    Json(serde_json::json!({
        "query": {
            "contig_count": query.contigs.len(),
            "has_md5": query.has_md5s(),
            "md5_coverage": query.md5_coverage(),
            "naming_convention": format!("{:?}", query.naming_convention),
        },
        "matches": results,
        "processing_info": {
            "detected_format": input_data.format.as_ref().map_or("unknown", super::format_detection::FileFormat::display_name),
            "processing_time_ms": processing_time,
            "configuration": {
                "score_threshold": config.score_threshold,
                "result_limit": config.result_limit,
                "scoring_weights": config.scoring_weights,
            }
        }
    }))
    .into_response()
}

/// Handle detailed response mode for contig breakdown
#[allow(
    clippy::cast_possible_truncation,
    clippy::unused_async,
    clippy::too_many_lines
)] // JSON indices; TODO: refactor
async fn handle_detailed_response(
    params: &DetailedQueryParams,
    matches: &[crate::matching::engine::MatchResult],
    query: &crate::core::header::QueryHeader,
    start_time: std::time::Instant,
    config: &ConfigurationInfo,
) -> Response {
    use crate::core::contig::Contig;

    // Get the specific match or default to first match
    let match_index = params.match_id.unwrap_or(0);
    let Some(selected_match) = matches.get(match_index) else {
        return (
            StatusCode::BAD_REQUEST,
            Json(create_safe_error_response(
                "invalid_match_id",
                "Invalid match ID specified",
                Some("Match index out of bounds"),
            )),
        )
            .into_response();
    };

    // Set up pagination parameters
    let query_page = params.query_page.unwrap_or(0);
    let query_page_size = params.query_page_size.unwrap_or(100).min(500);
    let ref_page = params.ref_page.unwrap_or(0);
    let ref_page_size = params.ref_page_size.unwrap_or(100).min(500);

    // Extract query contigs with pagination
    let total_query_contigs = query.contigs.len();
    let query_start = query_page * query_page_size;
    let query_end = (query_start + query_page_size).min(total_query_contigs);
    let query_contigs_page: Vec<&Contig> = if query_start < total_query_contigs {
        query.contigs[query_start..query_end].iter().collect()
    } else {
        Vec::new()
    };

    // Extract reference contigs with pagination
    let total_ref_contigs = selected_match.reference.contigs.len();
    let ref_start = ref_page * ref_page_size;
    let ref_end = (ref_start + ref_page_size).min(total_ref_contigs);
    let ref_contigs_page: Vec<&Contig> = if ref_start < total_ref_contigs {
        selected_match.reference.contigs[ref_start..ref_end]
            .iter()
            .collect()
    } else {
        Vec::new()
    };

    // Build detailed mapping information
    let mut exact_match_mappings = Vec::new();
    let mut renamed_match_mappings = Vec::new();
    let mut conflict_mappings = Vec::new();
    let mut query_only_indices = Vec::new();
    let mut reference_only_indices = Vec::new();

    // Create lookup maps for efficient indexing
    let query_name_to_index: std::collections::HashMap<&str, usize> = query
        .contigs
        .iter()
        .enumerate()
        .map(|(i, c)| (c.name.as_str(), i))
        .collect();

    let ref_name_to_index: std::collections::HashMap<&str, usize> = selected_match
        .reference
        .contigs
        .iter()
        .enumerate()
        .map(|(i, c)| (c.name.as_str(), i))
        .collect();

    // Process exact matches (need to map back to contigs since ContigMatch is empty)
    for (i, _) in selected_match.diagnosis.exact_matches.iter().enumerate() {
        // Since ContigMatch doesn't contain contig data, we need to reconstruct
        // the mapping by analyzing the query and reference contigs
        // This is a limitation of the current data structure
        exact_match_mappings.push(serde_json::json!({
            "type": "exact",
            "query_index": i, // This is approximate - we'd need better data structure
            "reference_index": i // This is approximate - we'd need better data structure
        }));
    }

    // Process renamed matches
    for rename in &selected_match.diagnosis.renamed_matches {
        if let (Some(&query_idx), Some(&ref_idx)) = (
            query_name_to_index.get(rename.query_name.as_str()),
            ref_name_to_index.get(rename.reference_name.as_str()),
        ) {
            renamed_match_mappings.push(serde_json::json!({
                "type": "renamed",
                "query_index": query_idx,
                "reference_index": ref_idx,
                "query_name": rename.query_name,
                "reference_name": rename.reference_name
            }));
        }
    }

    // Process conflicts
    for conflict in &selected_match.diagnosis.conflicts {
        if let Some(&query_idx) = query_name_to_index.get(conflict.query_contig.name.as_str()) {
            let ref_idx = conflict
                .expected
                .as_ref()
                .and_then(|expected| ref_name_to_index.get(expected.name.as_str()));

            conflict_mappings.push(serde_json::json!({
                "type": "conflict",
                "query_index": query_idx,
                "reference_index": ref_idx,
                "conflict_type": format!("{:?}", conflict.conflict_type),
                "description": conflict.description
            }));
        }
    }

    // Process query-only contigs
    for contig in &selected_match.diagnosis.query_only {
        if let Some(&index) = query_name_to_index.get(contig.name.as_str()) {
            query_only_indices.push(index);
        }
    }

    // Identify reference-only contigs (those not matched by any query contig)
    let mut matched_ref_indices = std::collections::HashSet::new();
    #[allow(clippy::cast_possible_truncation)] // Contig indices bounded by MAX_CONTIGS
    for mapping in &exact_match_mappings {
        if let Some(ref_idx) = mapping
            .get("reference_index")
            .and_then(serde_json::Value::as_u64)
        {
            matched_ref_indices.insert(ref_idx as usize);
        }
    }
    #[allow(clippy::cast_possible_truncation)]
    for mapping in &renamed_match_mappings {
        if let Some(ref_idx) = mapping
            .get("reference_index")
            .and_then(serde_json::Value::as_u64)
        {
            matched_ref_indices.insert(ref_idx as usize);
        }
    }
    #[allow(clippy::cast_possible_truncation)]
    for mapping in &conflict_mappings {
        if let Some(ref_idx) = mapping
            .get("reference_index")
            .and_then(serde_json::Value::as_u64)
        {
            matched_ref_indices.insert(ref_idx as usize);
        }
    }

    for (i, _) in selected_match.reference.contigs.iter().enumerate() {
        if !matched_ref_indices.contains(&i) {
            reference_only_indices.push(i);
        }
    }

    // Build response
    #[allow(clippy::cast_possible_truncation)] // Processing time won't exceed u64
    let processing_time = start_time.elapsed().as_millis() as u64;

    Json(serde_json::json!({
        "mode": "detailed",
        "match_id": match_index,
        "query": {
            "contigs": query_contigs_page.iter().enumerate().map(|(page_idx, contig)| {
                let global_idx = query_start + page_idx;
                // Determine match status for this contig
                let match_status = if query_only_indices.contains(&global_idx) {
                    "missing"
                } else if conflict_mappings.iter().any(|c| c.get("query_index").and_then(serde_json::Value::as_u64).map(|i| i as usize) == Some(global_idx)) {
                    "conflict"
                } else if renamed_match_mappings.iter().any(|r| r.get("query_index").and_then(serde_json::Value::as_u64).map(|i| i as usize) == Some(global_idx)) {
                    "renamed"
                } else if exact_match_mappings.iter().any(|e| e.get("query_index").and_then(serde_json::Value::as_u64).map(|i| i as usize) == Some(global_idx)) {
                    "exact"
                } else {
                    "unknown"
                };

                serde_json::json!({
                    "index": global_idx,
                    "name": contig.name,
                    "length": contig.length,
                    "md5": contig.md5,
                    "sequence_role": format!("{:?}", contig.sequence_role),
                    "aliases": contig.aliases,
                    "match_status": match_status
                })
            }).collect::<Vec<_>>(),
            "pagination": {
                "page": query_page,
                "page_size": query_page_size,
                "total_count": total_query_contigs,
                "total_pages": total_query_contigs.div_ceil(query_page_size)
            }
        },
        "reference": {
            "id": selected_match.reference.id.0,
            "display_name": selected_match.reference.display_name,
            "assembly": format!("{}", selected_match.reference.assembly),
            "contigs": ref_contigs_page.iter().enumerate().map(|(page_idx, contig)| {
                let global_idx = ref_start + page_idx;
                // Determine match status for this reference contig
                let match_status = if reference_only_indices.contains(&global_idx) {
                    "missing"
                } else if conflict_mappings.iter().any(|c| c.get("reference_index").and_then(serde_json::Value::as_u64).map(|i| i as usize) == Some(global_idx)) {
                    "conflict"
                } else if renamed_match_mappings.iter().any(|r| r.get("reference_index").and_then(serde_json::Value::as_u64).map(|i| i as usize) == Some(global_idx)) {
                    "renamed"
                } else if exact_match_mappings.iter().any(|e| e.get("reference_index").and_then(serde_json::Value::as_u64).map(|i| i as usize) == Some(global_idx)) {
                    "exact"
                } else {
                    "unknown"
                };

                serde_json::json!({
                    "index": global_idx,
                    "name": contig.name,
                    "length": contig.length,
                    "md5": contig.md5,
                    "sequence_role": format!("{:?}", contig.sequence_role),
                    "aliases": contig.aliases,
                    "match_status": match_status
                })
            }).collect::<Vec<_>>(),
            "pagination": {
                "page": ref_page,
                "page_size": ref_page_size,
                "total_count": total_ref_contigs,
                "total_pages": total_ref_contigs.div_ceil(ref_page_size)
            }
        },
        "mappings": {
            "exact_matches": exact_match_mappings,
            "renamed_matches": renamed_match_mappings,
            "conflicts": conflict_mappings,
            "query_only": query_only_indices,
            "reference_only": reference_only_indices
        },
        "match_summary": {
            "match_type": format!("{:?}", selected_match.diagnosis.match_type),
            "reordered": selected_match.diagnosis.reordered,
            "score": {
                "composite": selected_match.score.composite,
                "confidence": format!("{:?}", selected_match.score.confidence)
            }
        },
        "processing_info": {
            "processing_time_ms": processing_time,
            "configuration": {
                "score_threshold": config.score_threshold,
                "result_limit": config.result_limit,
                "scoring_weights": config.scoring_weights,
            }
        }
    }))
    .into_response()
}

/// Extract input data and configuration from multipart form
#[allow(clippy::too_many_lines)] // TODO: Refactor into smaller functions
async fn extract_request_data(
    multipart: &mut Multipart,
) -> Result<(InputData, ConfigurationInfo), Response> {
    let mut input_data = InputData {
        text_content: None,
        binary_content: None,
        filename: None,
        format: None,
    };

    let mut config = ConfigurationInfo {
        score_threshold: 0.1, // Default 10%
        result_limit: 10,
        scoring_weights: ScoringWeights::default(),
    };

    let mut fields_received = 0usize;
    let mut had_parse_error = false;

    // Process multipart fields
    loop {
        // Check field count limit before processing
        if fields_received >= MAX_MULTIPART_FIELDS {
            return Err((
                StatusCode::BAD_REQUEST,
                Json(ErrorResponse {
                    error: "Too many form fields".to_string(),
                    error_type: "field_limit_exceeded".to_string(),
                    details: None, // No internal details for security
                }),
            )
                .into_response());
        }

        match multipart.next_field().await {
            Ok(Some(field)) => {
                fields_received += 1;
                let name = field.name().unwrap_or_default().to_string();

                match name.as_str() {
                    "file" => {
                        let filename = field.file_name().map(std::string::ToString::to_string);

                        match field.bytes().await {
                            Ok(bytes) => {
                                // Validate field size before processing
                                if bytes.len() > MAX_FILE_FIELD_SIZE {
                                    return Err((
                                        StatusCode::PAYLOAD_TOO_LARGE,
                                        Json(ErrorResponse {
                                            error: "File size exceeds limit".to_string(),
                                            error_type: "file_too_large".to_string(),
                                            details: None,
                                        }),
                                    )
                                        .into_response());
                                }

                                // Detect format from filename for validation
                                let detected_format = if let Some(ref name) = filename {
                                    detect_binary_format(name).unwrap_or(FileFormat::Auto)
                                } else {
                                    FileFormat::Auto
                                };

                                // Use comprehensive validation function for security
                                match validate_upload(filename.as_deref(), &bytes, detected_format)
                                {
                                    Ok(validated_filename) => {
                                        input_data.filename = validated_filename;

                                        // Detect if content is binary or text
                                        if is_binary_content(&bytes) {
                                            input_data.binary_content = Some(bytes.to_vec());
                                            input_data.format = Some(detected_format);
                                        } else {
                                            input_data.text_content =
                                                Some(String::from_utf8_lossy(&bytes).to_string());
                                        }
                                    }
                                    Err(ValidationError::FilenameTooLong) => {
                                        return Err((
                                            StatusCode::BAD_REQUEST,
                                            Json(create_safe_error_response(
                                                "filename_too_long",
                                                "Filename exceeds maximum length limit",
                                                Some("Filename validation failed due to length constraints")
                                            )),
                                        ).into_response());
                                    }
                                    Err(ValidationError::InvalidFilename) => {
                                        return Err((
                                            StatusCode::BAD_REQUEST,
                                            Json(create_safe_error_response(
                                                "invalid_filename",
                                                "Filename contains invalid or dangerous characters",
                                                Some("Filename validation failed due to invalid characters")
                                            )),
                                        ).into_response());
                                    }
                                    Err(ValidationError::FormatValidationFailed) => {
                                        return Err((
                                            StatusCode::BAD_REQUEST,
                                            Json(create_safe_error_response(
                                                "format_mismatch",
                                                "File content does not match the expected format based on filename",
                                                Some("Format validation failed")
                                            )),
                                        ).into_response());
                                    }
                                    Err(ValidationError::InvalidFileContent) => {
                                        return Err((
                                            StatusCode::BAD_REQUEST,
                                            Json(create_safe_error_response(
                                                "invalid_content",
                                                "File content appears malformed or corrupted",
                                                None,
                                            )),
                                        )
                                            .into_response());
                                    }
                                    Err(_) => {
                                        return Err((
                                            StatusCode::BAD_REQUEST,
                                            Json(create_safe_error_response(
                                                "validation_failed",
                                                "File validation failed",
                                                None,
                                            )),
                                        )
                                            .into_response());
                                    }
                                }
                            }
                            Err(_) => had_parse_error = true,
                        }
                    }
                    "header_text" => match field.text().await {
                        Ok(text) => {
                            // Validate text field size
                            if text.len() > MAX_TEXT_FIELD_SIZE {
                                return Err((
                                    StatusCode::PAYLOAD_TOO_LARGE,
                                    Json(ErrorResponse {
                                        error: "Text field size exceeds limit".to_string(),
                                        error_type: "text_too_large".to_string(),
                                        details: None,
                                    }),
                                )
                                    .into_response());
                            }

                            if !text.trim().is_empty() {
                                input_data.text_content = Some(text);
                            }
                        }
                        Err(_) => had_parse_error = true,
                    },
                    "score_threshold" => {
                        if let Ok(text) = field.text().await {
                            if let Ok(threshold) = text.parse::<f64>() {
                                config.score_threshold = threshold.clamp(0.0, 1.0);
                            }
                        }
                    }
                    "result_limit" => {
                        if let Ok(text) = field.text().await {
                            if let Ok(limit) = text.parse::<usize>() {
                                config.result_limit = limit.clamp(1, 50); // Reasonable limits
                            }
                        }
                    }
                    "scoring_weights" => {
                        if let Ok(text) = field.text().await {
                            if let Ok(weights) = serde_json::from_str::<HashMap<String, f64>>(&text)
                            {
                                config.scoring_weights = parse_scoring_weights(&weights);
                            }
                        }
                    }
                    _ => {} // Ignore unknown fields
                }
            }
            Ok(None) => break, // No more fields
            Err(_) => {
                had_parse_error = true;
                break;
            }
        }
    }

    // Validate that we have some input
    if input_data.text_content.is_none() && input_data.binary_content.is_none() {
        let error_msg = if had_parse_error {
            "Failed to parse upload. Please check the file format."
        } else if fields_received == 0 {
            "No data received. Please upload a file or paste header text."
        } else {
            "No valid header data found in upload."
        };

        return Err((
            StatusCode::BAD_REQUEST,
            Json(create_safe_error_response(
                "missing_input",
                error_msg,
                None, // Never include details for consistency
            )),
        )
            .into_response());
    }

    Ok((input_data, config))
}

/// Parse input data using intelligent format detection
fn parse_input_data(
    input_data: &InputData,
) -> Result<crate::core::header::QueryHeader, Box<Response>> {
    if let Some(text_content) = &input_data.text_content {
        // Text-based parsing with format detection
        let Ok(detected_format) = detect_format(text_content, input_data.filename.as_deref())
        else {
            return Err(Box::new(
                (
                    StatusCode::BAD_REQUEST,
                    Json(create_safe_error_response(
                        "format_detection_failed",
                        "Unable to detect file format. Please check the file type and try again.",
                        Some("Format detection failed during parsing"),
                    )),
                )
                    .into_response(),
            ));
        };

        match parse_with_format(text_content, detected_format) {
            Ok(query) => Ok(query),
            Err(_) => Err(Box::new((
                StatusCode::BAD_REQUEST,
                Json(create_safe_error_response(
                    "parse_failed",
                    "Unable to process file content. Please check the file format and try again.",
                    Some("File parsing failed during content processing"),
                )),
            )
                .into_response())),
        }
    } else if let Some(binary_content) = &input_data.binary_content {
        // Binary file parsing
        let format = input_data.format.unwrap_or(FileFormat::Bam);

        match parse_binary_file(binary_content, format) {
            Ok(query) => Ok(query),
            Err(_) => Err(Box::new((
                StatusCode::BAD_REQUEST,
                Json(create_safe_error_response(
                    "binary_parse_failed",
                    "Unable to process binary file. Please verify the file format and try again.",
                    Some("Binary file parsing failed during processing"),
                )),
            )
                .into_response())),
        }
    } else {
        Err(Box::new(
            (
                StatusCode::INTERNAL_SERVER_ERROR,
                Json(ErrorResponse {
                    error: "Internal error: no input data".to_string(),
                    error_type: "internal_error".to_string(),
                    details: None,
                }),
            )
                .into_response(),
        ))
    }
}

/// Check if content appears to be binary
fn is_binary_content(bytes: &[u8]) -> bool {
    // Simple heuristic: if more than 1% of first 1024 bytes are non-printable, consider binary
    let sample_size = std::cmp::min(bytes.len(), 1024);

    // For very small samples, use a minimum threshold to avoid false positives
    if sample_size < 10 {
        return false; // Assume text for very small samples
    }

    let non_printable_count = bytes[..sample_size]
        .iter()
        .filter(|&&b| b < 9 || (b > 13 && b < 32) || b > 126)
        .count();

    // Use floating-point math to maintain consistent 1% threshold
    count_to_f64(non_printable_count) > (count_to_f64(sample_size) * 0.01)
}

/// Detect binary format from filename
fn detect_binary_format(filename: &str) -> Option<FileFormat> {
    let lower = filename.to_lowercase();
    if std::path::Path::new(&lower)
        .extension()
        .is_some_and(|ext| ext.eq_ignore_ascii_case("bam"))
    {
        Some(FileFormat::Bam)
    } else if std::path::Path::new(&lower)
        .extension()
        .is_some_and(|ext| ext.eq_ignore_ascii_case("cram"))
    {
        Some(FileFormat::Cram)
    } else {
        None
    }
}

/// Parse scoring weights from frontend format
fn parse_scoring_weights(weights: &HashMap<String, f64>) -> ScoringWeights {
    // Note: The frontend sends percentages (0-100), but the backend expects ratios (0-1)
    // New scoring model: contig_match, coverage, order, and conflict_penalty
    let contig_match = weights.get("contigMatch").unwrap_or(&70.0) / 100.0;
    let coverage = weights.get("coverage").unwrap_or(&20.0) / 100.0;
    let order = weights.get("orderScore").unwrap_or(&10.0) / 100.0;
    // Conflict penalty is a multiplier (0-1), not a weight percentage
    let conflict_penalty = weights.get("conflictPenalty").unwrap_or(&10.0) / 100.0;

    ScoringWeights {
        contig_match,
        coverage,
        order,
        conflict_penalty,
    }
}

/// Return list of references in catalog
async fn catalog_handler(State(state): State<Arc<AppState>>) -> Json<serde_json::Value> {
    let refs: Vec<serde_json::Value> = state
        .catalog
        .references
        .iter()
        .map(|r| {
            serde_json::json!({
                "id": r.id.0,
                "display_name": r.display_name,
                "assembly": format!("{}", r.assembly),
                "source": format!("{}", r.source),
                "contig_count": r.contigs.len(),
                "has_decoy": r.has_decoy(),
                "has_alt": r.has_alt(),
                "tags": r.tags,
            })
        })
        .collect();

    Json(serde_json::json!({
        "count": refs.len(),
        "references": refs,
    }))
}
