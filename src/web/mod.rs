//! Web server for browser-based reference identification.
//!
//! This module provides an interactive web interface using Axum and HTMX.
//! Users can paste SAM headers or upload files to identify references.
//!
//! ## Starting the Server
//!
//! ```text
//! # Start on default port 8080
//! ref-finder serve
//!
//! # Custom port and auto-open browser
//! ref-finder serve --port 3000 --open
//!
//! # Bind to all interfaces
//! ref-finder serve --address 0.0.0.0
//! ```
//!
//! ## API Endpoints
//!
//! - `GET /` - Main page with header input form
//! - `POST /api/identify` - Identify reference from header (multipart form)
//! - `GET /api/catalog` - List all references in the catalog

pub mod format_detection;
pub mod server;
