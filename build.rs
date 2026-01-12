use std::path::Path;

fn main() {
    let catalog_path = Path::new("catalogs/human_references.json");
    validate_catalog_file(catalog_path);
    set_build_dependencies();
}

fn validate_catalog_file(catalog_path: &Path) {
    // Ensure catalog exists at build time
    assert!(
        catalog_path.exists(),
        "\n\nCATALOG BUILD ERROR: File not found\n\
         Path: {}\n\
         Please create the catalog file before building.\n",
        catalog_path.display()
    );

    // Read catalog file
    let catalog_contents = std::fs::read_to_string(catalog_path).unwrap_or_else(|e| {
        panic!(
            "\n\nCATALOG BUILD ERROR: Failed to read file\n\
             Path: {}\n\
             Error: {e}\n",
            catalog_path.display()
        );
    });

    // Parse and validate JSON
    let catalog: serde_json::Value = serde_json::from_str(&catalog_contents).unwrap_or_else(|e| {
        panic!(
            "\n\nCATALOG BUILD ERROR: Invalid JSON\n\
             Path: {}\n\
             Error: {e}\n\
             Hint: Check for missing commas, brackets, or invalid syntax.\n",
            catalog_path.display()
        );
    });

    validate_catalog_structure(&catalog);
}

fn validate_catalog_structure(catalog: &serde_json::Value) {
    // Validate structure
    assert!(
        catalog.is_object(),
        "\n\nCATALOG BUILD ERROR: Root must be a JSON object\n\
         Got: {catalog}\n"
    );

    let references = catalog.get("references").unwrap_or_else(|| {
        panic!(
            "\n\nCATALOG BUILD ERROR: Missing 'references' field\n\
             The catalog must have a top-level 'references' array.\n"
        );
    });

    let refs = references.as_array().unwrap_or_else(|| {
        panic!(
            "\n\nCATALOG BUILD ERROR: 'references' must be an array\n\
             Got: {references}\n"
        );
    });

    // Validate each reference
    let total_contigs = validate_references(refs);

    println!(
        "cargo:warning=Validated catalog: {} references, {total_contigs} total contigs",
        refs.len()
    );
}

fn validate_references(refs: &[serde_json::Value]) -> usize {
    let mut total_contigs = 0;

    for (i, reference) in refs.iter().enumerate() {
        let ref_id = reference
            .get("id")
            .and_then(|v| v.as_str())
            .unwrap_or("<unknown>");

        validate_reference_fields(reference, ref_id, i);
        total_contigs += validate_reference_contigs(reference, ref_id);
    }

    total_contigs
}

fn validate_reference_fields(reference: &serde_json::Value, ref_id: &str, index: usize) {
    assert!(
        reference.get("id").is_some(),
        "\n\nCATALOG BUILD ERROR: Reference at index {index} missing 'id' field\n"
    );
    assert!(
        reference.get("display_name").is_some(),
        "\n\nCATALOG BUILD ERROR: Reference '{ref_id}' (index {index}) missing 'display_name' field\n"
    );
    assert!(
        reference.get("contigs").is_some(),
        "\n\nCATALOG BUILD ERROR: Reference '{ref_id}' (index {index}) missing 'contigs' field\n"
    );
}

fn validate_reference_contigs(reference: &serde_json::Value, ref_id: &str) -> usize {
    if let Some(contigs) = reference.get("contigs").and_then(|c| c.as_array()) {
        for (j, contig) in contigs.iter().enumerate() {
            validate_contig_fields(contig, ref_id, j);
        }
        contigs.len()
    } else {
        0
    }
}

fn validate_contig_fields(contig: &serde_json::Value, ref_id: &str, index: usize) {
    let contig_name = contig
        .get("name")
        .and_then(|v| v.as_str())
        .unwrap_or("<unknown>");

    assert!(
        contig.get("name").is_some(),
        "\n\nCATALOG BUILD ERROR: Reference '{ref_id}' contig {index} missing 'name' field\n"
    );

    let length = contig.get("length");
    assert!(
        length.is_some(),
        "\n\nCATALOG BUILD ERROR: Reference '{ref_id}' contig '{contig_name}' (index {index}) missing 'length' field\n"
    );

    // Validate length is positive
    if let Some(len) = length.and_then(serde_json::Value::as_u64) {
        assert!(
            len > 0,
            "\n\nCATALOG BUILD ERROR: Reference '{ref_id}' contig '{contig_name}' has zero length\n\
             Contigs must have length > 0.\n"
        );
    }
}

fn set_build_dependencies() {
    // Tell cargo to rerun if catalog changes
    println!("cargo:rerun-if-changed=catalogs/human_references.json");

    // Tell cargo to rerun if build.rs changes
    println!("cargo:rerun-if-changed=build.rs");
}
