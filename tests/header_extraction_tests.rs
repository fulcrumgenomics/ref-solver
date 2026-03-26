//! Integration tests for header-only parsing from binary uploads.
//!
//! Verifies that the server can parse headers from BAM/CRAM bytes
//! sent as multipart uploads, including truncated streams.

use std::io::Cursor;
use std::num::NonZeroUsize;

use noodles::bam;
use noodles::sam;
use noodles::sam::header::record::value::map::ReferenceSequence;
use noodles::sam::header::record::value::Map;

/// Create a minimal BAM file in memory with the given contigs.
fn create_test_bam(contigs: &[(&str, usize)]) -> Vec<u8> {
    let mut header = sam::Header::builder();
    for &(name, length) in contigs {
        header = header.add_reference_sequence(
            name,
            Map::<ReferenceSequence>::new(NonZeroUsize::new(length).unwrap()),
        );
    }
    let header = header.build();

    let mut buf = Vec::new();
    {
        let mut writer = bam::io::Writer::new(&mut buf);
        writer.write_header(&header).unwrap();
    }
    buf
}

#[test]
fn test_parse_bam_from_reader_roundtrip() {
    use ref_solver::parsing::sam::parse_bam_from_reader;

    let contigs = [
        ("chr1", 248_956_422),
        ("chr2", 242_193_529),
        ("chr3", 198_295_559),
        ("chrX", 156_040_895),
        ("chrY", 57_227_415),
        ("chrM", 16_569),
    ];

    let bam_bytes = create_test_bam(&contigs);
    let cursor = Cursor::new(&bam_bytes);
    let query = parse_bam_from_reader(cursor).unwrap();

    assert_eq!(query.contigs.len(), contigs.len());
    for (i, &(name, length)) in contigs.iter().enumerate() {
        assert_eq!(query.contigs[i].name, name);
        assert_eq!(query.contigs[i].length, length as u64);
    }
}

#[test]
fn test_parse_bam_from_reader_with_trailing_garbage() {
    use ref_solver::parsing::sam::parse_bam_from_reader;

    let mut bam_bytes = create_test_bam(&[("chr1", 248_956_422)]);
    // Simulate a truncated upload — header is complete but record data is garbage.
    // parse_bam_from_reader only reads the header, so this must still succeed.
    bam_bytes.extend_from_slice([0xDE, 0xAD, 0xBE, 0xEF].repeat(25).as_slice());
    let cursor = Cursor::new(&bam_bytes);
    let query = parse_bam_from_reader(cursor).unwrap();

    assert_eq!(query.contigs.len(), 1);
    assert_eq!(query.contigs[0].name, "chr1");
}
