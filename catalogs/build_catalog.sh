#!/bin/bash
# Build all reference catalog entries
# Run from the repository root directory
#
# Each reference is built from three inputs:
#   1. Dict file — provides contig names, lengths, MD5, and SAM metadata (AS, UR, SP, AN)
#   2. FASTA file — provides sha512t24u digests (and MD5 for validation)
#   3. NCBI assembly report — provides sequence roles, aliases, and metadata
#
# Usage: ./catalogs/build_catalog.sh <refs_dir>
#   refs_dir: directory containing per-reference subdirectories with dict and FASTA files
#
# Expected directory layout under <refs_dir>:
#   Homo_sapiens_assembly38/{Homo_sapiens_assembly38.dict,Homo_sapiens_assembly38.fasta}
#   GRCh38_reference_genome/{GRCh38_full_analysis_set_plus_decoy_hla.dict,.fa}
#   hs38DH/{hs38DH.dict,hs38DH.fa}
#   hg38/{hg38.dict,hg38.fa}
#   hg38_p12_ucsc/{hg38.p12.dict,hg38.p12.fa}
#   GRCh38.p13/{GRCh38.p13.dict,GRCh38.p13.fasta}
#   grch38_dragen/{hg38_alt_aware_nohla.dict,.fa}
#   grch38_dragen_altmasked/{hg38.dict,hg38.fa}
#   grch38_gdc/{GRCh38.d1.vd1.dict,.fa}
#   Homo_sapiens_assembly19/{Homo_sapiens_assembly19.dict,.fasta}
#   hs37d5/{hs37d5.dict,hs37d5.fa}
#   hg19/{hg19.dict,hg19.fasta}
#   grch37_ncbi/{grch37_ncbi.dict,grch37_ncbi.fna}
#   chm13v2/{chm13v2.dict,chm13v2.fna}

set -e

REFS="${1:?Usage: $0 <refs_dir>}"
BINARY="./target/release/ref-solver"
SOURCES="catalogs/sources"
OUTPUT="$SOURCES"

# Assembly report URLs
GRCH38_REPORT="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt"
GRCH38_P13_REPORT="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt"
GRCH38_P12_REPORT="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.38_GRCh38.p12/GCF_000001405.38_GRCh38.p12_assembly_report.txt"
GRCH37_REPORT="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_assembly_report.txt"
CHM13_REPORT="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_assembly_report.txt"

echo "Building reference catalog entries..."

# GRCh38 references
echo "Building grch38_broad_analysis_set..."
$BINARY catalog build \
  --id grch38_broad_analysis_set \
  --name "GRCh38 Broad Analysis Set" \
  --assembly grch38 \
  --source broad \
  --input "$REFS/Homo_sapiens_assembly38/Homo_sapiens_assembly38.dict" \
  --input "$REFS/Homo_sapiens_assembly38/Homo_sapiens_assembly38.fasta" \
  --input "$SOURCES/GCF_000001405.40_GRCh38.p14_assembly_report.txt" \
  --download-url "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta" \
  --assembly-report-url "$GRCH38_REPORT" \
  --description "Broad Institute GRCh38 analysis set for GATK Best Practices" \
  --species "Homo sapiens" \
  --output "$OUTPUT/grch38_broad_analysis_set.json"

echo "Building grch38_1kg_analysis..."
$BINARY catalog build \
  --id grch38_1kg_analysis \
  --name "GRCh38 1000 Genomes Analysis Set" \
  --assembly grch38 \
  --source 1kg \
  --input "$REFS/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.dict" \
  --input "$REFS/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
  --input "$SOURCES/GCF_000001405.40_GRCh38.p14_assembly_report.txt" \
  --download-url "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
  --assembly-report-url "$GRCH38_REPORT" \
  --description "1000 Genomes GRCh38 analysis set with decoy and HLA" \
  --species "Homo sapiens" \
  --output "$OUTPUT/grch38_1kg_analysis.json"

echo "Building hs38DH..."
$BINARY catalog build \
  --id hs38DH \
  --name "hs38DH (GRCh38 + ALT + decoy + HLA)" \
  --assembly grch38 \
  --source 1kg \
  --input "$REFS/hs38DH/hs38DH.dict" \
  --input "$REFS/hs38DH/hs38DH.fa" \
  --input "$SOURCES/GCF_000001405.40_GRCh38.p14_assembly_report.txt" \
  --download-url "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
  --assembly-report-url "$GRCH38_REPORT" \
  --description "hs38DH from bwakit - GRCh38 with ALT, decoy, and HLA contigs" \
  --species "Homo sapiens" \
  --output "$OUTPUT/hs38DH.json"

echo "Building hg38_ucsc..."
$BINARY catalog build \
  --id hg38_ucsc \
  --name "hg38 (UCSC)" \
  --assembly grch38 \
  --source ucsc \
  --input "$REFS/hg38/hg38.dict" \
  --input "$REFS/hg38/hg38.fa" \
  --input "$SOURCES/GCF_000001405.40_GRCh38.p14_assembly_report.txt" \
  --download-url "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz" \
  --assembly-report-url "$GRCH38_REPORT" \
  --description "UCSC hg38 reference genome" \
  --species "Homo sapiens" \
  --output "$OUTPUT/hg38_ucsc.json"

echo "Building hg38_p12_ucsc..."
$BINARY catalog build \
  --id hg38_p12_ucsc \
  --name "hg38.p12 (UCSC)" \
  --assembly grch38 \
  --source ucsc \
  --input "$REFS/hg38_p12_ucsc/hg38.p12.dict" \
  --input "$REFS/hg38_p12_ucsc/hg38.p12.fa" \
  --input "$SOURCES/GCF_000001405.38_GRCh38.p12_assembly_report.txt" \
  --download-url "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p12/hg38.p12.fa.gz" \
  --assembly-report-url "$GRCH38_P12_REPORT" \
  --description "UCSC hg38.p12 reference genome (GRCh38.p12)" \
  --species "Homo sapiens" \
  --output "$OUTPUT/hg38_p12_ucsc.json"

echo "Building grch38_ncbi..."
$BINARY catalog build \
  --id grch38_ncbi \
  --name "GRCh38.p13 (NCBI)" \
  --assembly grch38 \
  --source ncbi \
  --input "$REFS/GRCh38.p13/GRCh38.p13.dict" \
  --input "$REFS/GRCh38.p13/GRCh38.p13.fasta" \
  --input "$SOURCES/GCF_000001405.39_GRCh38.p13_assembly_report.txt" \
  --download-url "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz" \
  --assembly-report-url "$GRCH38_P13_REPORT" \
  --description "NCBI GRCh38.p13 no-ALT analysis set" \
  --species "Homo sapiens" \
  --output "$OUTPUT/grch38_ncbi.json"

echo "Building grch38_dragen..."
$BINARY catalog build \
  --id grch38_dragen \
  --name "GRCh38 DRAGEN" \
  --assembly grch38 \
  --source illumina \
  --input "$REFS/grch38_dragen/hg38_alt_aware_nohla.dict" \
  --input "$REFS/grch38_dragen/hg38_alt_aware_nohla.fa" \
  --input "$SOURCES/GCF_000001405.40_GRCh38.p14_assembly_report.txt" \
  --download-url "https://ilmn-dragen-giab-samples.s3.amazonaws.com/FASTA/hg38_alt_aware_nohla.fa" \
  --assembly-report-url "$GRCH38_REPORT" \
  --description "Illumina DRAGEN GRCh38 ALT-aware (no HLA)" \
  --species "Homo sapiens" \
  --output "$OUTPUT/grch38_dragen.json"

echo "Building grch38_dragen_altmasked..."
$BINARY catalog build \
  --id grch38_dragen_altmasked \
  --name "GRCh38 DRAGEN ALT-masked" \
  --assembly grch38 \
  --source illumina \
  --input "$REFS/grch38_dragen_altmasked/hg38.dict" \
  --input "$REFS/grch38_dragen_altmasked/hg38.fa" \
  --input "$SOURCES/GCF_000001405.40_GRCh38.p14_assembly_report.txt" \
  --download-url "https://ilmn-dragen-giab-samples.s3.amazonaws.com/FASTA/hg38.fa" \
  --assembly-report-url "$GRCH38_REPORT" \
  --description "Illumina DRAGEN GRCh38 with ALT contigs masked" \
  --species "Homo sapiens" \
  --output "$OUTPUT/grch38_dragen_altmasked.json"

echo "Building grch38_gdc..."
$BINARY catalog build \
  --id grch38_gdc \
  --name "GRCh38.d1.vd1 (NCI GDC)" \
  --assembly grch38 \
  --source gdc \
  --input "$REFS/grch38_gdc/GRCh38.d1.vd1.dict" \
  --input "$REFS/grch38_gdc/GRCh38.d1.vd1.fa" \
  --input "$SOURCES/GCF_000001405.40_GRCh38.p14_assembly_report.txt" \
  --download-url "https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834" \
  --assembly-report-url "$GRCH38_REPORT" \
  --description "NCI Genomic Data Commons GRCh38 with decoys and viral sequences" \
  --species "Homo sapiens" \
  --output "$OUTPUT/grch38_gdc.json"

# GRCh37 references
echo "Building b37_broad..."
$BINARY catalog build \
  --id b37_broad \
  --name "b37 (Broad)" \
  --assembly grch37 \
  --source broad \
  --input "$REFS/Homo_sapiens_assembly19/Homo_sapiens_assembly19.dict" \
  --input "$REFS/Homo_sapiens_assembly19/Homo_sapiens_assembly19.fasta" \
  --input "$SOURCES/GCF_000001405.25_GRCh37.p13_assembly_report.txt" \
  --download-url "https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta" \
  --assembly-report-url "$GRCH37_REPORT" \
  --description "Broad Institute b37 (GRCh37-based) for GATK" \
  --species "Homo sapiens" \
  --output "$OUTPUT/b37_broad.json"

echo "Building hs37d5..."
$BINARY catalog build \
  --id hs37d5 \
  --name "hs37d5 (1000 Genomes + decoy)" \
  --assembly grch37 \
  --source 1kg \
  --input "$REFS/hs37d5/hs37d5.dict" \
  --input "$REFS/hs37d5/hs37d5.fa" \
  --input "$SOURCES/GCF_000001405.25_GRCh37.p13_assembly_report.txt" \
  --download-url "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz" \
  --assembly-report-url "$GRCH37_REPORT" \
  --description "1000 Genomes GRCh37 with decoy sequences" \
  --species "Homo sapiens" \
  --output "$OUTPUT/hs37d5.json"

echo "Building hg19_ucsc..."
$BINARY catalog build \
  --id hg19_ucsc \
  --name "hg19 (UCSC)" \
  --assembly grch37 \
  --source ucsc \
  --input "$REFS/hg19/hg19.dict" \
  --input "$REFS/hg19/hg19.fasta" \
  --input "$SOURCES/GCF_000001405.25_GRCh37.p13_assembly_report.txt" \
  --download-url "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz" \
  --assembly-report-url "$GRCH37_REPORT" \
  --description "UCSC hg19 reference genome" \
  --species "Homo sapiens" \
  --output "$OUTPUT/hg19_ucsc.json"

echo "Building grch37_ncbi..."
$BINARY catalog build \
  --id grch37_ncbi \
  --name "GRCh37.p13 (NCBI)" \
  --assembly grch37 \
  --source ncbi \
  --input "$REFS/grch37_ncbi/grch37_ncbi.dict" \
  --input "$REFS/grch37_ncbi/grch37_ncbi.fna" \
  --input "$SOURCES/GCF_000001405.25_GRCh37.p13_assembly_report.txt" \
  --download-url "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.fna.gz" \
  --assembly-report-url "$GRCH37_REPORT" \
  --description "NCBI GRCh37.p13 reference genome" \
  --species "Homo sapiens" \
  --output "$OUTPUT/grch37_ncbi.json"

# T2T-CHM13
echo "Building chm13v2..."
$BINARY catalog build \
  --id chm13v2 \
  --name "T2T-CHM13v2.0" \
  --assembly chm13 \
  --source t2t \
  --input "$REFS/chm13v2/chm13v2.dict" \
  --input "$REFS/chm13v2/chm13v2.fna" \
  --input "$SOURCES/GCF_009914755.1_T2T-CHM13v2.0_assembly_report.txt" \
  --download-url "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz" \
  --assembly-report-url "$CHM13_REPORT" \
  --description "Telomere-to-Telomere CHM13v2.0 complete human genome" \
  --species "Homo sapiens" \
  --output "$OUTPUT/chm13v2.json"

echo ""
echo "Merging into catalog..."
jq -s '{
  "version": "1.0.0",
  "created_at": "'$(date -u +%Y-%m-%dT%H:%M:%SZ)'",
  "references": .
}' \
  "$OUTPUT/grch38_broad_analysis_set.json" \
  "$OUTPUT/grch38_1kg_analysis.json" \
  "$OUTPUT/hs38DH.json" \
  "$OUTPUT/hg38_ucsc.json" \
  "$OUTPUT/hg38_p12_ucsc.json" \
  "$OUTPUT/grch38_ncbi.json" \
  "$OUTPUT/grch38_dragen.json" \
  "$OUTPUT/grch38_dragen_altmasked.json" \
  "$OUTPUT/grch38_gdc.json" \
  "$OUTPUT/b37_broad.json" \
  "$OUTPUT/hs37d5.json" \
  "$OUTPUT/hg19_ucsc.json" \
  "$OUTPUT/grch37_ncbi.json" \
  "$OUTPUT/chm13v2.json" \
  > catalogs/human_references.json

echo "Done! Catalog written to catalogs/human_references.json"
echo ""
echo "To embed in binary, run: cargo build --release"
