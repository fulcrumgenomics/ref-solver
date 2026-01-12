#!/bin/bash
# Build catalog entries for all GRCh37 and GRCh38 patch versions
# Run this after downloading FASTAs and generating .dict files

set -e

PATCH_DIR="catalogs/sources/fasta/patches"
REPORT_DIR="catalogs/sources/ncbi_reports"

echo "Generating .dict files for all patch FASTAs..."
for f in "$PATCH_DIR"/*.fna.gz; do
    dict="${f%.fna.gz}.dict"
    if [ -f "$dict" ]; then
        echo "Already have $(basename "$dict"), skipping..."
    else
        echo "Generating $(basename "$dict")..."
        samtools dict -o "$dict" "$f"
    fi
done

echo ""
echo "Building GRCh38 patch references..."

# GRCh38 versions (26-40)
declare -A GRCH38_PATCHES=(
    [26]="GRCh38"
    [27]="GRCh38.p1"
    [28]="GRCh38.p2"
    [29]="GRCh38.p3"
    [30]="GRCh38.p4"
    [31]="GRCh38.p5"
    [32]="GRCh38.p6"
    [33]="GRCh38.p7"
    [34]="GRCh38.p8"
    [35]="GRCh38.p9"
    [36]="GRCh38.p10"
    [37]="GRCh38.p11"
    [38]="GRCh38.p12"
    [39]="GRCh38.p13"
    [40]="GRCh38.p14"
)

for ver in "${!GRCH38_PATCHES[@]}"; do
    patch="${GRCH38_PATCHES[$ver]}"
    id="grch38_p${patch#GRCh38.p}"
    [ "$patch" = "GRCh38" ] && id="grch38_original"

    dict="$PATCH_DIR/GCF_000001405.${ver}_${patch}_genomic.dict"
    report="$REPORT_DIR/GCF_000001405.${ver}_${patch}_assembly_report.txt"

    if [ ! -f "$dict" ]; then
        echo "Missing dict for $patch, skipping..."
        continue
    fi

    echo "Building $id ($patch)..."

    cargo run --release --quiet -- catalog build \
        --id "$id" \
        --name "$patch (NCBI RefSeq)" \
        --input "$dict" \
        --input "$report" \
        --assembly GRCh38 \
        --source NCBI \
        --download-url "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.${ver}_${patch}/GCF_000001405.${ver}_${patch}_genomic.fna.gz" \
        --assembly-report-url "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.${ver}_${patch}/GCF_000001405.${ver}_${patch}_assembly_report.txt" \
        --output "catalogs/patches/${id}.json"
done

echo ""
echo "Building GRCh37 patch references..."

# GRCh37 versions
declare -A GRCH37_PATCHES=(
    [13]="GRCh37"
    [14]="GRCh37.p2"
    [17]="GRCh37.p5"
    [21]="GRCh37.p9"
    [22]="GRCh37.p10"
    [23]="GRCh37.p11"
    [24]="GRCh37.p12"
    [25]="GRCh37.p13"
)

for ver in "${!GRCH37_PATCHES[@]}"; do
    patch="${GRCH37_PATCHES[$ver]}"
    id="grch37_p${patch#GRCh37.p}"
    [ "$patch" = "GRCh37" ] && id="grch37_original"

    dict="$PATCH_DIR/GCF_000001405.${ver}_${patch}_genomic.dict"
    report="$REPORT_DIR/GCF_000001405.${ver}_${patch}_assembly_report.txt"

    if [ ! -f "$dict" ]; then
        echo "Missing dict for $patch, skipping..."
        continue
    fi

    echo "Building $id ($patch)..."

    cargo run --release --quiet -- catalog build \
        --id "$id" \
        --name "$patch (NCBI RefSeq)" \
        --input "$dict" \
        --input "$report" \
        --assembly GRCh37 \
        --source NCBI \
        --download-url "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.${ver}_${patch}/GCF_000001405.${ver}_${patch}_genomic.fna.gz" \
        --assembly-report-url "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.${ver}_${patch}/GCF_000001405.${ver}_${patch}_assembly_report.txt" \
        --output "catalogs/patches/${id}.json"
done

echo ""
echo "Merging all patch references into catalog..."
mkdir -p catalogs/patches

# Merge all individual JSON files into the main catalog
for json in catalogs/patches/*.json; do
    if [ -f "$json" ]; then
        echo "Merging $(basename "$json")..."
        jq -s '.[0].references + [.[1]] | {version: "1.0", references: .}' \
            catalogs/human_references.json "$json" > catalogs/human_references.json.tmp
        mv catalogs/human_references.json.tmp catalogs/human_references.json
    fi
done

echo ""
echo "Done! Patch references added to catalog."
