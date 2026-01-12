#!/bin/bash
# Download all GRCh37 and GRCh38 patch version FASTAs
# These are needed to generate .dict files with MD5 checksums

set -e

DEST="catalogs/sources/fasta/patches"
mkdir -p "$DEST"

echo "Downloading GRCh38 patch FASTAs..."

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
  url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.${ver}_${patch}/GCF_000001405.${ver}_${patch}_genomic.fna.gz"
  outfile="$DEST/GCF_000001405.${ver}_${patch}_genomic.fna.gz"

  if [ -f "$outfile" ]; then
    echo "Already have $patch, skipping..."
  else
    echo "Downloading $patch (~3GB)..."
    curl -L -o "$outfile" "$url"
  fi
done

echo ""
echo "Downloading GRCh37 patch FASTAs..."

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
  url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.${ver}_${patch}/GCF_000001405.${ver}_${patch}_genomic.fna.gz"
  outfile="$DEST/GCF_000001405.${ver}_${patch}_genomic.fna.gz"

  if [ -f "$outfile" ]; then
    echo "Already have $patch, skipping..."
  else
    echo "Downloading $patch (~3GB)..."
    curl -L -o "$outfile" "$url"
  fi
done

echo ""
echo "All FASTAs downloaded to $DEST"
echo ""
echo "Next step: Generate .dict files with:"
echo "  for f in $DEST/*.fna.gz; do"
echo "    samtools dict -o \${f%.fna.gz}.dict \$f"
echo "  done"
