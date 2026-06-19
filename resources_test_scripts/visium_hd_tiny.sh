#!/bin/bash

set -eo pipefail

# Visium HD tiny test fixture.
#  * The cropped SpaceRanger outs come from the official 10x "Visium HD Tiny
#    3prime Mouse Brain" developer dataset (subset_visium_hd.py crops it to a
#    small tissue-straddling region, since its bin grid spans the whole slide).
#  * The raw FASTQs and microscope image — needed to test the ingestion workflow
#    end to end — are subset from the full "Visium HD 3prime Mouse Brain" dataset
#    (the same sample), mirroring visium_tiny.sh.
# Source: https://www.10xgenomics.com/datasets

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

DIR="$REPO_ROOT/resources_test/visium_hd"
ID="Visium_HD_Mouse_Brain"
TINY_OUTS_URL="https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_Tiny_3prime_Dataset/Visium_HD_Tiny_3prime_Dataset_outs.zip"
FULL_BASE="https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_3prime_Mouse_Brain"

# create tempdir for the (large) public downloads
MY_TEMP="${VIASH_TEMP:-/tmp}"
TMPDIR=$(mktemp -d "$MY_TEMP/visium_hd_tiny-XXXXXX")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT

mkdir -p "$DIR"

# 1. cropped SpaceRanger outs from the official tiny dataset
curl -fSL -o "$TMPDIR/outs.zip" "$TINY_OUTS_URL"
unzip -q "$TMPDIR/outs.zip" -d "$TMPDIR/outs"
python3 "$SCRIPT_DIR/subset_visium_hd.py" \
    "$TMPDIR/outs" \
    "$DIR/${ID}_tiny_spaceranger" \
    "$ID"

# 2. tiny FASTQ run folder (first 100,000 reads) from the full dataset
curl -fSL -o "$TMPDIR/fastqs.tar" "$FULL_BASE/Visium_HD_3prime_Mouse_Brain_fastqs.tar"
mkdir -p "$TMPDIR/fastqs" "$DIR/${ID}_tiny"
tar -xf "$TMPDIR/fastqs.tar" -C "$TMPDIR/fastqs"
for r in R1 R2; do
  src=$(find "$TMPDIR/fastqs" -name "*_L001_${r}_001.fastq.gz" | sort | head -1)
  gzip -cdf "$src" | head -n 400000 | gzip -c \
    > "$DIR/${ID}_tiny/${ID}_S1_L001_${r}_001.fastq.gz"
done

# 3. downsized microscope image from the full dataset
curl -fSL -o "$TMPDIR/image.tif" "$FULL_BASE/Visium_HD_3prime_Mouse_Brain_image.tif"
convert "$TMPDIR/image.tif" -resize 2000x2000 "$DIR/${ID}_image_tiny.jpg"

# Sync to S3 (dry-run; drop --dryrun to upload)
aws s3 sync \
    --profile di \
    "$DIR" \
    s3://openpipelines-bio/openpipeline_spatial/resources_test/visium_hd \
    --delete \
    --dryrun
