#!/bin/bash

set -eo pipefail

# Visium HD tiny test fixture, cropped from the official 10x "Visium HD Tiny
# 3prime Mouse Brain" developer dataset (Space Ranger 4.0.1). That dataset is
# downsampled in reads and image, but its bin grid still spans the whole slide
# (millions of 2 um bins), so it is cropped to one small tissue-straddling region
# (see subset_visium_hd.py). Source:
# https://www.10xgenomics.com/support/software/space-ranger/latest/resources/visium-hd-example-data

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

DIR="$REPO_ROOT/resources_test/visium_hd"
ID="Visium_HD_Mouse_Brain"
URL="https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_Tiny_3prime_Dataset/Visium_HD_Tiny_3prime_Dataset_outs.zip"

# create tempdir for the (large) public outs zip
MY_TEMP="${VIASH_TEMP:-/tmp}"
TMPDIR=$(mktemp -d "$MY_TEMP/visium_hd_tiny-XXXXXX")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT

mkdir -p "$DIR"

# download and extract the official tiny SpaceRanger outs
curl -fSL -o "$TMPDIR/outs.zip" "$URL"
unzip -q "$TMPDIR/outs.zip" -d "$TMPDIR/outs"

# crop to a tiny tissue-straddling region across all bin sizes
python3 "$SCRIPT_DIR/subset_visium_hd.py" \
    "$TMPDIR/outs" \
    "$DIR/${ID}_tiny_spaceranger" \
    "$ID"

# Sync to S3 (dry-run; drop --dryrun to upload)
aws s3 sync \
    --profile di \
    "$DIR" \
    s3://openpipelines-bio/openpipeline_spatial/resources_test/visium_hd \
    --delete \
    --dryrun
