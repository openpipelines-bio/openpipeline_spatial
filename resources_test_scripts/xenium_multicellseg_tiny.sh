#!/bin/bash

set -eo pipefail

# A multimodal-cell-segmentation Xenium tiny test fixture, meant to replace
# xenium_tiny.sh's fixture going forward: 10x's own "Xenium V1 MultiCellSeg
# Human Ovary tiny" XOA v4.0 example dataset (their site: "artificially subset
# to three 640 pixel square patches across two FOVs"), kept in the native raw
# Xenium bundle format exactly as downloaded (like xenium_tiny.sh's fixture),
# then run through the same processing suite xenium_tiny.sh does, so
# downstream components/tests can eventually switch over to it.
#
# Why this dataset, rather than xenium_tiny.sh's:
#  * xenium_tiny.sh's source (nf-core "Xenium_Prime_Mouse_Ileum_tiny_outs") and the
#    plain "Xenium V1 Human Ovary tiny" dataset are both segmented from nuclear
#    expansion alone -- a single DAPI channel. Neither exercises the multi-channel
#    (DAPI + protein/RNA boundary stains) multimodal cell segmentation path that
#    xeniumranger components need to test against.
#
# The bundle is deliberately not cropped: cropping keeps some files (e.g.
# metrics_summary.csv, analysis_summary.html) unchanged, so they would no
# longer match the rest of the bundle, and it relies on this dataset's patch
# layout. Components that need a smaller input can crop it themselves at test
# time with the `filter/subset_xenium` component.

REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

DIR="resources_test/xenium_multichannel"
ID="xenium_multicellseg_tiny"
OPENPIPELINE_SPATIAL_VERSION="v0.6.0"

MY_TEMP="${VIASH_TEMP:-/tmp}"
TMPDIR=$(mktemp -d "$MY_TEMP/$ID-XXXXXX")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT

mkdir -p "$DIR"

# fetch the source dataset (10x's own XOA v4.0 example data)
rm -rf "$DIR/$ID"
mkdir -p "$DIR/$ID"
curl -fSL -o "$TMPDIR/xenium_multicellseg_tiny.zip" \
    "https://cf.10xgenomics.com/samples/xenium/4.0.0/Xenium_V1_MultiCellSeg_Human_Ovary_tiny/Xenium_V1_MultiCellSeg_Human_Ovary_tiny_outs.zip"
unzip -q "$TMPDIR/xenium_multicellseg_tiny.zip" -d "$DIR/$ID"

echo "> Download complete"

# Sync to S3 (dry-run; drop --dryrun to upload)
aws s3 sync \
    --profile di \
    "$DIR" \
    s3://openpipelines-bio/openpipeline_spatial/resources_test/xenium_multichannel \
    --delete \
    --dryrun
