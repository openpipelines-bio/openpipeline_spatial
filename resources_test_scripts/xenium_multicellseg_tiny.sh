#!/bin/bash

set -eo pipefail

# A multimodal-cell-segmentation Xenium tiny test fixture, meant to replace
# xenium_tiny.sh's fixture going forward: 10x's own "Xenium V1 MultiCellSeg
# Human Ovary tiny" XOA v4.0 example dataset (their site: "artificially subset
# to three 640 pixel square patches across two FOVs"). Two raw Xenium bundles
# are stored under resources_test/xenium_multichannel/:
#  * xenium_multicellseg_tiny: the bundle exactly as downloaded (3 patches).
#  * xenium_multicellseg_tiny_cropped: the same bundle cropped to its largest
#    patch, with the local `filter/subset_xenium` component, for components
#    that need a smaller input (e.g. cellpose_sam).
#
# Why this dataset, rather than xenium_tiny.sh's:
#  * xenium_tiny.sh's source (nf-core "Xenium_Prime_Mouse_Ileum_tiny_outs") and the
#    plain "Xenium V1 Human Ovary tiny" dataset are both segmented from nuclear
#    expansion alone -- a single DAPI channel. Neither exercises the multi-channel
#    (DAPI + protein/RNA boundary stains) multimodal cell segmentation path that
#    xeniumranger components need to test against.
#
# Note that `subset_xenium` only works on bundles that already consist of
# separate patches of cells (as this dataset does), and leaves some files
# (e.g. metrics_summary.csv, analysis_summary.html) unchanged: see its
# description.

REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

DIR="resources_test/xenium_multichannel"
ID="xenium_multicellseg_tiny"
ID_CROPPED="${ID}_cropped"
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

# crop to one dense patch (the largest, by default) with the local
# filter/subset_xenium component, keeping the native raw bundle format.
rm -rf "$DIR/$ID_CROPPED"
viash run "$REPO_ROOT/src/filter/subset_xenium/config.vsh.yaml" -- \
    --input "$DIR/$ID" \
    --output "$DIR/$ID_CROPPED"

echo "> Cropping complete"

# Sync to S3 (dry-run; drop --dryrun to upload)
aws s3 sync \
    --profile di \
    "$DIR" \
    s3://openpipelines-bio/openpipeline_spatial/resources_test/xenium_multichannel \
    --delete \
    --dryrun
