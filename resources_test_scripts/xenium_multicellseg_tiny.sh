#!/bin/bash

set -eo pipefail

# A multimodal-cell-segmentation Xenium tiny test fixture: converts 10x's own
# "Xenium V1 MultiCellSeg Human Ovary tiny" XOA v4.0 example dataset (their site:
# "artificially subset to three 640 pixel square patches across two FOVs") and
# crops it down to one dense patch.
#
# Why this dataset and method, rather than xenium_tiny.sh's:
#  * xenium_tiny.sh's source (nf-core "Xenium_Prime_Mouse_Ileum_tiny_outs") and the
#    plain "Xenium V1 Human Ovary tiny" dataset are both segmented from nuclear
#    expansion alone -- a single DAPI channel. Neither exercises the multi-channel
#    (DAPI + protein/RNA boundary stains) multimodal cell segmentation path that
#    xeniumranger components need to test against.
#  * This script converts the *full* MultiCellSeg Ovary dataset first (all elements,
#    `cells.zarr.zip` parsed exactly once by spatialdata_io itself), then crops the
#    resulting SpatialData object with spatialdata's own bounding_box_query() --
#    The result keeps everything: images (all morphology_focus channels),
#    raster labels, boundary shapes, transcripts, and the cell annotation table, all
#    consistently cropped together.


REPO_ROOT=$(git rev-parse --show-toplevel)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

DIR="$REPO_ROOT/resources_test/xenium"
ID="xenium_multicellseg_tiny"

MY_TEMP="${VIASH_TEMP:-/tmp}"
TMPDIR=$(mktemp -d "$MY_TEMP/$ID-XXXXXX")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT

mkdir -p "$DIR"

# 1. fetch the source dataset (10x's own XOA v4.0 example data)
SRC_OUTS="$TMPDIR/Xenium_V1_MultiCellSeg_Human_Ovary_tiny_outs"
mkdir -p "$SRC_OUTS"
curl -fSL -o "$TMPDIR/xenium_multicellseg_tiny.zip" \
    "https://cf.10xgenomics.com/samples/xenium/4.0.0/Xenium_V1_MultiCellSeg_Human_Ovary_tiny/Xenium_V1_MultiCellSeg_Human_Ovary_tiny_outs.zip"
unzip -q "$TMPDIR/xenium_multicellseg_tiny.zip" -d "$SRC_OUTS"

# 2. convert the *full* dataset (all cells, all three patches, full imaging FOV,
#    all four morphology_focus channels) -- every element enabled (the converter's
#    defaults)
viash run "$REPO_ROOT/src/convert/from_xenium_to_spatialdata/config.vsh.yaml" -- \
    --input "$SRC_OUTS" \
    --output "$TMPDIR/full.zarr"

# 3. crop to one dense patch (the largest, by default) with bounding_box_query.
CONVERTER_IMAGE="ghcr.io/openpipelines-bio/openpipeline_spatial/convert/from_xenium_to_spatialdata:latest"
rm -rf "$DIR/$ID.zarr"
docker run --rm \
    -v "$TMPDIR:$TMPDIR" \
    -v "$DIR:$DIR" \
    -v "$SCRIPT_DIR:$SCRIPT_DIR" \
    "$CONVERTER_IMAGE" \
    python3 "$SCRIPT_DIR/crop_xenium_to_patch.py" \
        --input "$TMPDIR/full.zarr" \
        --output "$DIR/$ID.zarr"

# Sync to S3 (dry-run; drop --dryrun to upload)
aws s3 sync \
    --profile di \
    "$DIR" \
    s3://openpipelines-bio/openpipeline_spatial/resources_test/xenium \
    --delete \
    --dryrun
