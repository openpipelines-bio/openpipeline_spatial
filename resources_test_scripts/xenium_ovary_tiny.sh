#!/bin/bash

set -eo pipefail

# A denser Xenium tiny test fixture: converts 10x's own "Xenium V1 Human Ovary tiny"
# XOA v4.0 example dataset (their site: "artificially subset to three 640 pixel
# square patches across two FOVs") and crops it down to one dense patch.
#
# Why a different source and method than xenium_tiny.sh:
#  * xenium_tiny.sh's source (nf-core "Xenium_Prime_Mouse_Ileum_tiny_outs") only has
#    23 cells -- too thin to meaningfully test a segmentation tool.
#  * This script converts the *full* Ovary dataset first (all elements,
#    `cells.zarr.zip` parsed exactly once by spatialdata_io itself), then crops the
#    resulting SpatialData object with spatialdata's own bounding_box_query() --
#    see crop_xenium_to_patch.py for why that sidesteps the cells.zarr.zip problem
#    entirely. The result keeps everything: image, raster labels, boundary shapes,
#    transcripts, and the cell annotation table, all consistently cropped together.
#    ~300 cells in one patch vs. 23 in the whole old fixture.
#
# IMPORTANT: this script does NOT sync anything to the shared test-data bucket.
# Validate the output locally against downstream components first.

REPO_ROOT=$(git rev-parse --show-toplevel)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

DIR="$REPO_ROOT/resources_test/xenium"
ID="xenium_ovary_tiny"

MY_TEMP="${VIASH_TEMP:-/tmp}"
TMPDIR=$(mktemp -d "$MY_TEMP/$ID-XXXXXX")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT

# 1. fetch the source dataset (10x's own XOA v4.0 example data; no top-level
#    directory in the zip, unlike the nf-core tarball xenium_tiny.sh uses)
SRC_OUTS="$TMPDIR/Xenium_V1_Human_Ovary_tiny_outs"
mkdir -p "$SRC_OUTS"
curl -fSL -o "$TMPDIR/xenium_ovary_tiny.zip" \
    "https://cf.10xgenomics.com/samples/xenium/4.0.0/Xenium_V1_Human_Ovary_tiny/Xenium_V1_Human_Ovary_tiny_outs.zip"
unzip -q "$TMPDIR/xenium_ovary_tiny.zip" -d "$SRC_OUTS"

# 2. convert the *full* dataset (all cells, all three patches, full imaging FOV) --
#    every element enabled (the converter's defaults), so cells.zarr.zip is read
#    once, correctly, by spatialdata_io itself
viash run "$REPO_ROOT/src/convert/from_xenium_to_spatialdata/config.vsh.yaml" -- \
    --input "$SRC_OUTS" \
    --output "$TMPDIR/full.zarr"

# 3. crop to one dense patch (the largest, by default) with bounding_box_query.
#    Run inside the *same* pinned Docker image as step 2, not the host's local
#    Python -- writing the cropped SpatialData object needs a matching
#    spatialdata/anndata/pandas stack to what components will later read it with.
#    A local venv with newer pandas/anndata than this project's pins (e.g. pandas
#    3.x's new default nullable-string dtype for the AnnData table's obs index) can
#    silently write a fixture that a pinned-anndata component then fails to re-write.
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

echo ""
echo "Done. Output written locally to:"
echo "  $DIR/$ID.zarr"
echo ""
echo "This script does NOT sync to S3. Once you've validated the output against"
echo "downstream components, sync it manually (drop --dryrun to actually upload):"
echo "  aws s3 sync --profile di \"$DIR\" s3://openpipelines-bio/openpipeline_spatial/resources_test/xenium --dryrun"
