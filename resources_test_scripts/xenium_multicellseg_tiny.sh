#!/bin/bash

set -eo pipefail

# A multimodal-cell-segmentation Xenium tiny test fixture: crops 10x's own
# "Xenium V1 MultiCellSeg Human Ovary tiny" XOA v4.0 example dataset (their site:
# "artificially subset to three 640 pixel square patches across two FOVs") down
# to one dense patch, kept in the native raw Xenium bundle format (like
# xenium_tiny.sh's fixture), plus a SpatialData-converted copy alongside it.
#
# Why this dataset, rather than xenium_tiny.sh's:
#  * xenium_tiny.sh's source (nf-core "Xenium_Prime_Mouse_Ileum_tiny_outs") and the
#    plain "Xenium V1 Human Ovary tiny" dataset are both segmented from nuclear
#    expansion alone -- a single DAPI channel. Neither exercises the multi-channel
#    (DAPI + protein/RNA boundary stains) multimodal cell segmentation path that
#    xeniumranger components need to test against.
#
# Cropping uses the `filter/subset_xenium` component, which crops the raw bundle
# directly (cells.parquet, cell/nucleus boundaries, transcripts, the
# cell_feature_matrix.h5 CellRanger matrix, the morphology_focus OME-TIFF
# channels, and the cells.zarr.zip raster labels + metadata table), rather than
# converting first and cropping the converted representation. It's run locally
# since it was added alongside this test fixture and isn't part of a release yet.
#
# The cropped bundle is then also converted to SpatialData, using the
# openpipeline_spatial v0.6.0 release published on Viash Hub
# (https://www.viash-hub.com/packages/openpipeline_spatial), the same way
# cosmx_tiny.sh converts its own cropped/subset bundle.

REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

DIR="resources_test/xenium"
ID="xenium_multicellseg_tiny"
OPENPIPELINE_SPATIAL_VERSION="v0.6.0"

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

# 2. crop to one dense patch (the largest, by default) with the local
#    filter/subset_xenium component, keeping the native raw bundle format.
rm -rf "$DIR/$ID"
viash run "$REPO_ROOT/src/filter/subset_xenium/config.vsh.yaml" -- \
    --input "$SRC_OUTS" \
    --output "$DIR/$ID"

echo "> Cropping complete"

# 3. also publish a SpatialData-converted copy of the cropped bundle, using the
#    openpipeline_spatial $OPENPIPELINE_SPATIAL_VERSION release published on
#    Viash Hub.
cat > "$TMPDIR/convert_params.yaml" <<HERE
param_list:
- id: $ID
  input: "$DIR/$ID"
  output: "$ID.zarr"
HERE

rm -rf "$DIR/$ID.zarr"
nextflow run https://packages.viash-hub.com/vsh/openpipeline_spatial.git \
  -revision "$OPENPIPELINE_SPATIAL_VERSION" \
  -main-script target/nextflow/convert/from_xenium_to_spatialdata/main.nf \
  -params-file "$TMPDIR/convert_params.yaml" \
  -profile docker \
  -resume \
  -c src/workflows/utils/labels_ci.config \
  --publish_dir "$DIR"

echo "> Conversion to SpatialData complete"

# Sync to S3 (dry-run; drop --dryrun to upload)
aws s3 sync \
    --profile di \
    "$DIR" \
    s3://openpipelines-bio/openpipeline_spatial/resources_test/xenium \
    --delete \
    --dryrun
