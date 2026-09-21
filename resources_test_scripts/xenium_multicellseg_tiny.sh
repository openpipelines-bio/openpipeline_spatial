#!/bin/bash

set -eo pipefail

# A multimodal-cell-segmentation Xenium tiny test fixture, meant to replace
# xenium_tiny.sh's fixture going forward: crops 10x's own "Xenium V1
# MultiCellSeg Human Ovary tiny" XOA v4.0 example dataset (their site:
# "artificially subset to three 640 pixel square patches across two FOVs")
# down to one dense patch, kept in the native raw Xenium bundle format (like
# xenium_tiny.sh's fixture), then runs the same processing suite xenium_tiny.sh
# does, so downstream components/tests can eventually switch over to it.
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
# channels, and the cells.zarr.zip raster labels + metadata table). It's run locally
# since it was added alongside this test fixture and isn't part of a release yet.

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

# 3. convert the cropped bundle to SpatialData, then to h5mu, using the
#    openpipeline_spatial $OPENPIPELINE_SPATIAL_VERSION release on Viash Hub.
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

cat > "$TMPDIR/to_h5mu_params.yaml" <<HERE
param_list:
- id: $ID
  input: "$DIR/$ID.zarr"
  output: "$ID.h5mu"
HERE

nextflow run https://packages.viash-hub.com/vsh/openpipeline_spatial.git \
  -revision "$OPENPIPELINE_SPATIAL_VERSION" \
  -main-script target/nextflow/convert/from_spatialdata_to_h5mu/main.nf \
  -params-file "$TMPDIR/to_h5mu_params.yaml" \
  -profile docker \
  -resume \
  -c src/workflows/utils/labels_ci.config \
  --publish_dir "$DIR"

echo "> Conversion to h5mu complete"

# 4. spatial neighborhood graph on the raw h5mu, in place.
cat > "$TMPDIR/neighbors_params.yaml" <<HERE
param_list:
- id: $ID
  input: "$DIR/$ID.h5mu"
  output: "$ID.h5mu"
HERE

nextflow run https://packages.viash-hub.com/vsh/openpipeline_spatial.git \
  -revision "$OPENPIPELINE_SPATIAL_VERSION" \
  -main-script target/nextflow/neighbors/spatial_neighborhood_graph/main.nf \
  -params-file "$TMPDIR/neighbors_params.yaml" \
  -profile docker \
  -resume \
  -c src/workflows/utils/labels_ci.config \
  --publish_dir "$DIR"

echo "> Spatial neighborhood graph complete"

# 5. QC workflow, from the openpipelines-bio/openpipeline repo (as in xenium_tiny.sh).
cat > "$TMPDIR/qc.yaml" <<HERE
param_list:
  - id: $ID
    input: "$DIR/$ID.h5mu"
var_name_mitochondrial_genes: mitochondrial
var_name_ribosomal_genes: ribosomal
output: '\$id.qc.h5mu'
output_compression: gzip
publish_dir: "$DIR"
HERE

nextflow run openpipelines-bio/openpipeline \
  -latest \
  -r 2.1.0 \
  -main-script target/nextflow/workflows/qc/qc/main.nf \
  -profile docker \
  -params-file "$TMPDIR/qc.yaml" \
  -resume \
  -config src/workflows/utils/labels_ci.config

echo "> QC complete"

# 6. spatial neighborhood graph on the QC'd h5mu.
cat > "$TMPDIR/neighbors_qc_params.yaml" <<HERE
param_list:
- id: $ID
  input: "$DIR/$ID.qc.h5mu"
  output: "$ID.qc.neighbors.h5mu"
HERE

nextflow run https://packages.viash-hub.com/vsh/openpipeline_spatial.git \
  -revision "$OPENPIPELINE_SPATIAL_VERSION" \
  -main-script target/nextflow/neighbors/spatial_neighborhood_graph/main.nf \
  -params-file "$TMPDIR/neighbors_qc_params.yaml" \
  -profile docker \
  -resume \
  -c src/workflows/utils/labels_ci.config \
  --publish_dir "$DIR"

echo "> Post-QC spatial neighborhood graph complete"

# 7. PCA, from the openpipelines-bio/openpipeline repo (as in xenium_tiny.sh).
cat > "$TMPDIR/pca.yaml" <<HERE
param_list:
  - id: $ID
    input: "$DIR/${ID}.qc.neighbors.h5mu"
output: '\$id.qc.neighbors.pca.h5mu'
output_compression: gzip
publish_dir: "$TMPDIR"
HERE

nextflow run openpipelines-bio/openpipeline \
  -latest \
  -r v4.2.0 \
  -main-script target/nextflow/dimred/pca/main.nf \
  -profile docker \
  -params-file "$TMPDIR/pca.yaml" \
  -config src/workflows/utils/labels_ci.config \
  -resume

echo "> PCA complete"

# 8. find_neighbors, from the openpipelines-bio/openpipeline repo (as in xenium_tiny.sh).
cat > "$TMPDIR/find_neighbors.yaml" <<HERE
param_list:
  - id: $ID
    input: "$TMPDIR/$ID.qc.neighbors.pca.h5mu"
output: '\$id.qc.all_neighbors.pca.h5mu'
output_compression: gzip
publish_dir: "$TMPDIR"
HERE

nextflow run openpipelines-bio/openpipeline \
  -r v4.2.0 \
  -main-script target/nextflow/neighbors/find_neighbors/main.nf \
  -profile docker \
  -params-file "$TMPDIR/find_neighbors.yaml" \
  -config src/workflows/utils/labels_ci.config \
  -resume

echo "> find_neighbors complete"

# Move the final output to the destination directory
mv "$TMPDIR/$ID.qc.all_neighbors.pca.h5mu" "$DIR/$ID.qc.all_neighbors.pca.h5mu"

# Sync to S3 (dry-run; drop --dryrun to upload)
aws s3 sync \
    --profile di \
    "$DIR" \
    s3://openpipelines-bio/openpipeline_spatial/resources_test/xenium \
    --delete \
    --dryrun
