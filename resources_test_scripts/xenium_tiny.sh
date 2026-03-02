#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# Define absolute directory paths
DIR="$REPO_ROOT/resources_test/xenium"
ID="xenium_tiny"
OUT="$DIR/$ID"

# create tempdir
MY_TEMP="${VIASH_TEMP:-/tmp}"
TMPDIR=$(mktemp -d "$MY_TEMP/$ID-XXXXXX")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT

if [ ! -d "$OUT" ]; then
    tiny_dataset="https://raw.githubusercontent.com/nf-core/test-datasets/spatialxe/Xenium_Prime_Mouse_Ileum_tiny_outs.tar.gz"
    wget "$tiny_dataset" -O "$TMPDIR/xenium_tiny.tar.gz"

    mkdir -p "$TMPDIR/xenium_tiny"
    tar -xzf "$TMPDIR/xenium_tiny.tar.gz" -C "$TMPDIR/xenium_tiny"
    mkdir -p "$OUT"
    mv "$TMPDIR/xenium_tiny/Xenium_Prime_Mouse_Ileum_tiny_outs/"* "$OUT/"
fi

rm -rf "$DIR/$ID.zarr"
viash run "$REPO_ROOT/src/convert/from_xenium_to_spatialdata/config.vsh.yaml" -- \
    --input "$OUT" \
    --output "$DIR/$ID.zarr"

viash run "$REPO_ROOT/src/convert/from_spatialdata_to_h5mu/config.vsh.yaml" -- \
    --input "$DIR/$ID.zarr" \
    --output "$DIR/$ID.h5mu"

viash run "$REPO_ROOT/src/neighbors/spatial_neighborhood_graph/config.vsh.yaml" -- \
    --input "$DIR/$ID.h5mu" \
    --output "$DIR/${ID}_neighbors.h5mu"

# Run PCA via openpipeline on the existing xenium_tiny.qc.neighbors.h5mu,
# which already has the spatial neighborhood graph pre-computed.
cat > /tmp/pca.yaml <<EOF
param_list:
  - id: xenium_tiny
    input: "$DIR/${ID}_neighbors.h5mu"
output: '\$id.qc.neighbors.pca.h5mu'
output_compression: gzip
publish_dir: "$TMPDIR"
EOF

nextflow run openpipelines-bio/openpipeline \
  -latest \
  -r v4.0.3 \
  -main-script target/nextflow/dimred/pca/main.nf \
  -profile docker \
  -params-file /tmp/pca.yaml \
  -config src/workflows/utils/labels_ci.config \
  -resume

# Run find_neighbors to add expression connectivities to .obsp.
# The input already contains spatial_connectivities; find_neighbors adds
# connectivities (expression) and distances without overwriting spatial keys.
cat > /tmp/find_neighbors.yaml <<EOF
param_list:
  - id: xenium_tiny
    input: "$TMPDIR/xenium_tiny.qc.neighbors.pca.h5mu"
output: '\$id.qc.neighbors.h5mu'
output_compression: gzip
publish_dir: "$TMPDIR"
EOF

nextflow run openpipelines-bio/openpipeline \
  -latest \
  -r v4.0.3 \
  -main-script target/nextflow/neighbors/find_neighbors/main.nf \
  -profile docker \
  -params-file /tmp/find_neighbors.yaml \
  -config src/workflows/utils/labels_ci.config \
  -resume

# Move the final output to the destination directory
mv "$TMPDIR/xenium_tiny.qc.neighbors.h5mu" "$DIR/xenium_tiny.qc.neighbors.h5mu"

# Sync to S3
aws s3 sync \
    --profile di \
    "$DIR" \
    s3://openpipelines-bio/openpipeline_spatial/resources_test/xenium \
    --delete \
    --dryrun