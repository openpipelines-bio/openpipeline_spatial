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
    tiny_dataset="https://raw.githubusercontent.com/nf-core/test-datasets/spatialxe/Xenium_Prime_Mouse_Ileum_tiny_outs.zip"
    wget "$tiny_dataset" -O "$TMPDIR/xenium_tiny.zip"

    unzip -q "$TMPDIR/xenium_tiny.zip" -d "$TMPDIR/xenium_tiny"
    mkdir -p "$OUT"
    mv "$TMPDIR/xenium_tiny/Xenium_Prime_Mouse_Ileum_tiny_outs/"* "$OUT/"
fi

viash run "$REPO_ROOT/src/convert/from_xenium_to_spatialdata/config.vsh.yaml" -- \
    --input "$OUT" \
    --output "$DIR/$ID.zarr"

viash run "$REPO_ROOT/src/convert/from_spatialdata_to_h5mu/config.vsh.yaml" -- \
    --input "$DIR/$ID.zarr" \
    --output "$DIR/$ID.h5mu"

# Sync to S3
aws s3 sync \
    --profile di \
    "$DIR" \
    s3://openpipelines-bio/openpipeline_spatial/resources_test/xenium \
    --delete \
    --dryrun