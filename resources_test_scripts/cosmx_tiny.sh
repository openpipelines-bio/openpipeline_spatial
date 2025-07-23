#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

DIR="resources_test/cosmx"
ID="Lung5_Rep2"
OUT="$DIR/$ID/"

# create tempdir
MY_TEMP="${VIASH_TEMP:-/tmp}"
TMPDIR=$(mktemp -d "$MY_TEMP/$ID-XXXXXX")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT

if [ ! -d "$OUT" ]; then
    flat_dataset="https://nanostring-public-share.s3.us-west-2.amazonaws.com/SMI-Compressed/Lung5_Rep2/Lung5_Rep2+SMI+Flat+data.tar.gz"
    wget  "$flat_dataset" -O "$TMPDIR/Lung5_Rep2.tar.gz"
    mkdir -p "$TMPDIR/Lung5_Rep2"
    tar -xzf "$TMPDIR/Lung5_Rep2.tar.gz" -C "$TMPDIR/Lung5_Rep2"
    mkdir -p "$OUT"
    mv "$TMPDIR/Lung5_Rep2/Lung5_Rep2/Lung5_Rep2-Flat_files_and_images/"* "$OUT/"
fi

viash run src/filter/subset_cosmx/config.vsh.yaml -- \
    --input "$OUT" \
    --num_fovs 3 \
    --subset_transcripts_file True \
    --subset_polygons_file False \
    --output "${DIR}/${ID}_tiny"

viash run src/convert/from_cosmx_to_h5mu/config.vsh.yaml -- \
    --input ${DIR}/${ID}_tiny \
    --output "$DIR/${ID}_tiny.h5mu" \
    --output_compression "gzip"

rm -rf "$OUT"

# Sync to S3
aws s3 sync \
    --profile di \
    "$DIR" \
    s3://openpipelines-bio/openpipeline_spatial/resources_test/cosmx \
    --delete \
    --dryrun
