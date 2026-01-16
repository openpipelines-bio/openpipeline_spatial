#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

DIR="resources_test/cosmx"

# create tempdir
MY_TEMP="${VIASH_TEMP:-/tmp}"
TMPDIR=$(mktemp -d "$MY_TEMP/cosmx-XXXXXX")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT

if [ ! -d "$DIR/cosmx" ]; then

    flat_dataset_rep_1="https://nanostring-public-share.s3.us-west-2.amazonaws.com/SMI-Compressed/Lung5_Rep1/Lung5_Rep1+SMI+Flat+data.tar.gz"
    wget  "$flat_dataset_rep_1" -O "$TMPDIR/Lung5_Rep1.tar.gz"
    mkdir -p "$TMPDIR/Lung5_Rep1"
    tar -xzf "$TMPDIR/Lung5_Rep1.tar.gz" -C "$TMPDIR/Lung5_Rep1"
    mkdir -p "$DIR/Lung5_Rep1/"
    mv "$TMPDIR/Lung5_Rep1/Lung5_Rep1/Lung5_Rep1-Flat_files_and_images/"* "$DIR/Lung5_Rep1/"

    flat_dataset_rep_2="https://nanostring-public-share.s3.us-west-2.amazonaws.com/SMI-Compressed/Lung5_Rep2/Lung5_Rep2+SMI+Flat+data.tar.gz"
    wget  "$flat_dataset_rep_2" -O "$TMPDIR/Lung5_Rep2.tar.gz"
    mkdir -p "$TMPDIR/Lung5_Rep2"
    tar -xzf "$TMPDIR/Lung5_Rep2.tar.gz" -C "$TMPDIR/Lung5_Rep2"
    mkdir -p "$DIR/Lung5_Rep2/"
    mv "$TMPDIR/Lung5_Rep2/Lung5_Rep2/Lung5_Rep2-Flat_files_and_images/"* "$DIR/Lung5_Rep2/"
fi

echo "> Downloading of datasets complete"

# Subset dataset to make it tiny
cat > /tmp/params.yaml << HERE
param_list:
- id: Lung5_Rep1
  input: "$DIR/Lung5_Rep1/"
  output : "Lung5_Rep1_tiny"
- id: Lung5_Rep2
  input: "$DIR/Lung5_Rep2/"
  output : "Lung5_Rep2_tiny"
HERE

nextflow run https://packages.viash-hub.com/vsh/openpipeline_spatial.git \
  -revision v0.1.1 \
  -main-script target/nextflow/filter/subset_cosmx/main.nf \
  -params-file params.yaml \
  -profile docker \
  -resume \
  -c src/workflows/utils/labels_ci.config \
  --publish_dir "$DIR" \
  --num_fovs 3 \
  --subset_transcripts_file True \
  --subset_polygons_file False  


echo "> Subsetting complete"

# Convert to h5mu
cat > /tmp/params.yaml << HERE
param_list:
- id: Lung5_Rep1
  input: "$DIR/Lung5_Rep1_tiny/"
  output : "Lung5_Rep1_tiny.h5mu"
- id: Lung5_Rep2
  input: "$DIR/Lung5_Rep2_tiny/"
  output : "Lung5_Rep2_tiny.h5mu"
HERE

nextflow run https://packages.viash-hub.com/vsh/openpipeline_spatial.git \
  -revision v0.1.1 \
  -main-script target/nextflow/convert/from_cosmx_to_h5mu/main.nf \
  -params-file /tmp/params.yaml \
  -profile docker \
  -resume \
  -c src/workflows/utils/labels_ci.config \
  --publish_dir "$DIR" \
  --output_compression "gzip"

echo "> Conversion to H5MU complete"


# Spatial Neighborhood Graph Calculation
cat > /tmp/params.yaml << HERE
param_list:
- id: Lung5_Rep1
  input: "$DIR/Lung5_Rep1_tiny.h5mu"
  output : "Lung5_Rep1_tiny.h5mu"
- id: Lung5_Rep2
  input: "$DIR/Lung5_Rep2_tiny.h5mu"
  output : "Lung5_Rep2_tiny.h5mu"
HERE

nextflow run https://packages.viash-hub.com/vsh/openpipeline_spatial.git \
  -revision v0.1.1 \
  -main-script target/nextflow/neighbors/spatial_neighborhood_graph/main.nf \
  -params-file /tmp/params.yaml \
  -profile docker \
  -resume \
  -c src/workflows/utils/labels_ci.config \
  --publish_dir "$DIR" \
  --output_compression "gzip"

echo "> Spatial neighbor graph calculation complete"

rm -rf "$DIR"/Lung5_Rep1/
rm -rf "$DIR"/Lung5_Rep2/

# Process datasets
cat > /tmp/params.yaml << HERE
param_list:
- id: Lung5_Rep1
  input: "$DIR/Lung5_Rep1_tiny.h5mu"
- id: Lung5_Rep2
  input: "$DIR/Lung5_Rep2_tiny.h5mu"
HERE

nextflow run https://packages.viash-hub.com/vsh/openpipeline_spatial.git \
  -revision v0.1.1 \
  -main-script target/nextflow/workflows/multiomics/spatial_process_samples/main.nf \
  -params-file /tmp/params.yaml \
  -profile docker \
  -resume \
  -c src/workflows/utils/labels_ci.config \
  --publish_dir "$DIR" \
  --output Lung5_tiny_processed.h5mu \
  --output_compression "gzip"

echo "> Sample processing complete"

# Sync to S3
aws s3 sync \
    --profile di \
    "$DIR" \
    s3://openpipelines-bio/openpipeline_spatial/resources_test/cosmx \
    --delete \
    --dryrun
