#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# Define absolute directory path
DIR="$REPO_ROOT/resources_test/visium"
ID="Visium_FFPE_Human_Ovarian_Cancer_tiny"

# from https://www.10xgenomics.com/resources/datasets/human-ovarian-cancer-1-standard
mkdir -p "$DIR"

# Input Files - download to the specific directory
curl -o "$DIR/Visium_FFPE_Human_Ovarian_Cancer_fastqs.tar" https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Ovarian_Cancer/Visium_FFPE_Human_Ovarian_Cancer_fastqs.tar
curl -o "$DIR/Visium_FFPE_Human_Ovarian_Cancer_image.jpg" https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Ovarian_Cancer/Visium_FFPE_Human_Ovarian_Cancer_image.jpg
curl -o "$DIR/Visium_FFPE_Human_Ovarian_Cancer_probe_set.csv" https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Ovarian_Cancer/Visium_FFPE_Human_Ovarian_Cancer_probe_set.csv

# Slide layout file, so spaceranger does not have to download it at runtime
curl -o "$DIR/V10L13-020.gpr" https://s3.us-west-2.amazonaws.com/10x.spatial-slides/gpr/V10L13/V10L13-020.gpr

# Extract in the specific directory
tar xvf "$DIR/Visium_FFPE_Human_Ovarian_Cancer_fastqs.tar" -C "$DIR"

# Create subsampled dataset with ImageMagick
# https://imagemagick.org/index.php
mkdir -p "$DIR/Visium_FFPE_Human_Ovarian_Cancer_tiny"
convert "$DIR/Visium_FFPE_Human_Ovarian_Cancer_image.jpg" -resize 2000x2000 "$DIR/Visium_FFPE_Human_Ovarian_Cancer_image_tiny.jpg"
for f in "$DIR"/Visium_FFPE_Human_Ovarian_Cancer_fastqs/*L001*R*; do 
  gzip -cdf "$f" | head -n 40000 | gzip -c > "$DIR/Visium_FFPE_Human_Ovarian_Cancer_tiny/$(basename "$f")"; 
done

echo "> Downloading and subsampling of datasets complete"

# Run spaceranger
viash run src/mapping/spaceranger_count/config.vsh.yaml -- \
  --input "$DIR/Visium_FFPE_Human_Ovarian_Cancer_tiny" \
  --gex_reference "$REPO_ROOT/resources_test/GRCh38/" \
  --probe_set "$DIR/Visium_FFPE_Human_Ovarian_Cancer_probe_set.csv" \
  --image "$DIR/Visium_FFPE_Human_Ovarian_Cancer_image_tiny.jpg" \
  --slide "V10L13-020" \
  --area "D1" \
  --slidefile "$DIR/V10L13-020.gpr" \
  --create_bam "false" \
  --output "$DIR/Visium_FFPE_Human_Ovarian_Cancer_tiny_spaceranger"

echo "> Running spaceranger complete"

rm -rf "$DIR/Visium_FFPE_Human_Ovarian_Cancer_fastqs"
rm -f "$DIR/Visium_FFPE_Human_Ovarian_Cancer_image.jpg"

viash run "$REPO_ROOT/src/convert/from_spaceranger_to_h5mu/config.vsh.yaml" -- \
    --input "$DIR/Visium_FFPE_Human_Ovarian_Cancer_tiny_spaceranger" \
    --output "$DIR/$ID.h5mu"

echo "> Conversion to H5MU complete"

# Process dataset (QC, normalization, PCA) for use as a downstream-workflow test fixture
nextflow run . \
  -main-script target/nextflow/workflows/multiomics/spatial_process_samples/main.nf \
  -profile docker \
  -resume \
  -c src/workflows/utils/labels_ci.config \
  --id "visium_tiny" \
  --input "$DIR/$ID.h5mu" \
  --output "visium_tiny_processed.h5mu" \
  --publish_dir "$DIR" \
  --output_compression "gzip"

echo "> Sample processing complete"

aws s3 sync \
    --profile di \
    --exclude "*.yaml" \
    "$DIR" \
    s3://openpipelines-bio/openpipeline_spatial/resources_test/visium \
    --delete \
    --dryrun