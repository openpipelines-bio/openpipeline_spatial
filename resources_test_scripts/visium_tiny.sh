#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# Define absolute directory path
DIR="$REPO_ROOT/resources_test/visium"

# from https://www.10xgenomics.com/resources/datasets/human-ovarian-cancer-1-standard
mkdir -p "$DIR"

# Input Files - download to the specific directory
curl -o "$DIR/Visium_FFPE_Human_Ovarian_Cancer_fastqs.tar" https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Ovarian_Cancer/Visium_FFPE_Human_Ovarian_Cancer_fastqs.tar
curl -o "$DIR/Visium_FFPE_Human_Ovarian_Cancer_image.jpg" https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Ovarian_Cancer/Visium_FFPE_Human_Ovarian_Cancer_image.jpg
curl -o "$DIR/Visium_FFPE_Human_Ovarian_Cancer_probe_set.csv" https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Ovarian_Cancer/Visium_FFPE_Human_Ovarian_Cancer_probe_set.csv

# Extract in the specific directory
tar xvf "$DIR/Visium_FFPE_Human_Ovarian_Cancer_fastqs.tar" -C "$DIR"

# Create subsampled dataset with ImageMagick
# https://imagemagick.org/index.php
mkdir -p "$DIR/subsampled"
convert "$DIR/Visium_FFPE_Human_Ovarian_Cancer_image.jpg" -resize 2000x2000 "$DIR/subsampled/Visium_FFPE_Human_Ovarian_Cancer_image.jpg"
for f in "$DIR"/Visium_FFPE_Human_Ovarian_Cancer_fastqs/*L001*R*; do 
  gzip -cdf "$f" | head -n 40000 | gzip -c > "$DIR/subsampled/$(basename "$f")"; 
done

aws s3 sync \
    --profile di \
    "$DIR" \
    s3://openpipelines-bio/openpipeline_spatial/resources_test/visium \
    --delete \
    --dryrun