#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"
DIR="resources_test/visium"

# from https://www.10xgenomics.com/resources/datasets/human-ovarian-cancer-1-standard
mkdir -p $DIR && cd $DIR

# Input Files
curl -O https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Ovarian_Cancer/Visium_FFPE_Human_Ovarian_Cancer_fastqs.tar
curl -O https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Ovarian_Cancer/Visium_FFPE_Human_Ovarian_Cancer_image.jpg
curl -O https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Ovarian_Cancer/Visium_FFPE_Human_Ovarian_Cancer_probe_set.csv

# Extract
tar xvf Visium_FFPE_Human_Ovarian_Cancer_fastqs.tar

# Create subsampled dataset with ImageMagick
# https://imagemagick.org/index.php
mkdir subsampled
convert Visium_FFPE_Human_Ovarian_Cancer_image.jpg -resize 2000x2000 subsampled/Visium_FFPE_Human_Ovarian_Cancer_image.jpg
for f in Visium_FFPE_Human_Ovarian_Cancer_fastqs/*L001*R*; do gzip -cdf $f | head -n 40000 | gzip -c > subsampled/$(basename $f); done

aws s3 sync \
    --profile di \
    "$DIR" \
    s3://openpipelines-bio/openpipeline_spatial/resources_test/visium \
    --delete \
    --dryrun