#!/bin/bash

set -eo pipefail

# Visium HD (FFPE, probe-based) test fixture from nf-core/test-datasets (human
# lung cancer post-Xenium HD FFPE). Its reads are curated to a tiny GRCh38
# reference, so Space Ranger's Visium HD pipeline runs to completion. Used for
# the spaceranger_hd_mapping integration test.
# Source: https://github.com/nf-core/test-datasets/tree/spatialvi/testdata/human-lung-cancer-post-xenium_hd_ffpe

REPO_ROOT=$(git rev-parse --show-toplevel)
DIR="$REPO_ROOT/resources_test/visium_hd_ffpe"
# Pin to a commit so the fixture stays reproducible if the branch moves.
COMMIT="2fe24b16442c450d289a2bd65b3b579bb1e16010"
BASE="https://raw.githubusercontent.com/nf-core/test-datasets/$COMMIT/testdata/human-lung-cancer-post-xenium_hd_ffpe"
REF_URL="https://raw.githubusercontent.com/nf-core/test-datasets/$COMMIT/testdata/GRCh38.tar.gz"
PRE="Visium_HD_Human_Lung_Cancer_HD_Only_Experiment2"

mkdir -p "$DIR/${PRE}_fastqs" "$DIR/GRCh38"

# FASTQs (2 lanes), already following the 10x naming convention
for lane in L001 L002; do
  for read in R1 R2; do
    curl -fSL -o "$DIR/${PRE}_fastqs/${PRE}_S1_${lane}_${read}_001.fastq.gz" \
      "$BASE/${PRE}_S1_${lane}_${read}_001.fastq.gz"
  done
done

# Probe set (full human transcriptome panel)
curl -fSL -o "$DIR/${PRE}_probe_set.csv" "$BASE/${PRE}_probe_set.csv"

# CytAssist brightfield image (--cytaimage) and H&E microscope image (--image).
# nf-core's samplesheet maps image.tif -> cytaimage and tissue_image.btf -> image.
curl -fSL -o "$DIR/${PRE}_cytaimage.tif" "$BASE/${PRE}_image.tif"
curl -fSL -o "$DIR/${PRE}_image.btf" "$BASE/${PRE}_tissue_image.btf"

# Tiny GRCh38 reference (built with spaceranger mkref; a handful of genes matched
# to the reads above)
curl -fSL -o "$DIR/GRCh38.tar.gz" "$REF_URL"
tar xzf "$DIR/GRCh38.tar.gz" -C "$DIR/GRCh38"
rm -f "$DIR/GRCh38.tar.gz"

# Sync to S3 (dry-run; drop --dryrun to upload)
aws s3 sync \
    --profile di \
    "$DIR" \
    s3://openpipelines-bio/openpipeline_spatial/resources_test/visium_hd_ffpe \
    --delete \
    --dryrun
