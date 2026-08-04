#!/bin/bash

set -eo pipefail

# Visium HD (FFPE, probe-based) test fixture, taken from the nf-core/test-datasets
# "spatialvi" branch (human lung cancer post-Xenium HD FFPE). Unlike a naive tiny
# fixture, this dataset's reads are curated to match a tiny GRCh38 reference
# (only a handful of genes), so Space Ranger produces non-empty counts and its
# Visium HD pipeline runs to completion -- which a whole-transcriptome-reads vs
# tiny-reference mismatch does not. Used for the spaceranger_hd_mapping live
# integration test.
#
# Source: https://github.com/nf-core/test-datasets/tree/spatialvi/testdata/human-lung-cancer-post-xenium_hd_ffpe
#   reference: https://github.com/nf-core/test-datasets/blob/spatialvi/testdata/GRCh38.tar.gz
# nf-core runs Space Ranger HD on this via conf/test_spaceranger_hd.config.

REPO_ROOT=$(git rev-parse --show-toplevel)
DIR="$REPO_ROOT/resources_test/visium_hd_ffpe"
BASE="https://raw.githubusercontent.com/nf-core/test-datasets/spatialvi/testdata/human-lung-cancer-post-xenium_hd_ffpe"
REF_URL="https://raw.githubusercontent.com/nf-core/test-datasets/spatialvi/testdata/GRCh38.tar.gz"
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
