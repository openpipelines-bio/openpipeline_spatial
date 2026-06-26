#!/bin/bash

set -eo pipefail

# Tiny mm10 Space Ranger reference for the Visium HD ingestion test, kept as a
# pre-staged artifact on S3 (same approach as reference_tiny.sh for GRCh38).
#
# The artifact was built by subsetting the official 10x "mm10-2020-A" reference
# to chr19 (smallest mouse autosome) + chrM and rebuilding the reference package
# with (genome name must stay "mm10" to match the mm10-2020-A mouse probe set):
#
#   spaceranger mkref --genome=mm10 --fasta=tiny.fa --genes=tiny.gtf --nthreads=4 --memgb=6
#
# TODO: replace this manual mkref step with a `viash run` of a spaceranger_mkref
# component (to be implemented), so the reference can be rebuilt here directly.

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"
DIR="resources_test/mm10"

mkdir -p $DIR

aws s3 sync \
    --profile di \
    s3://openpipelines-bio/openpipeline_spatial/resources_test/mm10 \
    "$DIR" \
    --delete \
    --dryrun
