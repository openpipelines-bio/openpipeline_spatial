#!/bin/bash

set -eo pipefail

# Tiny mm10 Space Ranger reference for the Visium HD ingestion test, kept as a
# pre-staged artifact on S3 (same approach as reference_tiny.sh for GRCh38). This
# script syncs it down; it is not rebuilt here.
#
# Provenance (how the artifact on S3 was made): the official 10x "mm10-2020-A"
# reference was subset to chr19 (the smallest mouse autosome) + chrM, then the
# Space Ranger reference package was rebuilt with `spaceranger mkref --genome=mm10`
# (the genome name must stay "mm10" so it matches the mm10-2020-A mouse probe set,
# otherwise Space Ranger refuses to run). This keeps the STAR index small enough
# to run Space Ranger locally while still aligning a useful fraction of the reads.

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
