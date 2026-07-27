#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"
DIR="resources_test/GRCh38"

mkdir -p $DIR

aws s3 sync \
    --profile di \
    s3://openpipelines-bio/openpipeline_spatial/resources_test/GRCh38 \
    "$DIR" \
    --delete \
    --dryrun