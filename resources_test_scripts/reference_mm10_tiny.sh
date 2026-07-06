#!/bin/bash

set -eo pipefail

# Tiny mm10 Space Ranger reference for the Visium HD ingestion test, kept as a
# pre-staged artifact on S3 (same approach as reference_tiny.sh for GRCh38).
#
# Source: the official 10x Genomics "mm10-2020-A" gene expression reference
#   https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz
#
# The artifact was built by subsetting that reference to chr19 (smallest mouse
# autosome) + chrM and rebuilding the reference package. The genome name must
# stay "mm10" to match the mm10-2020-A mouse probe set. Reproduction recipe:
#
#   # Download and unpack the full 10x mm10-2020-A reference
#   curl -fSL -o refdata-gex-mm10-2020-A.tar.gz \
#     https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz
#   tar xzf refdata-gex-mm10-2020-A.tar.gz
#
#   # Subset the genome FASTA to chr19 + chrM
#   samtools faidx refdata-gex-mm10-2020-A/fasta/genome.fa chr19 chrM > tiny.fa
#
#   # Subset the gene annotation to the same contigs (keep header lines)
#   awk '/^#/ || $1 == "chr19" || $1 == "chrM"' \
#     refdata-gex-mm10-2020-A/genes/genes.gtf > tiny.gtf
#
#   # Rebuild the reference package (produces ./mm10)
#   spaceranger mkref --genome=mm10 --fasta=tiny.fa --genes=tiny.gtf \
#     --nthreads=4 --memgb=6
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
