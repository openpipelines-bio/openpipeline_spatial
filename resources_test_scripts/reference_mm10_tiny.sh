#!/bin/bash

set -eo pipefail

# Tiny mm10 reference for the Visium HD ingestion test, built by subsetting the
# official 10x "mm10-2020-A" reference to a single small chromosome (chr19, the
# smallest mouse autosome) plus chrM, then rebuilding the Space Ranger reference
# package with `spaceranger mkref`. This keeps the STAR index small enough to run
# Space Ranger locally while still aligning a useful fraction of the mouse brain
# reads. The genome name MUST stay "mm10" so it matches the `reference_genome`
# field of the mm10-2020-A mouse probe set (see visium_hd_tiny.sh), otherwise
# Space Ranger refuses to run ("reference genome ... and probe set must be identical").
# Source: https://www.10xgenomics.com/support/software/space-ranger/downloads

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

DIR="$REPO_ROOT/resources_test/mm10"
GENOME="mm10"
CHROMS=("chr19" "chrM")
FULL_REF_URL="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz"
SPACERANGER_IMAGE="ghcr.io/data-intuitive/spaceranger:3.1"

# create tempdir for the (large) public download
MY_TEMP="${VIASH_TEMP:-/tmp}"
TMPDIR=$(mktemp -d "$MY_TEMP/reference_mm10_tiny-XXXXXX")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT

# 1. download and unpack the full mm10-2020-A reference
curl -fSL -o "$TMPDIR/mm10.tar.gz" "$FULL_REF_URL"
tar -xzf "$TMPDIR/mm10.tar.gz" -C "$TMPDIR"
FULL="$TMPDIR/refdata-gex-mm10-2020-A"

# 2. subset the genome FASTA and gene GTF to the selected chromosomes
chrom_re=$(printf "%s|" "${CHROMS[@]}")
chrom_re="${chrom_re%|}"
awk -v re="^($chrom_re)$" 'BEGIN{RS=">";ORS=""} NR>1{
    n=substr($0,1,index($0,"\n")-1); split(n,a," ");
    if(a[1] ~ re) print ">"$0
}' "$FULL/fasta/genome.fa" > "$TMPDIR/tiny.fa"
awk -F'\t' -v re="^($chrom_re)$" '/^#/ || $1 ~ re' "$FULL/genes/genes.gtf" > "$TMPDIR/tiny.gtf"

# 3. rebuild the Space Ranger reference package (keep genome name "mm10")
docker run --rm -v "$TMPDIR:/work" -w /work --entrypoint bash "$SPACERANGER_IMAGE" -c \
  "spaceranger mkref --genome=$GENOME --fasta=tiny.fa --genes=tiny.gtf --nthreads=4 --memgb=6"

rm -rf "$DIR"
mkdir -p "$(dirname "$DIR")"
mv "$TMPDIR/$GENOME" "$DIR"

# Sync to S3 (dry-run; drop --dryrun to upload)
aws s3 sync \
    --profile di \
    "$DIR" \
    s3://openpipelines-bio/openpipeline_spatial/resources_test/mm10 \
    --delete \
    --dryrun
