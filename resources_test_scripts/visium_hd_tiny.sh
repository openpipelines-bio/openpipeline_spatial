#!/bin/bash

set -eo pipefail

# Visium HD tiny test fixture, all derived from the official 10x "Visium HD Tiny
# 3prime Mouse Brain" developer dataset (a downsampled corner of the tissue):
#  * the cropped SpaceRanger outs (subset_visium_hd.py crops them further to a
#    small tissue-straddling region, since the bin grid spans the whole slide);
#  * the raw FASTQs, reconstructed from the outs BAM with 10x bamtofastq, so they
#    are consistent with the outs and let the ingestion workflow run Space Ranger;
#  * the microscope image, from the full "Visium HD 3prime Mouse Brain" dataset.
# Source: https://www.10xgenomics.com/datasets

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

DIR="$REPO_ROOT/resources_test/visium_hd"
ID="Visium_HD_Mouse_Brain"
TINY_OUTS_URL="https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_Tiny_3prime_Dataset/Visium_HD_Tiny_3prime_Dataset_outs.zip"
FULL_BASE="https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_3prime_Mouse_Brain"
BAMTOFASTQ_VERSION="1.4.1"

# create tempdir for the (large) public downloads
MY_TEMP="${VIASH_TEMP:-/tmp}"
TMPDIR=$(mktemp -d "$MY_TEMP/visium_hd_tiny-XXXXXX")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT

mkdir -p "$DIR"

# 1. cropped SpaceRanger outs from the official tiny dataset
curl -fSL -o "$TMPDIR/outs.zip" "$TINY_OUTS_URL"
unzip -q "$TMPDIR/outs.zip" -d "$TMPDIR/outs"
python3 "$SCRIPT_DIR/subset_visium_hd.py" \
    --input "$TMPDIR/outs" \
    --output "$DIR/${ID}_tiny_spaceranger" \
    --dataset-id "$ID"

# 2. tiny FASTQ run folder, reconstructed from the cropped corner's BAM with 10x
#    bamtofastq. This keeps the reads consistent with the SpaceRanger outs above
#    and avoids downloading the full (~12 GB) FASTQ tar.
case "$(uname -s)" in
  Darwin) bamtofastq_asset="bamtofastq_macos" ;;
  *) bamtofastq_asset="bamtofastq_linux" ;;
esac
curl -fSL -o "$TMPDIR/bamtofastq" \
  "https://github.com/10XGenomics/bamtofastq/releases/download/v${BAMTOFASTQ_VERSION}/${bamtofastq_asset}"
chmod +x "$TMPDIR/bamtofastq"
"$TMPDIR/bamtofastq" --nthreads=4 \
  "$TMPDIR/outs/possorted_genome_bam.bam" "$TMPDIR/bamfq"
mkdir -p "$DIR/${ID}_tiny"
for r in R1 R2; do
  src=$(find "$TMPDIR/bamfq" -name "*_${r}_001.fastq.gz" | sort | head -1)
  cp "$src" "$DIR/${ID}_tiny/${ID}_S1_L001_${r}_001.fastq.gz"
done

# 3. downsized microscope image from the full dataset
curl -fSL -o "$TMPDIR/image.tif" "$FULL_BASE/Visium_HD_3prime_Mouse_Brain_image.tif"
convert "$TMPDIR/image.tif" -resize 2000x2000 "$DIR/${ID}_image_tiny.jpg"

# Sync to S3 (dry-run; drop --dryrun to upload)
aws s3 sync \
    --profile di \
    "$DIR" \
    s3://openpipelines-bio/openpipeline_spatial/resources_test/visium_hd \
    --delete \
    --dryrun
