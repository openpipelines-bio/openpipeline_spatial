#!/bin/bash

set -eo pipefail

## VIASH START
meta_executable="target/native/spaceranger/spaceranger_count/spaceranger_count"
meta_resources_dir="resources_test"
## VIASH END

test_data="$meta_resources_dir/visium"

echo "> Default test run"
"$meta_executable" \
    --id test_spaceranger \
    --transcriptome "$test_data/GRCh38" \
    --fastqs "$test_data/subsampled" \
    --probe_set "$test_data/Visium_FFPE_Human_Ovarian_Cancer_probe_set.csv" \
    --image "$test_data/subsampled/Visium_FFPE_Human_Ovarian_Cancer_image.jpg" \
    --unknown_slide visium-1 \
    --create_bam false

echo "> Checking outputs..."

# Define output directory
OUT_DIR="test_spaceranger/outs"

# Function to check if file exists and is non-empty
check_file() {
    local file=$1
    local description=$2
    echo -n "Checking $description... "
    if [ ! -f "$file" ]; then
        echo "FAIL (file not found)"
        exit 1
    elif [ ! -s "$file" ]; then
        echo "FAIL (file is empty)"
        exit 1
    else
        echo "OK"
    fi
}

# Check essential files
check_file "$OUT_DIR/web_summary.html" "web summary"
check_file "$OUT_DIR/metrics_summary.csv" "metrics summary"

echo "> All tests passed successfully!"