#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

DIR="resources_test"
ID="xenium_tiny"
OUT="resources_test/$ID"

[ ! -d "$DIR" ] && mkdir -p "$DIR"

# tiny_dataset="https://raw.githubusercontent.com/nf-core/test-datasets/spatialxe/Xenium_Prime_Mouse_Ileum_tiny_outs.zip"
# wget "$tiny_dataset" -O "$DIR/xenium_tiny.zip"

unzip -j "$DIR/xenium_tiny.zip" -d "$OUT"
rm "$DIR/xenium_tiny.zip"