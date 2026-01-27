#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

ID=aviti
DIR=resources_test/$ID/
OUT=$DIR/teton_cells2stats_tiny/

# Create directories
[ -d "$DIR" ] || mkdir -p "$DIR"
[ -d "$OUT" ] || mkdir -p "$OUT"

echo "> Downloading Aviti Teton data"
wget "https://go.elementbiosciences.com/l/938263/28kddnj7/d59cp" -O "${DIR}/PLUT-0105.tar.gz"
tar -xzf "${DIR}/PLUT-0105.tar.gz" -C "$DIR"
rm "${DIR}/PLUT-0105.tar.gz"

echo "> Processing and subsetting Aviti Teton data"
python <<HEREDOC
import os
import shutil
import pandas as pd
import glob
import json

src_dir = "${DIR}/PLUT-0105"
dest_dir = "${OUT}"
subset_image_dirs = False
wells_to_keep = ["A1"]
max_cells_per_well = 1000

os.makedirs(dest_dir, exist_ok=True)

print(f"Processing data from {src_dir} to {dest_dir}")

# Copy images
if subset_image_dirs:
    image_dirs = ["CellSegmentation", "Projection"]
    for image_dir in image_dirs:
        image_dir_path = os.path.join(src_dir, image_dir)
        if not os.path.exists(image_dir_path):
            print(f"Warning: Image directory not found: {image_dir_path}")
            continue
        if not os.path.isdir(image_dir_path):
            print(f"Warning: Path exists but is not a directory: {image_dir_path}")
            continue
        print(f"Processing image directory: {image_dir}")
        
        for well in wells_to_keep:
            dest_path = f"{dest_dir}/{image_dir}/Well{well}"
            os.makedirs(dest_path, exist_ok=True)
            src_path = glob.glob(os.path.join(src_dir, image_dir, f"Well{well}"))
            if len(src_path) != 1:
                print(f"Warning: Expected 1 path for Well{well}, found {len(src_path)}")
                continue
            shutil.copytree(src_path[0], os.path.join(dest_path), dirs_exist_ok=True)

# Copy count matrix
src_path = os.path.join(src_dir, "Cytoprofiling", "Instrument", "RawCellStats.parquet")
if os.path.exists(src_path):
    print(f"Processing count matrix: {src_path}")
    df = pd.read_parquet(src_path)
    print(f"Original data: {len(df)} rows")
    
    # Filter by wells
    df = df[df["Well"].isin(wells_to_keep)]
    print(f"After well filtering: {len(df)} rows")
    
    if max_cells_per_well:
        # Limit the number of cells per well
        df = df.head(max_cells_per_well)
        print(f"After cell limit: {len(df)} rows")

    dest_path = os.path.join(dest_dir, "Cytoprofiling", "Instrument")
    os.makedirs(dest_path, exist_ok=True)
    dest_file = os.path.join(dest_path, "RawCellStats.parquet")
    df.to_parquet(dest_file, engine="pyarrow")
    print(f"Saved processed count matrix to {dest_file}")
else:
    print(f"Warning: Count matrix not found at {src_path}")

# Copy Panel Metadata
panel_src_path = os.path.join(src_dir, "Panel.json")
if os.path.exists(panel_src_path):
    panel_dest_path = os.path.join(dest_dir, "Panel.json")
    shutil.copy2(panel_src_path, panel_dest_path)
    print(f"Copied Panel.json")
else:
    print(f"Warning: Panel.json not found at {panel_src_path}")
print("Processing complete!")
HEREDOC

echo "> Removing original aviti_teton folder"
rm -rf "$DIR/PLUT-0105"

echo "> Aviti Teton tiny dataset created successfully at $OUT"

viash run src/convert/from_cells2stats_to_h5mu/config.vsh.yaml -- \
    --input "${OUT}" \
    --output "$DIR/aviti_teton_tiny.h5mu" \
    --output_compression "gzip"

echo "> Conversion to H5MU complete"

viash run src/neighbors/spatial_neighborhood_graph/config.vsh.yaml -- \
    --input "$DIR/aviti_teton_tiny.h5mu" \
    --output "$DIR/aviti_teton_tiny.h5mu"

echo "> Spatial neighbor graph calculation complete"

aws s3 sync \
    --profile di \
    "$DIR" \
    s3://openpipelines-bio/openpipeline_spatial/resources_test/aviti \
    --delete \
    --dryrun