#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

ID=aviti
OUT=resources_test/$ID/
TINY_OUT=$OUT/aviti_teton_tiny/

# Create directories
[ -d "$OUT" ] || mkdir -p "$OUT"
[ -d "$TINY_OUT" ] || mkdir -p "$TINY_OUT"

echo "> Downloading Aviti Teton data"
# TODO: Replace with actual download URL when available
# wget "https://example.com/aviti_teton_data.tar.gz?download=1" -O "${OUT}/PLUT-0105.tar.gz"
# tar -xzf "${OUT}/PLUT-0105.tar.gz" -C "$OUT" --strip-components=1
# rm "${OUT}/PLUT-0105.tar.gz"

echo "> Processing and subsetting Aviti Teton data"
python <<HEREDOC
import os
import shutil
import pandas as pd
import glob
import json

src_dir = "${OUT}/PLUT-0105"
dest_dir = "${TINY_OUT}/"
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
    dest_file = os.path.join(dest_path, "RawCellStats_subset.parquet")
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

# Subset & Copy Run Manifest Metadata
manifest_src_path = os.path.join(src_dir, "RunManifest.json")
if os.path.exists(manifest_src_path):
    with open(manifest_src_path, 'r') as f:
        manifest_data = json.load(f)

    filtered_wells = [
        well for well in manifest_data["Wells"] 
        if well["WellLocation"] in wells_to_keep
    ]

    filtered_manifest = {
        "Wells": filtered_wells,
        "RunValues": manifest_data.get("RunValues", {})
    }
    manifest_dest_path = os.path.join(dest_dir, "RunManifest.json")
    with open(manifest_dest_path, 'w') as f:
        json.dump(filtered_manifest, f, indent=4)
    
    print(f"Filtered manifest: kept {len(filtered_wells)} wells out of {len(manifest_data['Wells'])}")
else:
    print(f"Warning: RunManifest.json not found at {manifest_src_path}")

# Subset & Copy Run Parameters metadata
run_params_src_path = os.path.join(src_dir, "RunParameters.json")
if os.path.exists(run_params_src_path):
    with open(run_params_src_path, 'r') as f:
        run_params_data = json.load(f)

    # Filter wells to keep only those specified in wells_to_keep
    filtered_wells = [
        well for well in run_params_data["Wells"] 
        if well["WellLocation"] in wells_to_keep
    ]

    # Extract tile names from filtered wells
    tiles_to_keep = []
    for well in filtered_wells:
        tiles_to_keep.extend([tile["Name"] for tile in well["Tiles"]])

    # Filter Tiles and RawTiles arrays
    filtered_tiles = [tile for tile in run_params_data["Tiles"] if tile in tiles_to_keep]
    filtered_raw_tiles = [tile for tile in run_params_data["RawTiles"] if tile in tiles_to_keep]

    # Create new filtered run parameters
    filtered_run_params = run_params_data.copy()
    filtered_run_params["Wells"] = filtered_wells
    filtered_run_params["Tiles"] = filtered_tiles
    filtered_run_params["RawTiles"] = filtered_raw_tiles

    # Save filtered run parameters
    run_params_dest_path = os.path.join(dest_dir, "RunParameters.json")
    with open(run_params_dest_path, 'w') as f:
        json.dump(filtered_run_params, f, indent=2)
    
    print(f"Filtered RunParameters: kept {len(filtered_wells)} wells, {len(filtered_tiles)} tiles")
else:
    print(f"Warning: RunParameters.json not found at {run_params_src_path}")

print("Processing complete!")
HEREDOC

# echo "> Removing original aviti_teton folder"
# rm -rf "$OUT"

echo "> Aviti Teton tiny dataset created successfully at $TINY_OUT"
