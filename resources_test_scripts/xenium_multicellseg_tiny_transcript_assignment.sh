#!/bin/bash

# Derive a Baysor-style transcript assignment (segmentation.csv) and viz polygons
# (segmentation_polygons_2d.json) from the XOA segmentation in xenium_multicellseg_tiny_raw,
# for testing xeniumranger import-segmentation --transcript-assignment / --viz-polygons.
#
# Requires a python3 with pandas and pyarrow (override with PYTHON=/path/to/python).

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# Define absolute directory paths
DIR="$REPO_ROOT/resources_test/xenium"
SRC="$DIR/xenium_multicellseg_tiny_raw"
ID="xenium_multicellseg_tiny_transcript_assignment"
OUT="$DIR/$ID"
PYTHON="${PYTHON:-python3}"

if [ ! -d "$SRC" ]; then
    echo "Source bundle '$SRC' not found." >&2
    exit 1
fi

if ! "$PYTHON" -c "import pandas, pyarrow" 2>/dev/null; then
    echo "'$PYTHON' needs pandas and pyarrow installed (or set PYTHON=/path/to/python)." >&2
    exit 1
fi

mkdir -p "$OUT"

SRC="$SRC" OUT="$OUT" "$PYTHON" - <<'EOF'
import json
import os
from pathlib import Path

import pandas as pd

src = Path(os.environ["SRC"])
out = Path(os.environ["OUT"])

tx = pd.read_parquet(src / "transcripts.parquet")
bounds = pd.read_parquet(src / "cell_boundaries.parquet")

# Baysor-style cell ids ("cell-N"), as expected by xeniumranger's cell id parser
assigned = tx["cell_id"] != "UNASSIGNED"
cells = sorted(tx.loc[assigned, "cell_id"].unique())
assert set(cells) == set(bounds["cell_id"]), "every assigned cell needs a polygon"
hq_cells = set(tx.loc[assigned & (tx["qv"] >= 20), "cell_id"])
assert hq_cells == set(cells), "every cell needs >= 1 high-quality (qv >= 20) transcript"
cell_map = {c: f"cell-{i + 1}" for i, c in enumerate(cells)}

# Transcript assignment CSV (Baysor segmentation.csv layout).
# x/y stay in microns: xeniumranger requires --units microns with --transcript-assignment.
seg = pd.DataFrame(
    {
        "transcript_id": tx["transcript_id"],
        "x": tx["x_location"].round(4),
        "y": tx["y_location"].round(4),
        "z": tx["z_location"],
        "gene": tx["feature_name"],
        "qv": tx["qv"],
        "cell": tx["cell_id"].map(cell_map).fillna(""),
        "is_noise": (~assigned).map({True: "true", False: "false"}),
    }
)
seg.to_csv(out / "segmentation.csv", index=False)

# Viz polygons: GeoJSON FeatureCollection, top-level "id" matching the CSV "cell"
features = []
for cell_id, grp in bounds.groupby("cell_id", sort=True):
    ring = grp[["vertex_x", "vertex_y"]].to_numpy(dtype=float).round(4).tolist()
    if ring[0] != ring[-1]:
        ring.append(ring[0])
    assert len(ring) >= 4, f"polygon for {cell_id} has fewer than 4 vertices"
    features.append(
        {
            "type": "Feature",
            "id": cell_map[cell_id],
            "geometry": {"type": "Polygon", "coordinates": [ring]},
            "properties": {"cell_id": cell_map[cell_id]},
        }
    )
with open(out / "segmentation_polygons_2d.json", "w") as f:
    json.dump({"type": "FeatureCollection", "features": features}, f)

print(f"{len(seg)} transcripts ({(~assigned).sum()} noise), {len(features)} cells -> {out}")
EOF
