import sys
from pathlib import Path
import scanpy as sc
import pandas as pd
import mudata as mu
import zipfile
import json
import os

## VIASH START
par = {
    "input": "test/xenium_tiny.zip",
    "output": "xenium_tiny_test.h5mu",
    "output_compression": "gzip",
    "obsm_coordinates": "spatial",
    "uns_experiment": "xenium_experiment",
    "uns_metrics": "xenium_metrics",
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from unzip_archived_folder import extract_selected_files_from_zip

logger = setup_logger()


def _retrieve_input_data(xenium_output_bundle):
    # Expected folder structure (showing only relevant files):
    # ├── cell_feature_matrix.h5
    # ├── cells.parquet
    # ├── experiment.xenium
    # └── metrics_summary.csv

    required_files = [
        "cell_feature_matrix.h5",
        "cells.parquet",
        "experiment.xenium",
        "metrics_summary.csv",
    ]

    if zipfile.is_zipfile(xenium_output_bundle):
        xenium_output_bundle = extract_selected_files_from_zip(
            xenium_output_bundle, members=required_files
        )

    assert os.path.isdir(xenium_output_bundle), (
        "Input is expected to be a (compressed) directory."
    )
    input_dir = Path(xenium_output_bundle)

    input_data = dict(
        zip(
            ["count_matrix", "cells_metadata", "experiment", "metrics_summary"],
            [input_dir / file for file in required_files],
        )
    )

    assert all([file.exists() for file in input_data.values()]), (
        f"Not all required input files are found. Make sure that {par['input']} contains {input_data.values()}."
    )
    return input_data


def _format_cell_id_column(cell_id_column: pd.Series) -> pd.Series:
    """Convert cell IDs to string format, decoding bytes if necessary."""
    return cell_id_column.apply(
        lambda x: x.decode("utf-8") if isinstance(x, bytes) else str(x)
    )


# Read data from Xenium output bundle
logger.info("Reading input data...")

input_data = _retrieve_input_data(par["input"])

adata = sc.read_10x_h5(input_data["count_matrix"])
metadata = pd.read_parquet(input_data["cells_metadata"], engine="pyarrow")
with open(input_data["experiment"], "r") as f:
    specs = json.load(f)
metrics_summary = pd.read_csv(
    input_data["metrics_summary"], decimal=".", quotechar='"', thousands=","
)

# Extract and format required columns
cell_ids = _format_cell_id_column(metadata["cell_id"])
coordinates = metadata[["x_centroid", "y_centroid"]].to_numpy()
metadata.drop(["cell_id", "x_centroid", "y_centroid"], axis=1, inplace=True)

# Updata AnnData with metadata
adata.obs = metadata
adata.obs_names = cell_ids
adata.obsm[par["obsm_coordinates"]] = coordinates
adata.uns[par["uns_experiment"]] = specs
adata.uns[par["uns_metrics"]] = metrics_summary

# Write output MuData
mdata = mu.MuData({"rna": adata})
mdata.write_h5mu(par["output"], compression=par["output_compression"])
