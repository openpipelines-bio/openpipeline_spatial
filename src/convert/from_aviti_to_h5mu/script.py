import sys
from pathlib import Path
import scanpy as sc
import pandas as pd
import mudata as mu
import json

## VIASH START
par = {
    "input": "./resources_test/aviti/aviti_teton_tiny",
    "output": "xenium_tiny_test.h5mu",
    "output_compression": "gzip",
    "obsm_coordinates": "spatial",
    "uns_gene_panel": "aviti_gene_panel",
    "uns_run_manifest": "aviti_run_manifest",
    "uns_run_parameters": "aviti_run_paraemters",
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

# Expected folder structure (showing only relevant files):
# ├── Cytoprofiling/
# │   └── Instrument/
# │       └── RawCellStats_subset.parquet
# ├── Panel.json
# ├── RunManifest.json
# └── RunParameters.json
input_dir = Path(par["input"])
input_data = {
    "count_matrix": input_dir / "Cytoprofiling" / "Instrument" / "RawCellStats.parquet",
    "panel_metadata": input_dir / "Panel.json",
    "manifest_metadata": input_dir / "RunManifest.json",
    "parameter_metadata": input_dir / "RunParameters.json",
}

def _format_cell_id_column(cell_id_column: pd.Series) -> pd.Series:
    """Convert cell IDs to string format, decoding bytes if necessary."""
    return cell_id_column.apply(
        lambda x: x.decode("utf-8") if isinstance(x, bytes) else str(x)
    )


# Read data from Xenium output bundle
logger.info("Reading input data...")

assert all([file.exists() for file in input_data.values()]), (
    f"Not all required input files are found. Make sure that {par['input']} contains {input_data.values()}."
)

df = pd.read_parquet(input_data["count_matrix"], engine="pyarrow")
with open(input_data["manifest_metadata"], "r") as f:
    manifest_metadata = json.load(f)
with open(input_data["parameter_metadata"], "r") as f:
    parameter_metadata = json.load(f)
with open(input_data["panel_metadata"], "r") as f:
    panel_metadata = json.load(f)

# Process obs data
obs_columns = [
    "Area", "AreaUm", "Cell", "NuclearArea", "NuclearAreaUm", 
    "Tile", "Well", "WellLabel", "X", "Xum", "Y", "Yum"
]
# Extract and format required columns
cell_ids = _format_cell_id_column(metadata["cell_id"])
coordinates = df[["x_centroid", "y_centroid"]].to_numpy()
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
