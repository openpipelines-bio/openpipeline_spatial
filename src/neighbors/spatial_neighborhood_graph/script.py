import sys
import squidpy as sq
import mudata as mu

## VIASH START
par = {
    # Inputs
    "input": "resources_test/cosmx/Lung5_Rep2_tiny.h5mu",
    "modality": "rna",
    "input_obsm_spatial_coords": "spatial",
    "input_obs_library_key": None,
    ## Spatial neighbor calculation
    "n_spatial_neighbors": 4,
    "coord_type": "generic",
    "delaunay": False,
    "output": "foo.h5mu",
}

meta = {"resources_dir": "src/utils/"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

## Read in data
adata = mu.read_h5ad(par["input"], mod=par["modality"])

## Validate the library key
library_key = par["input_obs_library_key"]
if library_key:
    if library_key not in adata.obs.columns:
        raise ValueError(
            f"--input_obs_library_key '{library_key}' not found in .obs of modality '{par['modality']}'."
        )
    # squidpy requires the library key to be a categorical column
    if adata.obs[library_key].dtype.name != "category":
        logger.info(f"Converting .obs['{library_key}'] to categorical...")
        adata.obs[library_key] = adata.obs[library_key].astype("category")

## Compute spatial neighbor graph
logger.info("Computing spatial neighbor graph...")
sq.gr.spatial_neighbors(
    adata,
    coord_type=par["coord_type"],
    spatial_key=par["input_obsm_spatial_coords"],
    library_key=library_key,
    n_neighs=par["n_spatial_neighbors"],
    delaunay=par["delaunay"],
)

# Making the connectivity matrix symmetric
logger.info("Making the connectivity matrix symmetric...")
adata.obsp["spatial_connectivities"] = adata.obsp["spatial_connectivities"].maximum(
    adata.obsp["spatial_connectivities"].T
)

## Save model and data
logger.info("Saving output data...")
mdata = mu.MuData({par["modality"]: adata})
mdata.write_h5mu(par["output"], compression=par["output_compression"])
