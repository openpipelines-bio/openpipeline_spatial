import sys
import squidpy as sq
import mudata as mu

## VIASH START
par = {
    # Inputs
    "input": "resources_test/cosmx/Lung5_Rep2_tiny.h5mu",
    "modality": "rna",
    "layer": None,
    "input_gp_mask": "resources_test/niche/prior_knowledge_gp_mask.json",
    "input_obsm_spatial_coords": "spatial",
    ## Spatial neighbor calculation
    "n_spatial_neighbors": 4,
    "coord_type": "generic",
    "delaunay": False,
    "output": "foo.h5mu"
}

meta = {"resources_dir": "src/utils/"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

## Read in data
adata = mu.read_h5ad(par["input"], mod=par["modality"])

## Compute spatial neighbor graph
logger.info("Computing spatial neighbor graph...")
sq.gr.spatial_neighbors(
    adata,
    coord_type=par["coord_type"],
    spatial_key=par["input_obsm_spatial_coords"],
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
mdata.write_h5mu(par["output"])
