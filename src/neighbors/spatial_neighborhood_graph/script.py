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
    "percentile": None,
    "output": "foo.h5mu",
    "output_compression": None,
}

meta = {"resources_dir": "src/utils/"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

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

## Validate the percentile
if par["percentile"] is not None and par["coord_type"] != "generic":
    raise ValueError("--percentile is only supported when --coord_type is 'generic'.")

## Compute spatial neighbor graph
logger.info("Computing spatial neighbor graph...")
sq.gr.spatial_neighbors(
    adata,
    coord_type=par["coord_type"],
    spatial_key=par["input_obsm_spatial_coords"],
    library_key=library_key,
    n_neighs=par["n_spatial_neighbors"],
    delaunay=par["delaunay"],
    percentile=par["percentile"],
)

# Making the connectivity matrix symmetric
logger.info("Making the connectivity matrix symmetric...")
adata.obsp["spatial_connectivities"] = adata.obsp["spatial_connectivities"].maximum(
    adata.obsp["spatial_connectivities"].T
)

## Save data
# Write the modality back into a copy of the input, so that the other modalities are preserved
logger.info("Saving output data...")
write_h5ad_to_h5mu_with_compression(
    par["output"], par["input"], par["modality"], adata, par["output_compression"]
)
