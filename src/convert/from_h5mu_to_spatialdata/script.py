import sys
import spatialdata as sd
import mudata as mu
import logging

## VIASH START
par = {
    "input": "./resources_test/xenium/xenium_tiny.h5mu",
    "input_spatialdata": "./resources_test/xenium/xenium_tiny.zarr",
    "output": "./resources_test/xenium/xenium_tiny_from_h5mu.zarr",
    "modality": "rna",
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

logger.info("Starting conversion from H5MU to SpatialData...")

logger.info(f"Reading input H5MU file from {par['input']}...")
mdata = mu.read_h5mu(par["input"])

logger.info("Extracting modality from MuData object...")
mod = mdata.mod[par["modality"]]

if (par.get("input_spatialdata", None) is not None):
    logger.info(f"Reading existing SpatialData from {par['input_spatialdata']}...")
    
    # Disable logger messages from spatialdata when reading
    logger.setLevel(logging.WARNING)
    sdata_existing = sd.read_zarr(par["input_spatialdata"])
    logger.setLevel(logging.INFO)

    logger.info("Checking modality matches existing SpatialData table...")
    if not mod.n_obs == sdata_existing["table"].n_obs:
        raise ValueError(
            "The number of observations in the selected modality does not match the existing SpatialData table."
        )
    if not mod.obs_names.equals(sdata_existing["table"].obs_names):
        raise ValueError(
            "The observation names in the selected modality do not match the existing SpatialData table."
        )

logger.info("Creating SpatialData object...")
if (par.get("input_spatialdata", None) is not None):
    logger.info("Using existing SpatialData...")
    sdata = sdata_existing
    sdata["table"] = mod
else:
    logger.info("Creating new SpatialData...")
    sdata = sd.SpatialData(tables={"table": mod})

logger.info(f"Writing output SpatialData Zarr store to {par['output']}...")
sdata.write(par["output"], overwrite=True)

logger.info("Done!")
