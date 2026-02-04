import sys
import spatialdata as sd
import mudata as mu

## VIASH START
par = {
    "input": "./resources_test/xenium/xenium_tiny.h5mu",
    "output": "./resources_test/xenium/xenium_tiny.zarr",
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
if ("spatialdata_attrs" in mod.uns):
    logger.info("Removing existing SpatialData attributes from modality...")
    del mod.uns["spatialdata_attrs"]

logger.info("Creating SpatialData object...")
sdata = sd.SpatialData(tables={"table": mod})

logger.info(f"Writing output SpatialData Zarr store to {par['output']}...")
sdata.write(par["output"], overwrite=True)

logger.info("Done!")
