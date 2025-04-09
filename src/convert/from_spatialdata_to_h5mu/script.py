import sys
import spatialdata as sd
import mudata as mu

## VIASH START
par = {
    "input": "./resources_test/xenium/xenium_tiny.zarr",
    "output": "./resources_test/xenium/xenium_tiny.h5mu",
    "modality": "rna",
    "output_compression": None
}
meta ={
    "resources_dir": "src/utils"
}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

logger.info("Reading in Xenium data...")
sdata = sd.read_zarr(par["input"])

logger.info("Fetching AnnData table from SpatialData object...")
adata = sdata.tables["table"]

logger.info("Writing output MuData object...")
mdata = mu.MuData({par["modality"]: adata})
mdata.write_h5mu(par["output"], compression=par["output_compression"])
