import sys
import spatialdata as sd
import mudata as mu

## VIASH START
par = {
    "input": "./resources_test/xenium/xenium_tiny.zarr",
    "input_table": "table",
    "output": "./resources_test/xenium/xenium_tiny.h5mu",
    "modality": "rna",
    "output_compression": None,
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

logger.info("Reading in Xenium data...")
sdata = sd.read_zarr(par["input"])

input_table = par["input_table"]
if input_table not in sdata.tables:
    raise ValueError(
        f"Table '{input_table}' not found in SpatialData. "
        f"Available tables: {list(sdata.tables.keys())}"
    )

logger.info(f"Fetching AnnData table '{input_table}' from SpatialData object...")
adata = sdata.tables[input_table]

logger.info("Writing output MuData object...")
mdata = mu.MuData({par["modality"]: adata})
mdata.write_h5mu(par["output"], compression=par["output_compression"])
