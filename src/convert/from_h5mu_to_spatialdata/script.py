import logging
import sys

import mudata as mu
import spatialdata as sd

## VIASH START
par = {
    "input": "./resources_test/xenium/xenium_tiny.h5mu",
    "input_spatialdata": "./resources_test/xenium/xenium_tiny.zarr",
    "input_table": "table",
    "output": "./resources_test/xenium/xenium_tiny_from_h5mu.zarr",
    "output_table": "table",
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

if par.get("input_spatialdata", None) is not None:
    logger.info(f"Reading existing SpatialData from {par['input_spatialdata']}...")

    # Disable logger messages from spatialdata when reading
    logger.setLevel(logging.WARNING)
    sdata_existing = sd.read_zarr(par["input_spatialdata"])
    logger.setLevel(logging.INFO)

    input_table = par["input_table"]
    if input_table not in sdata_existing.tables:
        raise ValueError(
            f"Table '{input_table}' not found in existing SpatialData. "
            f"Available tables: {list(sdata_existing.tables.keys())}"
        )

    logger.info(
        f"Checking modality matches existing SpatialData table '{input_table}'..."
    )
    if not mod.n_obs == sdata_existing[input_table].n_obs:
        raise ValueError(
            "The number of observations in the selected modality does not match the existing SpatialData table."
        )
    if not mod.obs_names.equals(sdata_existing[input_table].obs_names):
        raise ValueError(
            "The observation names in the selected modality do not match the existing SpatialData table."
        )

logger.info("Creating SpatialData object...")
if par.get("input_spatialdata", None) is not None:
    logger.info("Using existing SpatialData...")
    sdata = sdata_existing

    # Make sure mod is a compatible SpatialData table
    attrs = mod.uns["spatialdata_attrs"]
    mod = sd.models.TableModel.parse(
        mod,
        region=attrs["region"].tolist()
        if hasattr(attrs["region"], "tolist")
        else attrs["region"],
        region_key=attrs["region_key"],
        instance_key=attrs["instance_key"],
        overwrite_metadata=True,
    )

    # The input table is kept unless the modality is stored under the same name
    output_table = par["output_table"]
    if output_table != input_table and output_table in sdata.tables:
        logger.warning(
            f"Output table '{output_table}' already exists and will be overwritten"
        )
    sdata[output_table] = mod
else:
    logger.info("Creating new SpatialData...")
    sdata = sd.SpatialData(tables={par["output_table"]: mod})

logger.info(f"Writing output SpatialData Zarr store to {par['output']}...")
sdata.write(par["output"], overwrite=True)

logger.info("Done!")
