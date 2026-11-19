import os
import shutil
import sys
import spatialdata as sd

## VIASH START
par = {
    "input": "input.zarr",
    "input_table": "table",
    "input_layer": None,
    "output": "output.zarr",
    "output_layer": None,
    "delete_input_layer": True,
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()


def move_layer(table, input_layer, output_layer, delete_input):
    """Move a layer within a SpatialData table."""
    if input_layer:
        layer_data = table.layers[input_layer]
        if delete_input:
            logger.info(f"Deleting input layer '{input_layer}'...")
            del table.layers[input_layer]
    else:
        layer_data = table.X
        if delete_input:
            logger.info("Deleting input X matrix...")
            table.X = None

    if output_layer:
        table.layers[output_layer] = layer_data
    else:
        table.X = layer_data


logger.info(f"Move layer (spatialdata v{sd.__version__})")

# Copy the store and then load and rewrite only the table. This avoids loading
# the other elements into memory
if os.path.exists(par["output"]):
    logger.info(f"Removing existing output '{par['output']}'...")
    shutil.rmtree(par["output"])
logger.info(f"Copying SpatialData from '{par['input']}' to '{par['output']}'...")
shutil.copytree(par["input"], par["output"])

logger.info(f"Loading table from '{par['output']}'...")
sdata = sd.read_zarr(par["output"], selection=("tables",))
logger.info(f"SpatialData: {sdata}")

if par["input_table"] not in sdata.tables:
    raise ValueError(
        f"Table '{par['input_table']}' not found in SpatialData. Available tables: {list(sdata.tables.keys())}"
    )

table = sdata[par["input_table"]]
logger.info(f"Table '{par['input_table']}': {table}")

input_layer_name = par["input_layer"] if par["input_layer"] else "X"
output_layer_name = par["output_layer"] if par["output_layer"] else "X"

if input_layer_name == output_layer_name:
    raise ValueError(
        f"Input layer '{input_layer_name}' and output layer '{output_layer_name}' are the same, aborting"
    )

if par["input_layer"] and input_layer_name not in table.layers:
    raise ValueError(
        f"Input layer '{input_layer_name}' not found in table '{par['input_table']}'. Available layers: {[key for key in table.layers.keys() if key is not None]}"
    )

if par["output_layer"] and output_layer_name in table.layers:
    logger.warning(
        f"Output layer '{output_layer_name}' already exists and will be overwritten"
    )

logger.info(f"Moving layer '{input_layer_name}' to '{output_layer_name}'...")
move_layer(table, par["input_layer"], par["output_layer"], par["delete_input_layer"])

logger.info(f"SpatialData: {sdata}")
logger.info(f"Table '{par['input_table']}': {table}")

logger.info(f"Writing table to '{par['output']}'...")
# The store can't be overwritten because it is in use so we manually delete and
# write the table
sdata.delete_element_from_disk(par["input_table"])
sdata.write_element(par["input_table"])
sdata.write_consolidated_metadata()

logger.info("Done!")
