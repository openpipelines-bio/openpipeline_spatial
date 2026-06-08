import sys
import spatialdata as sd

## VIASH START
par = {
    "input": "input.zarr",
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


def move_layer(sdata, input_layer, output_layer, delete_input):
    """Move a layer within SpatialData object."""
    if input_layer:
        layer_data = sdata["table"].layers[input_layer]
        if delete_input:
            logger.info(f"Deleting input layer '{input_layer}'...")
            del sdata["table"].layers[input_layer]
    else:
        layer_data = sdata["table"].X
        if delete_input:
            logger.info("Deleting input X matrix...")
            sdata["table"].X = None

    if output_layer:
        sdata["table"].layers[output_layer] = layer_data
    else:
        sdata["table"].X = layer_data


logger.info(f"Move layer (spatialdata v{sd.__version__})")

logger.info(f"Loading SpatialData from '{par['input']}'...")
sdata = sd.read_zarr(par["input"])
logger.info(f"SpatialData: {sdata}")
logger.info(f"Table: {sdata['table']}")

input_layer_name = par["input_layer"] if par["input_layer"] else "X"
output_layer_name = par["output_layer"] if par["output_layer"] else "X"

if input_layer_name == output_layer_name:
    raise ValueError(
        f"Input layer '{input_layer_name}' and output layer '{output_layer_name}' are the same, aborting"
    )

if par["input_layer"] and input_layer_name not in sdata["table"].layers:
    raise ValueError(
        f"Input layer '{input_layer_name}' not found in SpatialData. Available layers: {list(sdata['table'].layers.keys())}"
    )

if par["output_layer"] and output_layer_name in sdata["table"].layers:
    logger.warning(
        f"Output layer '{output_layer_name}' already exists and will be overwritten"
    )

logger.info(f"Moving layer '{input_layer_name}' to '{output_layer_name}'...")
move_layer(sdata, par["input_layer"], par["output_layer"], par["delete_input_layer"])
logger.info(f"SpatialData: {sdata}")
logger.info(f"Table: {sdata['table']}")

logger.info(f"Writing output to '{par['output']}'...")
sdata.write(par["output"], overwrite=True)

logger.info("Done!")
