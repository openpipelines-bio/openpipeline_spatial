import sys
import spatialdata as sd

################################################################################
# VIASH
################################################################################

## VIASH START
par = {
    "input": "input.zarr",
    "input_layer": None,
    "output": "output.zarr",
    "output_layer": None,
    "delete_input_layer": True,
}
## VIASH END

################################################################################
# FUNCTIONS
################################################################################


def move_layer(sdata, input_layer, output_layer, delete_input):
    """Move a layer within SpatialData object."""
    if input_layer:
        layer_data = sdata["table"].layers[input_layer]
        if delete_input:
            print(f"Deleting input layer '{input_layer}'...", flush=True)
            del sdata["table"].layers[input_layer]
    else:
        layer_data = sdata["table"].X
        if delete_input:
            print("Deleting input X matrix...", flush=True)
            sdata["table"].X = None

    if output_layer:
        sdata["table"].layers[output_layer] = layer_data
    else:
        sdata["table"].X = layer_data


################################################################################
# MAIN
################################################################################


def main(par):
    print(f"====== Move layer (spatialdata v{sd.__version__}) ======", flush=True)

    print(f"\n>>> Loading SpatialData from '{par['input']}'...", flush=True)
    sdata = sd.read_zarr(par["input"])
    print(sdata, flush=True)
    print(sdata["table"], flush=True)

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
        print(
            f"WARNING: Output layer '{output_layer_name}' already exists and will be overwritten",
            flush=True,
        )

    print(
        f"\n>>> Moving layer '{input_layer_name}' to '{output_layer_name}'...",
        flush=True,
    )
    move_layer(
        sdata, par["input_layer"], par["output_layer"], par["delete_input_layer"]
    )
    print(sdata, flush=True)
    print(sdata["table"], flush=True)

    print(f"\n>>> Writing output to '{par['output']}'...", flush=True)
    sdata.write(par["output"])
    print("Done!", flush=True)


if __name__ == "__main__":
    main(par)
