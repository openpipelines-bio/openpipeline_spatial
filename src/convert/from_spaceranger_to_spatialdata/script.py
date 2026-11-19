import sys
from pathlib import Path

import spatialdata_io

## VIASH START
par = {
    "input": "resources_test/visium/Visium_FFPE_Human_Ovarian_Cancer_tiny_spaceranger",
    "output": "output.zarr",
    "output_table": "table",
    "output_type": "filtered",
    "output_layer": None,
    "dataset_id": None,
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger  # noqa: E402

logger = setup_logger()


def main(par):
    input_path = Path(par["input"])
    dataset_id = par["dataset_id"] or input_path.name
    counts_file = f"{par['output_type']}_feature_bc_matrix.h5"

    logger.info("Reading Space Ranger output from '%s'", input_path)
    logger.info("spatialdata-io version: %s", spatialdata_io.__version__)
    logger.info("dataset_id: %s", dataset_id)
    logger.info("counts_file: %s", counts_file)

    sdata = spatialdata_io.visium(
        path=input_path,
        counts_file=counts_file,
        dataset_id=dataset_id,
    )

    # visium() names the spot table "table"
    if par["output_table"] != "table":
        logger.info("Renaming table 'table' to '%s'", par["output_table"])
        sdata[par["output_table"]] = sdata["table"]
        del sdata["table"]

    logger.info("SpatialData object:\n%s", sdata)

    if par["output_layer"]:
        logger.info("Storing counts as layer '%s'", par["output_layer"])
        sdata[par["output_table"]].layers[par["output_layer"]] = sdata[
            par["output_table"]
        ].X.copy()

    logger.info("Writing SpatialData to '%s'", par["output"])
    sdata.write(par["output"], overwrite=True)

    logger.info("Done")


if __name__ == "__main__":
    main(par)
