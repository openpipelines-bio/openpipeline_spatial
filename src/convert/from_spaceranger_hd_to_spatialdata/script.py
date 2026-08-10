import sys
from pathlib import Path

import spatialdata_io

## VIASH START
par = {
    "input": "resources_test/visium_hd/Visium_HD_Mouse_Brain_tiny_spaceranger",
    "output": "output.zarr",
    "mode": "segmented_cells",
    "bin_size": 8,
    "output_type": "filtered",
    "output_layer": None,
    "dataset_id": None,
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger  # noqa: E402

logger = setup_logger()

# Space Ranger names the segmented-cell table "cell_segmentations"
# (spatialdata_io VisiumHDKeys.CELL_SEG_KEY_HD).
CELL_TABLE_KEY = "cell_segmentations"


def main(par):
    input_path = Path(par["input"])
    filtered = par["output_type"] == "filtered"

    logger.info("Reading Space Ranger HD output from '%s'", input_path)
    logger.info("spatialdata-io version: %s", spatialdata_io.__version__)
    logger.info("mode: %s", par["mode"])
    logger.info("dataset_id: %s", par["dataset_id"] or "(inferred from feature slice)")
    logger.info("output_type: %s", par["output_type"])

    if par["mode"] == "bins":
        logger.info("bin_size: %s um", par["bin_size"])
        sdata = spatialdata_io.visium_hd(
            path=input_path,
            dataset_id=par["dataset_id"],
            bin_size=[par["bin_size"]],
            filtered_counts_file=filtered,
            load_segmentations_only=False,
        )
        # visium_hd() names the binned table after the bin size, e.g. "square_008um".
        table_key = f"square_{par['bin_size']:03d}um"
    else:  # segmented_cells
        sdata = spatialdata_io.visium_hd(
            path=input_path,
            dataset_id=par["dataset_id"],
            filtered_counts_file=filtered,
            load_segmentations_only=True,
            load_nucleus_segmentations=False,
        )
        table_key = CELL_TABLE_KEY

    # Rename the loaded table to "table" so downstream components stay
    # technology-agnostic (the shapes and images keep their dataset-id names).
    sdata["table"] = sdata[table_key]
    del sdata[table_key]

    logger.info("SpatialData object:\n%s", sdata)

    if par["output_layer"]:
        logger.info("Storing counts as layer '%s'", par["output_layer"])
        sdata["table"].layers[par["output_layer"]] = sdata["table"].X.copy()

    logger.info("Writing SpatialData to '%s'", par["output"])
    sdata.write(par["output"], overwrite=True)

    logger.info("Done")


if __name__ == "__main__":
    main(par)
