import os
import shutil
import sys

import numpy as np
import pandas as pd
import spatialdata as sd

## VIASH START
par = {
    "input": "input.zarr",
    "input_table": "table",
    "output": "output.zarr",
    "output_obsm_spatial": "spatial",
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()


def get_shapes_region(sdata, table_name):
    """Get the shapes element and instance key that a table annotates."""
    attrs = sdata[table_name].uns.get("spatialdata_attrs", {})
    region = attrs.get("region")
    instance_key = attrs.get("instance_key")
    regions = region if isinstance(region, list) else [region]
    if len(regions) > 1:
        raise ValueError(
            f"Table '{table_name}' annotates multiple regions ({region}). "
            "This component only supports a single region."
        )
    region = regions[0] if regions else None
    if region is None or instance_key is None:
        raise ValueError(
            f"Table '{table_name}' has no spatialdata_attrs region/instance_key. "
            "Cannot determine which shapes element it annotates."
        )
    if region not in sdata.shapes:
        raise ValueError(
            f"Region '{region}' annotated by table '{table_name}' is not a shapes "
            f"element. Available shapes: {list(sdata.shapes)}"
        )

    return region, instance_key


logger.info(f"Derive .obsm spatial from shapes (spatialdata v{sd.__version__})")

# Copy the store and then load and rewrite only the table. This avoids loading
# images, labels and points into memory. If the output is the input, the table
# is updated in place.
if os.path.realpath(par["input"]) == os.path.realpath(par["output"]):
    logger.info(f"Output is the same as input, updating '{par['output']}' in place")
else:
    if os.path.exists(par["output"]):
        logger.info(f"Removing existing output '{par['output']}'...")
        shutil.rmtree(par["output"])
    logger.info(f"Copying SpatialData from '{par['input']}' to '{par['output']}'...")
    shutil.copytree(par["input"], par["output"])

logger.info(f"Loading tables and shapes from '{par['output']}'...")
sdata = sd.read_zarr(par["output"], selection=("tables", "shapes"))
logger.info(f"SpatialData: {sdata}")

if par["input_table"] not in sdata.tables:
    raise ValueError(
        f"Table '{par['input_table']}' not found in '{par['input']}'. "
        f"Available tables: {list(sdata.tables)}"
    )

region, instance_key = get_shapes_region(sdata, par["input_table"])
logger.info(f"Region: '{region}', instance key: '{instance_key}'")

table = sdata[par["input_table"]]
shapes = sdata.shapes[region]
instance_values = table.obs[instance_key].to_numpy()

missing = ~pd.Index(instance_values).isin(shapes.index)
if missing.any():
    raise ValueError(
        f"{missing.sum()} table instances not found in shapes '{region}': "
        f"{list(instance_values[missing][:5])}..."
    )

centroids = shapes.loc[instance_values].geometry.centroid
table.obsm[par["output_obsm_spatial"]] = np.column_stack([centroids.x, centroids.y])
logger.info(
    f"Wrote {table.obsm[par['output_obsm_spatial']].shape} centroids to "
    f".obsm['{par['output_obsm_spatial']}']"
)

logger.info(f"Writing table to '{par['output']}'...")
# The store can't be overwritten because it is in use so we manually delete and
# write the table
sdata.delete_element_from_disk(par["input_table"])
sdata.write_element(par["input_table"])
sdata.write_consolidated_metadata()

logger.info("Done!")
