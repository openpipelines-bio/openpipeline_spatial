import sys
from spatialdata_io import xenium

## VIASH START
par = {
    "input": "./resources_test/xenium_tiny",
    "output": "./test/xenium_tiny.zarr",
    "cells_boundaries": True,
    "nucleus_boundaries": True,
    "cells_as_circles": None,
    "cells_labels": True,
    "nucleus_labels": True,
    "transcripts": True,
    "morphology_mip": True,
    "morphology_focus": True,
    "aligned_images": True,
    "cells_table": True,
    "n_jobs": 1,
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

logger.info("Reading in Xenium data...")
sdata = xenium(
    par["input"],
    cells_boundaries=par["cells_boundaries"],
    nucleus_boundaries=par["nucleus_boundaries"],
    cells_as_circles=par["cells_as_circles"],
    cells_labels=par["cells_labels"],
    nucleus_labels=par["nucleus_labels"],
    transcripts=par["transcripts"],
    morphology_mip=par["morphology_mip"],  # only available in version < 2.0.0
    morphology_focus=par["morphology_focus"],
    aligned_images=par["aligned_images"],
    cells_table=par["cells_table"],
    n_jobs=par["n_jobs"],
)


logger.info("Writing out SpatialData object to Zarr...")
sdata.write(par["output"], overwrite=True)
