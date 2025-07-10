import sys
import os
import squidpy as sq
import mudata as mu

## VIASH START
par = {
    "input": "./resources_test/cosmx/Lung5_Rep2_tiny",
    "output": "./resources_test/cosmx/Lung5_Rep2_tiny.h5mu",
    "dataset_id": "Lung5_Rep2",
    "modality": "rna",
    "output_compression": None,
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

counts_file = f"{par['dataset_id']}_exprMat_file.csv"
fov_file = f"{par['dataset_id']}_fov_positions_file.csv"
meta_file = f"{par['dataset_id']}_metadata_file.csv"

for file in [counts_file, fov_file, meta_file]:
    assert os.path.isfile(os.path.join(par["input"], file)), (
        f"File does not exist: {file}"
    )

logger.info("Reading in CosMx data...")
adata = sq.read.nanostring(
    path=par["input"], counts_file=counts_file, meta_file=meta_file, fov_file=fov_file
)

logger.info("Writing output MuData object...")
mdata = mu.MuData({par["modality"]: adata})
mdata.write_h5mu(par["output"], compression=par["output_compression"])
