import sys
import os
import squidpy as sq
import mudata as mu
import glob

## VIASH START
par = {
    "input": "./resources_test/cosmx/Lung5_Rep2_tiny",
    "output": "./resources_test/cosmx/Lung5_Rep2_tiny.h5mu",
    "modality": "rna",
    "output_compression": None,
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

def find_matrix_file(suffix):
    pattern = os.path.join(par["input"], f"*{suffix}")
    files = glob.glob(pattern)
    assert len(files) == 1, f"Only one file matching pattern {pattern} should be present"
    return files[0]

counts_file = find_matrix_file("exprMat_file.csv")
fov_file = find_matrix_file("fov_positions_file.csv")
meta_file = find_matrix_file("metadata_file.csv")

logger.info("Reading in CosMx data...")
adata = sq.read.nanostring(
    path=par["input"], counts_file=counts_file, meta_file=meta_file, fov_file=fov_file
)

logger.info("Writing output MuData object...")
mdata = mu.MuData({par["modality"]: adata})
mdata.write_h5mu(par["output"], compression=par["output_compression"])
