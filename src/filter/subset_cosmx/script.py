import os
import shutil
import pandas as pd
import glob
import sys


## VIASH START
par = {
    "input": "./resources_test/cosmx/Lung5_Rep2",
    "output": "./resources_test/cosmx/Lung5_Rep2_tiny/",
    "subset_transcripts_file": True,
    "subset_polygons_file": False,
    "num_fovs": 5,
}
meta = {"resources_dir": "src/utils"}
## VIASH END


sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()


def find_matrix_file(suffix):
    pattern = os.path.join(par["input"], f"*{suffix}")
    files = glob.glob(pattern)
    assert len(files) == 1, (
        f"Only one file matching pattern {pattern} should be present"
    )
    return files[0]


kept_fovs = list(range(1, par["num_fovs"] + 1))

os.makedirs(par["output"], exist_ok=True)

# Images
image_dirs = ["CellComposite", "CellLabels", "CellOverlay", "CompartmentLabels"]

for image_dir in image_dirs:
    logger.info(f"Subsetting {image_dir}, keeping fovs {kept_fovs}")
    os.makedirs(f"{par['output']}/{image_dir}", exist_ok=True)
    for fov in kept_fovs:
        fov_str = f"{image_dir}_F{fov:03d}.*"

        file_path = glob.glob(os.path.join(par["input"], image_dir, fov_str))
        assert len(file_path) == 1
        shutil.copy2(file_path[0], os.path.join(par["output"], image_dir))

# Matrices
counts_file = find_matrix_file("exprMat_file.csv")
fov_file = find_matrix_file("fov_positions_file.csv")
meta_file = find_matrix_file("metadata_file.csv")

matrices = [counts_file, fov_file, meta_file]
if par["subset_transcripts_file"]:
    tx_file = find_matrix_file("tx_file.csv")
    matrices.append(tx_file)
if par["subset_polygons_file"]:
    polygons_file = find_matrix_file("polygons.csv")
    matrices.append(polygons_file)

for matrix in matrices:
    logger.info(f"Subsetting {matrix}, keeping fovs {kept_fovs}")
    data = pd.read_csv(matrix)
    data_tiny = data[data["fov"].isin(kept_fovs)]
    data_tiny.to_csv(os.path.join(par["output"], os.path.basename(matrix)), index=False)
