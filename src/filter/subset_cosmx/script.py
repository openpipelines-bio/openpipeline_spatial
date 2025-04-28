import os
import shutil
import pandas as pd
import glob
import sys


## VIASH START
par = {
    "input": "./resources_test/cosmx/Lung5_Rep2",
    "output": "./resources_test/cosmx/Lung5_Rep2_tiny/",
    "dataset_id": "Lung5_Rep2",
    "num_fovs": 5
}
meta ={
    "resources_dir": "src/utils"
}
## VIASH END


sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

counts_file = f"{par['dataset_id']}_exprMat_file.csv"
fov_file = f"{par['dataset_id']}_fov_positions_file.csv"
meta_file = f"{par['dataset_id']}_metadata_file.csv"
tx_file = f"{par['dataset_id']}_tx_file.csv"

for file in [counts_file, fov_file, meta_file]:
    assert os.path.isfile(os.path.join(par["input"], file)), f"File does not exist: {file}"

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
matrices = [counts_file, fov_file, meta_file, tx_file]
for matrix in matrices:
    logger.info(f"Subsetting {matrix}, keeping fovs {kept_fovs}")
    data = pd.read_csv(os.path.join(par["input"], matrix))
    data_tiny = data[data["fov"].isin(kept_fovs)]
    data_tiny.to_csv(os.path.join(par["output"], matrix), index=False)
