import sys
import os
import squidpy as sq
import mudata as mu
import zipfile_inflate64 as zipfile
from pathlib import Path

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
from unzip_archived_folder import extract_selected_files_from_zip

logger = setup_logger()


def retrieve_input_data(cosmx_output_bundle):
    # Expected folder structure (showing only relevant files):
    # ├── *_exprMat_file.csv
    # ├── *_fov_positions_file.csv
    # └── *_metadata_file.csv

    required_file_patterns = {
        "counts_file": "**/*exprMat_file.csv",
        "fov_file": "**/*fov_positions_file.csv",
        "meta_file": "**/*metadata_file.csv",
    }
    if zipfile.is_zipfile(cosmx_output_bundle):
        cosmx_output_bundle = extract_selected_files_from_zip(
            cosmx_output_bundle, members=required_file_patterns.values()
        )
    else:
        cosmx_output_bundle = Path(cosmx_output_bundle)

    assert os.path.isdir(cosmx_output_bundle), (
        "Input is expected to be a (compressed) directory."
    )

    input_data = {}
    for key, pattern in required_file_patterns.items():
        file = list(cosmx_output_bundle.glob(pattern))
        assert len(file) == 1, f"Expected one file for {key}, found {len(file)}."
        input_data[key] = file[0]

    parent_dirs = {file.parent for file in input_data.values()}
    assert len(parent_dirs) == 1, (
        f"Input files are expected to be in the same directory."
        f"Found files in {len(parent_dirs)} different directories: {parent_dirs}"
    )

    return input_data


def main():
    logger.info("Reading in CosMx data...")
    input_files = retrieve_input_data(par["input"])

    adata = sq.read.nanostring(
        path=input_files["counts_file"].parent,
        counts_file=input_files["counts_file"].name,
        meta_file=input_files["meta_file"].name,
        fov_file=input_files["fov_file"].name,
    )

    logger.info("Writing output MuData object...")
    mdata = mu.MuData({par["modality"]: adata})
    mdata.write_h5mu(par["output"], compression=par["output_compression"])


if __name__ == "__main__":
    main()
