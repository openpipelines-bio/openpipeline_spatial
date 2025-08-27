import sys
import os
import squidpy as sq
import mudata as mu
import glob
import zipfile

## VIASH START
par = {
    "input": "./resources_test/cosmx/Lung5_Rep2_tiny.zip",
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


def find_cosmx_files(cosmx_output_bundle, suffix):
    pattern = os.path.join(cosmx_output_bundle, f"*{suffix}")
    files = glob.glob(pattern)
    assert len(files) == 1, (
        f"Only one file matching pattern {pattern} should be present"
    )
    return files[0]


def retrieve_input_data(cosmx_output_bundle):
    # Expected folder structure (showing only relevant files):
    # ├── *_exprMat_file.csv
    # ├── *_fov_positions_file.csv
    # └── *_metadata_file.csv

    expected_file_patterns = [
        "exprMat_file.csv",
        "fov_positions_file.csv",
        "metadata_file.csv",
    ]
    if zipfile.is_zipfile(cosmx_output_bundle):
        cosmx_output_bundle = extract_selected_files_from_zip(
            cosmx_output_bundle, members=["*" + file for file in expected_file_patterns]
        )

    assert os.path.isdir(cosmx_output_bundle), (
        "Input is expected to be a (compressed) directory."
    )

    input_data = dict(
        zip(
            ["counts_file", "fov_file", "meta_file"],
            [
                find_cosmx_files(cosmx_output_bundle, glob_pattern)
                for glob_pattern in expected_file_patterns
            ],
        )
    )

    return input_data


def main():
    logger.info("Reading in CosMx data...")
    input_data = retrieve_input_data(par["input"])

    adata = sq.read.nanostring(
        path=par["input"],
        counts_file=input_data["counts_file"],
        meta_file=input_data["meta_file"],
        fov_file=input_data["fov_file"],
    )

    logger.info("Writing output MuData object...")
    mdata = mu.MuData({par["modality"]: adata})
    mdata.write_h5mu(par["output"], compression=par["output_compression"])


if __name__ == "__main__":
    main()
