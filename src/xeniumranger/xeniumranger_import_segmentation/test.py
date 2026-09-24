import json
import os
import sys
from pathlib import Path
import filecmp
import subprocess
import pytest
import re
import shutil


import pandas as pd
import scanpy as sc

## VIASH START
meta = {"name": "xeniumranger_import_segmentation", "resources_dir": "resources_test/xenium"}
## VIASH END

# 1. Nuclear-expansion fixture (full raw XOA bundle), ne = nuclear expansion
input_ne = meta["resources_dir"] + "/xenium_tiny/"
id_ne = "xenium_tiny_import_segmentation"

# 2. Multichannel fixture (full raw XOA bundle), mm = multichannel
input_mm = meta["resources_dir"] + "/xenium_multicellseg_tiny_raw/"
id_mm = "xenium_multicellseg_tiny_import_segmentation"
# Baysor-style transcript assignment + viz polygons derived from the mm bundle (micron units)
input_ta = meta["resources_dir"] + "/xenium_multicellseg_tiny_transcript_assignment/"

# 3. General arguments
units = "pixels"
expansion_distance = "5"

changed_files_ne = [
    "cells.csv.gz",
    "cells.parquet",
    "cells.zarr.zip",
    "cell_boundaries.csv.gz",
    "cell_boundaries.parquet",
    "nucleus_boundaries.csv.gz",
    "nucleus_boundaries.parquet",
    "cell_feature_matrix.h5",
    "cell_feature_matrix.zarr.zip",
    "transcripts.parquet",
    "transcripts.zarr.zip",
    "analysis.zarr.zip",
    "metrics_summary.csv",
    "analysis_summary.html",
    "experiment.xenium",
]
input_not_in_output_ne = [
    "aux_outputs.tar.gz",
    "analysis.tar.gz",
    "cell_feature_matrix.tar.gz",
]
output_not_in_input_ne = [
    "cell_id_map.csv.gz",
]

# xenium_multicellseg_tiny_raw is a full raw XOA bundle, same shape as xenium_tiny
changed_files_mm = [
    "cells.csv.gz",
    "cells.parquet",
    "cells.zarr.zip",
    "cell_boundaries.csv.gz",
    "cell_boundaries.parquet",
    "nucleus_boundaries.csv.gz",
    "nucleus_boundaries.parquet",
    "cell_feature_matrix.h5",
    "cell_feature_matrix.zarr.zip",
    "transcripts.parquet",
    "transcripts.zarr.zip",
    "analysis.zarr.zip",
    "metrics_summary.csv",
    "analysis_summary.html",
    "experiment.xenium",
]
input_not_in_output_mm = [
    "aux_outputs.tar.gz",
    "analysis.tar.gz",
    "cell_feature_matrix.tar.gz",
]
output_not_in_input_mm = [
    "cell_id_map.csv.gz",
]


def assert_outputs_exists(input, output, input_not_in_output, output_not_in_input):
    assert Path(output).is_dir() and len(list(Path(output).iterdir())) > 0, (
        "Output directory exists and is non-empty"
    )

    expected_output_dirs = ["analysis", "cell_feature_matrix", "morphology_focus"]

    input_files = sorted(
        f
        for f in os.listdir(input)
        if os.path.isfile(os.path.join(input, f)) and f not in input_not_in_output
    )
    output_files = sorted(
        f
        for f in os.listdir(output)
        if os.path.isfile(os.path.join(output, f)) and f not in output_not_in_input
    )

    assert input_files == output_files, (
        "Input and output directories should contain the same (type of) files"
    )

    output_dirs = sorted(
        d for d in os.listdir(output) if os.path.isdir(os.path.join(output, d))
    )

    assert output_dirs == expected_output_dirs, (
        "Expected output directory names should should match output directory names"
    )

    for d in output_dirs:
        current_dir = Path(output) / d
        assert current_dir.is_dir() and len(list(current_dir.iterdir())) > 0, (
            f"{current_dir} should exist and is non-empty"
        )

    assert all((Path(output) / f).stat().st_size > 0 for f in output_files), (
        "All output files should be non-empty"
    )


def assert_valid_files(output):
    transcripts = pd.read_parquet(Path(output) / "transcripts.parquet")
    assert len(transcripts) > 0, "transcripts.parquet should contain at least one row"

    cells = sc.read_10x_h5(Path(output) / "cell_feature_matrix.h5")
    assert cells.n_obs > 0, "cell_feature_matrix.h5 should contain at least one cell"
    assert cells.n_vars > 0, "cell_feature_matrix.h5 should contain at least one gene"

    return transcripts, cells


def assert_import_segmentation_used(output):
    with open(Path(output) / "experiment.xenium") as f:
        exp = json.load(f)

    xenium_ranger = exp.get("xenium_ranger")
    assert xenium_ranger, (
        "experiment.xenium should carry a xenium_ranger provenance block after import-segmentation"
    )
    assert "import-segmentation" in xenium_ranger.get("command_line", ""), (
        "xenium_ranger.command_line should record that import-segmentation was run"
    )
    assert xenium_ranger.get("original_analysis_uuid"), (
        "xenium_ranger.original_analysis_uuid should reference the pre-import analysis"
    )


def assert_identical(input, output, skip_files, input_not_in_output, output_not_in_input):
    input_files = sorted(
        f
        for f in os.listdir(input)
        if os.path.isfile(os.path.join(input, f)) and f not in input_not_in_output
    )

    output_files = sorted(
        f
        for f in os.listdir(output)
        if os.path.isfile(os.path.join(output, f)) and f not in output_not_in_input
    )

    assert input_files == output_files, (
        "Input and output directories should contain the same files"
    )

    for f in input_files:
        if f not in skip_files:
            assert filecmp.cmp(
                os.path.join(input, f), os.path.join(output, f), shallow=False
            ), f"{f} should be identical between input and output"


# Nuclear expansion path 
def test_basic_execution(run_component, random_path):
    nuclei = input_ne + "cells.zarr.zip"
    output = random_path()
    run_component(
        [
            "--xenium_bundle",
            input_ne,
            "--id",
            id_ne,
            "--nuclei",
            nuclei,
            "--units",
            units,
            "--expansion_distance",
            expansion_distance,
            "--output",
            output,
        ]
    )

    assert_outputs_exists(input_ne, output, input_not_in_output_ne, output_not_in_input_ne)
    assert_valid_files(output)
    assert_identical(
        input_ne, output, changed_files_ne, input_not_in_output_ne, output_not_in_input_ne
    )
    assert_import_segmentation_used(output)


def test_valid_id(run_component, random_path):
    output = random_path()
    malformed_id = ", ,"
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--xenium_bundle",
                input_ne,
                "--id",
                malformed_id,
                "--output",
                output,
            ]
        )
    assert re.search(
        r"invalid value.*--id",
        err.value.stdout.decode("utf-8"),
        re.IGNORECASE,
    )


def test_missing_file(run_component, random_path, tmp_path):
    output = random_path()
    incomplete_bundle = tmp_path / "incomplete_bundle"
    shutil.copytree(input_ne, incomplete_bundle)
    os.remove(incomplete_bundle / "experiment.xenium")
    with pytest.raises(subprocess.CalledProcessError):
        run_component(
            [
                "--xenium_bundle",
                str(incomplete_bundle),
                "--id",
                id_ne,
                "--output",
                output,
            ]
        )


# Multimodal path 
def test_basic_execution_multichannel(run_component, random_path):
    nuclei = input_mm + "cells.zarr.zip"
    output = random_path()
    run_component(
        [
            "--xenium_bundle",
            input_mm,
            "--id",
            id_mm,
            "--nuclei",
            nuclei,
            "--units",
            units,
            "--expansion_distance",
            expansion_distance,
            "--output",
            output,
        ]
    )

    assert_outputs_exists(input_mm, output, input_not_in_output_mm, output_not_in_input_mm)
    assert_valid_files(output)
    assert_identical(
        input_mm, output, changed_files_mm, input_not_in_output_mm, output_not_in_input_mm
    )
    assert_import_segmentation_used(output)

def test_cells_nuclei(run_component, random_path):
    nuclei = input_mm + "cells.zarr.zip"
    cells = input_mm + "cells.zarr.zip" 
    output = random_path()
    run_component(
        [
            "--xenium_bundle",
            input_mm,
            "--id",
            id_mm,
            "--nuclei",
            nuclei,
            "--cells", 
            cells,
            "--units",
            units,
            "--expansion_distance",
            expansion_distance,
            "--output",
            output,
        ]
    )

    assert_outputs_exists(input_mm, output, input_not_in_output_mm, output_not_in_input_mm)
    assert_valid_files(output)
    assert_identical(
        input_mm, output, changed_files_mm, input_not_in_output_mm, output_not_in_input_mm
    )
    assert_import_segmentation_used(output)

def test_transcript_assignment(run_component, random_path):
    transcript_assignment = input_ta + "segmentation.csv"
    viz_polygons = input_ta + "segmentation_polygons_2d.json"
    output = random_path()
    run_component(
        [
            "--xenium_bundle",
            input_mm,
            "--id",
            id_mm,
            "--transcript_assignment",
            transcript_assignment,
            "--viz_polygons",
            viz_polygons,
            "--units",
            "microns",
            "--output",
            output,
        ]
    )

    assert_outputs_exists(input_mm, output, input_not_in_output_mm, output_not_in_input_mm)
    assert_valid_files(output)
    assert_identical(
        input_mm, output, changed_files_mm, input_not_in_output_mm, output_not_in_input_mm
    )
    assert_import_segmentation_used(output)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
