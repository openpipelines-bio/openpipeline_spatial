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
meta = {"name": "xeniumranger_resegment", "resources_dir": "resources_test/xenium"}
## VIASH END

# 1. Nuclear expansion arguments, ne = nulcear expansion
id = "xeniun_tiny_resegment_ne"
boundary_stain_ne = "disable"
interior_stain_ne = "disable"

# 2. Multimodal segmentation arguments, mm = multimodal
id = "xeniun_tiny_resegment_mm"
boundary_stain_mm = "ATP1A1/CD45/E-Cadherin"
interior_stain_mm = "18S"


# 3. General arguments
dapi_filter = "15"
expansion_distance = "0"
changed_files = [
    "cells.csv.gz",
    "cells.parquet",
    "cells.zarr.zip",
    "cell_boundaries.csv.gz",
    "cell_boundaries.parquet",
    "nucleus_boundaries.csv.gz",
    "nucleus_boundaries.parquet",
    "cell_feature_matrix.h5",
    "cell_feature_matrix.tar.gz",
    "cell_feature_matrix.zarr.zip",
    "transcripts.parquet",
    "transcripts.zarr.zip",
    "analysis.tar.gz",
    "analysis.zarr.zip",
    "metrics_summary.csv",
    "analysis_summary.html",
    "experiment.xenium",
]


def assert_outputs_exists(input, output):
    assert Path(output).is_dir() and len(list(Path(output).iterdir())) > 0, (
        "Output directory exists and is non-empty"
    )

    expected_output_dirs = ["analysis", "cell_feature_matrix", "morphology_focus"]

    input_not_in_output = [
        "aux_outputs.tar.gz",
        "analysis.tar.gz",
        "cell_feature_matrix.tar.gz",
    ]

    input_files = sorted(
        f
        for f in os.listdir(input)
        if os.path.isfile(os.path.join(input, f)) and f not in input_not_in_output
    )
    output_files = sorted(
        f for f in os.listdir(output) if os.path.isfile(os.path.join(output, f))
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

    with open(Path(output) / "gene_panel.json") as f:
        panel_data = json.load(f)
    assert panel_data, "gene_panel.json should not be empty"

    return transcripts, cells, panel_data


def assert_stain_segmentation_used(output, boundary_stain, interior_stain):
    with open(Path(output) / "experiment.xenium") as f:
        exp = json.load(f)

    assert exp["segmentation_stain"], (
        "segmentation_stain should be non-empty when --boundary_stain/--interior_stain are set"
    )
    assert (
        exp["segmented_cell_boundary_frac"] > 0
        or exp["segmented_cell_interior_frac"] > 0
    ), (
        "stain-based segmentation fractions should be non-zero when a  boundary/interior stain is configured"
    )
    assert exp["segmented_cell_nuc_expansion_frac"] < 1.0, (
        "not every cell should fall back to pure nucleus-expansion when stain-based segmentation is active"
    )


def assert_identical(input, output, skip_files):
    input_not_in_output = [
        "aux_outputs.tar.gz",
        "analysis.tar.gz",
        "cell_feature_matrix.tar.gz",
    ]
    input_files = sorted(
        f
        for f in os.listdir(input)
        if os.path.isfile(os.path.join(input, f)) and f not in input_not_in_output
    )

    output_files = sorted(
        f for f in os.listdir(output) if os.path.isfile(os.path.join(output, f))
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
    input = meta["resources_dir"] + "/xenium_tiny/"
    output = random_path()
    run_component(
        [
            "--xenium_bundle",
            input,
            "--id",
            id,
            "--boundary_stain",
            boundary_stain_ne,
            "--interior_stain",
            interior_stain_ne,
            "--output",
            output,
        ]
    )

    assert_outputs_exists(input, output)
    assert_valid_files(output)
    assert_identical(input, output, changed_files)


def test_valid_id(run_component, random_path):
    input = meta["resources_dir"] + "/xenium_tiny/"
    output = random_path()
    malformed_id = ", ,"
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--xenium_bundle",
                input,
                "--id",
                malformed_id,
                "--boundary_stain",
                boundary_stain_ne,
                "--interior_stain",
                interior_stain_ne,
                "--output",
                output,
            ]
        )
    assert re.search(
        r"invalid value.*--id",
        err.value.stdout.decode("utf-8"),
        re.IGNORECASE,
    )


def test_resegment_nuclei(run_component, random_path):
    input = meta["resources_dir"] + "/xenium_tiny/"
    output = random_path()
    run_component(
        [
            "--xenium_bundle",
            input,
            "--id",
            id,
            "--boundary_stain",
            boundary_stain_ne,
            "--interior_stain",
            interior_stain_ne,
            "--resegment_nuclei",
            "--output",
            output,
        ]
    )

    assert_outputs_exists(input, output)
    assert_valid_files(output)
    assert_identical(input, output, changed_files)

    input_nb_parquet = pd.read_parquet(Path(input) / "nucleus_boundaries.parquet")
    output_nb_parquet = pd.read_parquet(Path(output) / "nucleus_boundaries.parquet")

    assert not input_nb_parquet.equals(output_nb_parquet), (
        "--resegment_nuclei should alter nuclear boundaries "
    )


def test_missing_file(run_component, random_path, tmp_path):
    output = random_path()
    input = meta["resources_dir"] + "/xenium_tiny/"
    incomplete_bundle = tmp_path / "incomplete_bundle"
    shutil.copytree(input, incomplete_bundle)
    os.remove(incomplete_bundle / "experiment.xenium")
    with pytest.raises(subprocess.CalledProcessError):
        run_component(
            [
                "--xenium_bundle",
                str(incomplete_bundle),
                "--id",
                id,
                "--boundary_stain",
                boundary_stain_ne,
                "--interior_stain",
                interior_stain_ne,
                "--output",
                output,
            ]
        )


# 2. Multimodal path
def test_multimodal_fixtures(run_component, random_path):
    input = meta["resources_dir"] + "/xenium_multicellseg_tiny_raw/"
    output = random_path()
    run_component(
        [
            "--xenium_bundle",
            input,
            "--id",
            id,
            "--boundary_stain",
            boundary_stain_mm,
            "--interior_stain",
            interior_stain_mm,
            "--output",
            output,
        ]
    )
    assert_outputs_exists(input, output)
    assert_valid_files(output)
    assert_identical(input, output, changed_files)
    assert_stain_segmentation_used(output, boundary_stain_mm, interior_stain_mm)


def test_segment_large_cells(run_component, random_path):
    input = meta["resources_dir"] + "/xenium_multicellseg_tiny_raw/"
    baseline_output = random_path()
    run_component(
        [
            "--xenium_bundle",
            input,
            "--id",
            id,
            "--boundary_stain",
            boundary_stain_mm,
            "--interior_stain",
            interior_stain_mm,
            "--output",
            baseline_output,
        ]
    )

    output = random_path()
    run_component(
        [
            "--xenium_bundle",
            input,
            "--id",
            id,
            "--boundary_stain",
            boundary_stain_mm,
            "--interior_stain",
            interior_stain_mm,
            "--segment_large_cells",
            "--output",
            output,
        ]
    )
    assert_outputs_exists(input, output)
    assert_valid_files(output)
    assert_identical(input, output, changed_files)
    assert_stain_segmentation_used(output, boundary_stain_mm, interior_stain_mm)

    baseline_cb_parquet = pd.read_parquet(
        Path(baseline_output) / "cell_boundaries.parquet"
    )
    output_cb_parquet = pd.read_parquet(Path(output) / "cell_boundaries.parquet")

    assert not baseline_cb_parquet.equals(output_cb_parquet), (
        "--segment_large_cells should alter cell boundaries relative to a resegment run without it"
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
