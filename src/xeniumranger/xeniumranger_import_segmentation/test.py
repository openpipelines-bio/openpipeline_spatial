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
meta = {
    "name": "xeniumranger_import_segmentation",
    "resources_dir": "resources_test/xenium_multichannel",
}
## VIASH END

# Raw XOA bundle cropped to a single patch of cells
input_bundle = meta["resources_dir"] + "/xenium_multicellseg_tiny_cropped/"
# Baysor-style transcript assignment + viz polygons derived from the cropped bundle (micron units)
input_ta = (
    meta["resources_dir"] + "/xenium_multicellseg_tiny_cropped_transcript_assignment/"
)

# General arguments
units = "pixels"
expansion_distance = "5"

changed_files = [
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
input_not_in_output = [
    ".end-of-run",
    "aux_outputs.tar.gz",
    "analysis.tar.gz",
    "cell_feature_matrix.tar.gz",
]
output_not_in_input = [
    "cell_id_map.csv.gz",
]


def assert_outputs_exists(
    input_bundle, output, input_not_in_output, output_not_in_input
):
    assert Path(output).is_dir() and len(list(Path(output).iterdir())) > 0, (
        "Output directory exists and is non-empty"
    )

    expected_output_dirs = ["analysis", "cell_feature_matrix", "morphology_focus"]

    input_files = sorted(
        f
        for f in os.listdir(input_bundle)
        if os.path.isfile(os.path.join(input_bundle, f))
        and f not in input_not_in_output
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


def assert_transcript_assignment_imported(output, transcript_assignment):
    assignment = pd.read_csv(transcript_assignment, keep_default_na=False)
    assigned = assignment[assignment["cell"] != ""]
    imported_cells = set(assigned["cell"])

    cell_id_map = pd.read_csv(Path(output) / "cell_id_map.csv.gz")
    assert set(cell_id_map["Imported cell ID"]) == imported_cells, (
        "cell_id_map.csv.gz should map exactly the imported cell IDs"
    )
    new_cells = set(cell_id_map["Xenium Ranger new cell ID"])
    assert len(new_cells) == len(imported_cells), (
        "Each imported cell should map to a unique Xenium Ranger cell ID"
    )

    cells = sc.read_10x_h5(Path(output) / "cell_feature_matrix.h5")
    assert cells.n_obs == len(imported_cells), (
        "cell_feature_matrix.h5 should contain one cell per imported cell"
    )
    assert set(cells.obs_names) == new_cells, (
        "cell_feature_matrix.h5 barcodes should match the cell_id_map output IDs"
    )

    transcripts = pd.read_parquet(Path(output) / "transcripts.parquet")
    assert (transcripts["cell_id"] != "UNASSIGNED").sum() == len(assigned), (
        "Only transcripts assigned in the input CSV should be assigned to cells"
    )

    with open(Path(output) / "experiment.xenium") as f:
        exp = json.load(f)
    assert exp["num_cells"] == len(imported_cells)
    assert exp["imported_cell_frac"] == 1.0, "All cells should be imported"


def assert_identical(
    input_bundle, output, skip_files, input_not_in_output, output_not_in_input
):
    input_files = sorted(
        f
        for f in os.listdir(input_bundle)
        if os.path.isfile(os.path.join(input_bundle, f))
        and f not in input_not_in_output
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
                os.path.join(input_bundle, f), os.path.join(output, f), shallow=False
            ), f"{f} should be identical between input and output"


def test_basic_execution(run_component, random_path):
    nuclei = input_bundle + "cells.zarr.zip"
    output = random_path()
    run_component(
        [
            "--xenium_bundle",
            input_bundle,
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

    assert_outputs_exists(
        input_bundle, output, input_not_in_output, output_not_in_input
    )
    assert_valid_files(output)
    assert_identical(
        input_bundle,
        output,
        changed_files,
        input_not_in_output,
        output_not_in_input,
    )
    assert_import_segmentation_used(output)


def test_missing_file(run_component, random_path, tmp_path):
    output = random_path()
    incomplete_bundle = tmp_path / "incomplete_bundle"
    shutil.copytree(input_bundle, incomplete_bundle)
    os.remove(incomplete_bundle / "experiment.xenium")
    nuclei = input_bundle + "cells.zarr.zip"
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--xenium_bundle",
                str(incomplete_bundle),
                "--nuclei",
                nuclei,
                "--output",
                output,
            ]
        )
    stdout = err.value.stdout.decode("utf-8")
    assert re.search(r"experiment\.xenium", stdout, re.IGNORECASE)


def test_cells_nuclei(run_component, random_path):
    nuclei = input_bundle + "cells.zarr.zip"
    cells = input_bundle + "cells.zarr.zip"
    output = random_path()
    run_component(
        [
            "--xenium_bundle",
            input_bundle,
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

    assert_outputs_exists(
        input_bundle, output, input_not_in_output, output_not_in_input
    )
    assert_valid_files(output)
    assert_identical(
        input_bundle,
        output,
        changed_files,
        input_not_in_output,
        output_not_in_input,
    )
    assert_import_segmentation_used(output)


def test_transcript_assignment(run_component, random_path):
    transcript_assignment = input_ta + "segmentation.csv"
    viz_polygons = input_ta + "segmentation_polygons_2d.json"
    output = random_path()
    run_component(
        [
            "--xenium_bundle",
            input_bundle,
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

    assert_outputs_exists(
        input_bundle, output, input_not_in_output, output_not_in_input
    )
    assert_valid_files(output)
    assert_identical(
        input_bundle,
        output,
        changed_files,
        input_not_in_output,
        output_not_in_input,
    )
    assert_import_segmentation_used(output)
    assert_transcript_assignment_imported(output, transcript_assignment)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
