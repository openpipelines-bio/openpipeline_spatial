import json
import os
import sys
from pathlib import Path
import filecmp
import subprocess
import pytest
import re
import random
import string

import pandas as pd
import scanpy as sc

## VIASH START
meta = {"name": "xeniumranger_resegment", "resources_dir": "resources_test/xenium"}
## VIASH END

input = meta["resources_dir"] + "/xenium_tiny/"
id = "xeniun_tiny_resegment"
boundary_stain = "ATP1A1/CD45/E-Cadherin"
interior_stain = "18S"
dapi_filter = 5
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
]

def assert_outputs_exists(input, output):
    assert Path(output).is_dir() and len(list(Path(output).iterdir())) > 0, (
        "Output directory exists and is non-empty"
    )

    expected_output_dirs = ["morphology_focus"]

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

def _read_experiment_xenium_pair(input, output):
    with open(Path(input) / "experiment.xenium") as f:
        input_exp = json.load(f)
    with open(Path(output) / "experiment.xenium") as f:
        output_exp = json.load(f)
    return input_exp, output_exp


def _read_metrics_summary_pair(input, output):
    input_metrics = pd.read_csv(Path(input) / "metrics_summary.csv").iloc[0]
    output_metrics = pd.read_csv(Path(output) / "metrics_summary.csv").iloc[0]
    return input_metrics, output_metrics


def _read_analysis_summary_html_pair(input, output):
    input_html = (Path(input) / "analysis_summary.html").read_text()
    output_html = (Path(output) / "analysis_summary.html").read_text()
    return input_html, output_html


def test_basic_execution(run_component, random_path):
    output = random_path()
    run_component(
        [
            "--xenium_bundle",
            input,
            "--id",
            id,
            "--boundary_stain",
            boundary_stain, 
            "--interior_stain",
            interior_stain,
            "--segment_large_cells"
            "--expansion_distance", 
            "--dapi-filter", 
            dapi_filter, 
            "--resegment_nuclei",
            "--output",
            output,
        ]
    )

    assert_outputs_exists(input, output)
    assert_valid_files(output)
    assert_identical(input, output, changed_files)
