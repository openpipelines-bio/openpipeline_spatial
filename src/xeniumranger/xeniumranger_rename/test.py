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
meta = {"name": "xeniumranger_rename", "resources_dir": "resources_test/xenium"}
## VIASH END

input = meta["resources_dir"] + "/xenium_tiny/"
id = "xeniun_tiny_rename"
input_region_name = "region"
input_cassette_name = "cassette"
output_region_name = "region_renamed"
output_cassette_name = "cassette_renamed"


def assert_outputs_exists(input, output):
    assert Path(output).is_dir() and len(list(Path(output).iterdir())) > 0, (
        "Output exists and is non-empty"
    )

    input_files = sorted(
        f for f in os.listdir(input) if os.path.isfile(os.path.join(input, f))
    )
    output_files = sorted(
        f for f in os.listdir(output) if os.path.isfile(os.path.join(output, f))
    )

    assert input_files == output_files, (
        "Input and output directories should contain the same (type of) files"
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

def assert_identical(input, output, skip_files = ["experiment.xenium", "metrics_summary.csv", "analysis_summary.html"]):
    input_files = sorted(
        f for f in os.listdir(input) if os.path.isfile(os.path.join(input, f))
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


def assert_region_renaming(input, output):
    input_exp, output_exp = _read_experiment_xenium_pair(input, output)
    input_metrics, output_metrics = _read_metrics_summary_pair(input, output)
    input_html, output_html = _read_analysis_summary_html_pair(input, output)

    assert input_exp["region_name"] == input_region_name, (
        "experiment.xenium in the input bundle should carry the original region_name"
    )
    assert output_exp["region_name"] == output_region_name, (
        "experiment.xenium in the output bundle should carry the new region_name"
    )
    assert input_metrics["region_name"] == input_region_name, (
        "metrics_summary.csv in the input bundle should carry the original region_name"
    )
    assert output_metrics["region_name"] == output_region_name, (
        "metrics_summary.csv in the output bundle should carry the new region_name"
    )
    assert output_region_name not in input_html, (
        "analysis_summary.html in the input bundle should not yet contain the new region_name"
    )
    assert output_region_name in output_html, (
        "analysis_summary.html in the output bundle should contain the new region_name"
    )


def assert_cassette_renaming(input, output):
    input_exp, output_exp = _read_experiment_xenium_pair(input, output)
    input_metrics, output_metrics = _read_metrics_summary_pair(input, output)
    input_html, output_html = _read_analysis_summary_html_pair(input, output)

    assert input_exp["cassette_name"] == input_cassette_name, (
        "experiment.xenium in the input bundle should carry the original cassette_name"
    )
    assert output_exp["cassette_name"] == output_cassette_name, (
        "experiment.xenium in the output bundle should carry the new cassette_name"
    )
    assert input_metrics["cassette_name"] == input_cassette_name, (
        "metrics_summary.csv in the input bundle should carry the original cassette_name"
    )
    assert output_metrics["cassette_name"] == output_cassette_name, (
        "metrics_summary.csv in the output bundle should carry the new cassette_name"
    )
    assert output_cassette_name not in input_html, (
        "analysis_summary.html in the input bundle should not yet contain the new cassette_name"
    )
    assert output_cassette_name in output_html, (
        "analysis_summary.html in the output bundle should contain the new cassette_name"
    )


def test_basic_execution_both(run_component, random_path):
    output = random_path()
    run_component(
        [
            "--xenium_bundle",
            input,
            "--id",
            id,
            "--region_name",
            output_region_name,
            "--cassette_name",
            output_cassette_name,
            "--output",
            output,
        ]
    )

    assert_outputs_exists(input, output)
    assert_valid_files(output)
    assert_region_renaming(input, output)
    assert_cassette_renaming(input, output)
    assert_identical(input,output)


def test_basic_execution_region_only(run_component, random_path):
    output = random_path()
    run_component(
        [
            "--xenium_bundle",
            input,
            "--id",
            id,
            "--region_name",
            output_region_name,
            "--output",
            output,
        ]
    )

    assert_outputs_exists(input, output)
    assert_valid_files(output)
    assert_region_renaming(input, output)
    assert_identical(input,output)



def test_basic_execution_cassette_only(run_component, random_path):
    output = random_path()
    run_component(
        [
            "--xenium_bundle",
            input,
            "--id",
            id,
            "--cassette_name",
            output_cassette_name,
            "--output",
            output,
        ]
    )

    assert_outputs_exists(input, output)
    assert_valid_files(output)
    assert_cassette_renaming(input, output)
    assert_identical(input,output)

def test_no_name_given(run_component, random_path):
    output = random_path()
    run_component(
        [
            "--xenium_bundle",
            input,
            "--id",
            id + "_no_name_given",
            "--output",
            output,
        ]
    )

    assert_outputs_exists(input, output)
    assert_valid_files(output)
    assert_identical(input, output, skip_files=[])

def test_valid_id(run_component, random_path):
    output = random_path()
    malformed_id = ", ,"
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--xenium_bundle",
                input,
                "--id",
                malformed_id,
                "--output",
                output,
            ]
        )
    assert re.search(
        r"invalid id",
        err.value.stdout.decode("utf-8"),
        re.IGNORECASE,
    )

def test_length_region_name(run_component, random_path):
    output = random_path()
    long_region = "".join(random.choices(string.ascii_letters, k=65))
    with pytest.raises(subprocess.CalledProcessError):
        run_component(
            [
                "--xenium_bundle",
                input,
                "--id",
                id + "_bad_region",
                "--region_name",
                long_region,
                "--output",
                output,
            ]
        )

def test_length_cassette_name(run_component, random_path):
    output = random_path()
    long_cassette = "".join(random.choices(string.ascii_letters, k=33))
    with pytest.raises(subprocess.CalledProcessError):
        run_component(
            [
                "--xenium_bundle",
                input,
                "--id",
                id + "_bad_cassette",
                "--cassette_name",
                long_cassette,
                "--output",
                output,
            ]
        )

def test_valid_region(run_component, random_path):
    output = random_path()
    malformed_region = ", ,"
    with pytest.raises(subprocess.CalledProcessError):
        run_component(
            [
                "--xenium_bundle",
                input,
                "--id",
                id + "_malformed_region",
                "--region_name",
                malformed_region,
                "--output",
                output,
            ]
        )

def test_valid_cassette(run_component, random_path):
    output = random_path()
    malformed_cassette = ", ,"
    with pytest.raises(subprocess.CalledProcessError):
        run_component(
            [
                "--xenium_bundle",
                input,
                "--id",
                id + "_malformed_cassette",
                "--cassette_name",
                malformed_cassette,
                "--output",
                output,
            ]
        )

if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))