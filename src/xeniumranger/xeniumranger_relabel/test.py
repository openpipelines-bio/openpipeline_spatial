import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import filecmp
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
import tifffile

## VIASH START
meta = {"name": "xeniumranger_relabel", "resources_dir": "resources_test"}
## VIASH END

input = meta["resources_dir"] + "/xenium/xenium_tiny/"
panel = meta["resources_dir"] + "/xenium/xenium_tiny/gene_panel.json"
id = "xeniun_tiny_relabel"


def _read_parquet_pair(input, output, filename, columns=None):
    original = pd.read_parquet(Path(input) / filename, columns=columns)
    relabeled = pd.read_parquet(Path(output) / filename, columns=columns)
    return original, relabeled


def _assert_same_keys(original, relabeled, key, label):
    assert len(original) == len(relabeled), (
        f"{label} row count should be unchanged by relabeling"
    )
    merged = original.merge(relabeled, on=key, suffixes=("_before", "_after"))
    assert len(merged) == len(original), (
        f"the same set of {key}s should be present in {label} before and after relabeling"
    )
    return merged


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


def modified_panel(panel_path, tmp_path, old_name="A1cf", new_name="A1cf_renamed"):
    with open(panel_path) as f:
        data = json.load(f)

    matches = [
        t for t in data["payload"]["targets"] if t["type"]["data"]["name"] == old_name
    ]
    assert len(matches) == 1, (
        f"expected exactly one target named {old_name!r}, found {len(matches)}"
    )
    matches[0]["type"]["data"]["name"] = new_name

    modified_path = tmp_path / "modified_gene_panel.json"
    with open(modified_path, "w") as f:
        json.dump(data, f)

    return modified_path


def assert_relabeling(input, output, old_name, new_name):
    original, relabeled = _read_parquet_pair(
        input, output, "transcripts.parquet", columns=["transcript_id", "feature_name"]
    )
    merged = _assert_same_keys(
        original, relabeled, "transcript_id", "transcripts.parquet"
    )

    was_old = merged["feature_name_before"] == old_name
    assert was_old.any(), (
        f"no transcripts were originally assigned to {old_name!r} -- nothing to check"
    )

    assert (merged.loc[was_old, "feature_name_after"] == new_name).all(), (
        f"all transcripts originally assigned to {old_name!r} should now carry {new_name!r}"
    )
    assert not (merged["feature_name_after"] == old_name).any(), (
        f"no transcript should still carry the old name {old_name!r}"
    )


def assert_geometry_unchanged(input, output):
    coord_columns = ["x_location", "y_location", "z_location", "qv"]
    original, relabeled = _read_parquet_pair(
        input, output, "transcripts.parquet", columns=["transcript_id"] + coord_columns
    )
    merged = _assert_same_keys(
        original, relabeled, "transcript_id", "transcripts.parquet"
    )
    for col in coord_columns:
        assert np.allclose(merged[f"{col}_before"], merged[f"{col}_after"]), (
            f"{col} should be unchanged by relabeling"
        )

    cell_columns = ["cell_id", "x_centroid", "y_centroid"]
    original_cells, relabeled_cells = _read_parquet_pair(
        input, output, "cells.parquet", columns=cell_columns
    )
    merged_cells = _assert_same_keys(
        original_cells, relabeled_cells, "cell_id", "cells.parquet"
    )
    for col in ["x_centroid", "y_centroid"]:
        assert np.allclose(
            merged_cells[f"{col}_before"], merged_cells[f"{col}_after"]
        ), f"{col} should be unchanged by relabeling"

    # cell_boundaries/nucleus_boundaries have no unique per-row key (many vertices share a
    # cell_id), so compare positionally instead of via merge -- segmentation isn't touched by
    # relabel, so row order is expected to carry through unchanged.
    for boundary_file in ["cell_boundaries.parquet", "nucleus_boundaries.parquet"]:
        original_boundaries, relabeled_boundaries = _read_parquet_pair(
            input, output, boundary_file
        )

        assert len(original_boundaries) == len(relabeled_boundaries), (
            f"{boundary_file} row count should be unchanged by relabeling"
        )
        assert (
            original_boundaries["cell_id"] == relabeled_boundaries["cell_id"]
        ).all(), f"{boundary_file} cell_id order should be unchanged by relabeling"
        assert np.allclose(
            original_boundaries["vertex_x"], relabeled_boundaries["vertex_x"]
        ), f"{boundary_file} vertex_x should be unchanged by relabeling"
        assert np.allclose(
            original_boundaries["vertex_y"], relabeled_boundaries["vertex_y"]
        ), f"{boundary_file} vertex_y should be unchanged by relabeling"


def assert_morphology_unchanged(input, output):
    morphology_files = [
        "morphology.ome.tif",
        "morphology_focus/morphology_focus_0000.ome.tif",
    ]
    for rel_path in morphology_files:
        original_path = Path(input) / rel_path
        relabeled_path = Path(output) / rel_path
        assert relabeled_path.is_file(), (
            f"{rel_path} should be present in relabeled output"
        )

        with tifffile.TiffFile(original_path) as original_tif:
            original_shapes = [s.shape for s in original_tif.series]
        with tifffile.TiffFile(relabeled_path) as relabeled_tif:
            relabeled_shapes = [s.shape for s in relabeled_tif.series]

        assert original_shapes == relabeled_shapes, (
            f"{rel_path} should have the same shape (pages x height x width) after relabeling"
        )


def test_basic_execution(run_component, random_path):
    output = random_path()
    run_component(
        [
            "--xenium_bundle",
            input,
            "--panel",
            panel,
            "--id",
            id,
            "--output",
            output,
        ]
    )

    assert_outputs_exists(input, output)
    assert_valid_files(output)


def test_relabelling(run_component, random_path, tmp_path):
    output = random_path()
    old_gene = "Defa5"
    new_gene = "Defa5_renamed"
    altered_panel = modified_panel(panel, tmp_path)

    run_component(
        [
            "--xenium_bundle",
            input,
            "--panel",
            altered_panel,
            "--id",
            id,
            "--output",
            output,
        ]
    )

    assert_outputs_exists(input, output)
    _, cells, _ = assert_valid_files(output)

    output_panel_path = Path(output) / "gene_panel.json"
    assert filecmp.cmp(altered_panel, output_panel_path, shallow=False), (
        "Output gene_panel.json should match the altered input gene panel content"
    )

    assert new_gene in cells.var_names, (
        "New gene name should be present in cell feature matrix"
    )
    assert old_gene not in cells.var_names, (
        "Old gene name should not be present in cell feature matrix"
    )

    assert_relabeling(input, output, old_gene, new_gene)

    assert_geometry_unchanged(input, output)

    assert_morphology_unchanged(input, output)


def test_valid_id(run_component, random_path):
    output = random_path()
    malformed_id = ", ,"
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--xenium_bundle",
                input,
                "--panel",
                panel,
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


def test_valid_panel(run_component, random_path, tmp_path):
    output = random_path()

    malformed_panel = tmp_path / "malformed_panel.json"
    malformed_data = {"gene": "not_a_gene", "transcript": "not_a_transcript"}
    with open(malformed_panel, "w") as file:
        json.dump(malformed_data, file)

    with pytest.raises(subprocess.CalledProcessError):
        run_component(
            [
                "--xenium_bundle",
                input,
                "--panel",
                malformed_panel,
                "--id",
                id,
                "--output",
                output,
            ]
        )


def test_incomplete_bundle(run_component, random_path, tmp_path):
    output = random_path()

    incomplete_bundle = tmp_path / "incomplete_bundle"
    shutil.copytree(input, incomplete_bundle)
    (incomplete_bundle / "cell_feature_matrix.h5").unlink()

    with pytest.raises(subprocess.CalledProcessError):
        run_component(
            [
                "--xenium_bundle",
                incomplete_bundle,
                "--panel",
                panel,
                "--id",
                id,
                "--output",
                output,
            ]
        )


def test_relative_paths(run_component, tmp_path, monkeypatch):
    work_dir = tmp_path / "workdir"
    work_dir.mkdir()
    shutil.copytree(input, work_dir / "bundle")
    monkeypatch.chdir(work_dir)

    run_component(
        [
            "--xenium_bundle",
            "bundle",
            "--panel",
            "bundle/gene_panel.json",
            "--id",
            id,
            "--output",
            "out",
        ]
    )

    assert_outputs_exists(work_dir / "bundle", work_dir / "out")


def test_repeated_id_isolation(run_component, random_path):
    output_first = random_path()
    output_second = random_path()

    # same --id run twice in one session: mktemp -d + cd in script.sh should give each
    # invocation its own pipestance dir, so the second run must not collide with the first
    run_component(
        [
            "--xenium_bundle",
            input,
            "--panel",
            panel,
            "--id",
            id,
            "--output",
            output_first,
        ]
    )
    run_component(
        [
            "--xenium_bundle",
            input,
            "--panel",
            panel,
            "--id",
            id,
            "--output",
            output_second,
        ]
    )

    assert_outputs_exists(input, output_first)
    assert_outputs_exists(input, output_second)
