import json
import subprocess
import sys
import zipfile

import h5py
import pandas as pd
import pytest
import tifffile
import zarr


def _read_zarr_zip(path):
    import tempfile
    from pathlib import Path

    with tempfile.TemporaryDirectory() as tmpdir:
        with zipfile.ZipFile(path) as zf:
            zf.extractall(tmpdir)
        z = zarr.open(Path(tmpdir), mode="r")
        return {
            "number_cells": int(z.attrs["number_cells"]),
            "n_cell_id_rows": z["cell_id"].shape[0],
            "n_cell_summary_rows": z["cell_summary"].shape[0],
        }


def test_basic_crop(run_component, tmp_path):
    input_dir = meta["resources_dir"] + "/xenium_multicellseg_tiny"
    output_dir = tmp_path / "output"

    run_component(
        [
            "--input",
            input_dir,
            "--output",
            str(output_dir),
        ]
    )

    assert output_dir.is_dir(), "Output bundle was not created"

    for fname in [
        "cells.parquet",
        "cell_boundaries.parquet",
        "nucleus_boundaries.parquet",
        "transcripts.parquet",
        "cell_feature_matrix.h5",
        "cells.zarr.zip",
        "experiment.xenium",
        "metrics_summary.csv",
    ]:
        assert (output_dir / fname).exists(), f"{fname} missing from output bundle"

    source_cells = pd.read_parquet(f"{input_dir}/cells.parquet")
    cells = pd.read_parquet(output_dir / "cells.parquet")
    assert 0 < len(cells) <= len(source_cells), (
        "Cropped cells.parquet should be non-empty and no larger than the source"
    )

    with h5py.File(output_dir / "cell_feature_matrix.h5") as f:
        shape = f["matrix/shape"][:]
        assert shape[1] == len(cells), (
            "cell_feature_matrix.h5 barcode count should match cells.parquet row count"
        )
        barcodes = [b.decode() for b in f["matrix/barcodes"][:]]
        assert barcodes == cells["cell_id"].tolist(), (
            "cell_feature_matrix.h5 barcodes should match cells.parquet cell_id values, in order"
        )

    zarr_info = _read_zarr_zip(output_dir / "cells.zarr.zip")
    assert zarr_info["n_cell_id_rows"] == len(cells), (
        "cells.zarr.zip cell_id row count should match cells.parquet"
    )
    assert zarr_info["n_cell_summary_rows"] == len(cells), (
        "cells.zarr.zip cell_summary row count should match cells.parquet"
    )
    assert zarr_info["number_cells"] == len(cells)

    cell_boundaries = pd.read_parquet(output_dir / "cell_boundaries.parquet")
    assert set(cell_boundaries["cell_id"]).issubset(set(cells["cell_id"])), (
        "cell_boundaries.parquet should only reference kept cells"
    )

    morphology_dir = output_dir / "morphology_focus"
    channel_files = sorted(morphology_dir.glob("*.ome.tif"))
    assert len(channel_files) == 4, "Expected 4 cropped morphology_focus channel files"
    shapes = set()
    for tif_path in channel_files:
        image = tifffile.imread(tif_path, key=0)
        shapes.add(image.shape)
        ome_xml = tifffile.tiffcomment(tif_path)
        assert f'SizeX="{image.shape[1]}"' in ome_xml
        assert f'SizeY="{image.shape[0]}"' in ome_xml
    assert len(shapes) == 1, "All channel files should be cropped to the same shape"

    with open(output_dir / "experiment.xenium") as f:
        experiment = json.load(f)
    assert "pixel_size" in experiment


def test_smaller_patch_with_tighter_eps(run_component, tmp_path):
    """A tighter clustering radius should split the dense fixture into smaller sub-patches."""
    input_dir = meta["resources_dir"] + "/xenium_multicellseg_tiny"
    default_output = tmp_path / "default_output"
    tight_output = tmp_path / "tight_output"

    run_component(
        [
            "--input",
            input_dir,
            "--output",
            str(default_output),
        ]
    )
    run_component(
        [
            "--input",
            input_dir,
            "--output",
            str(tight_output),
            "--eps_um",
            "8",
        ]
    )

    default_cells = pd.read_parquet(default_output / "cells.parquet")
    tight_cells = pd.read_parquet(tight_output / "cells.parquet")
    assert len(tight_cells) < len(default_cells), (
        "A tighter clustering radius should select a smaller sub-patch"
    )


def test_fail_patch_rank_out_of_range(run_component, tmp_path):
    input_dir = meta["resources_dir"] + "/xenium_multicellseg_tiny"
    output_dir = tmp_path / "output"

    with pytest.raises(subprocess.CalledProcessError):
        run_component(
            [
                "--input",
                input_dir,
                "--output",
                str(output_dir),
                "--patch_rank",
                "1000",
            ]
        )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
