import pytest
import os
import sys
import spatialdata as sd
import subprocess


def test_simple_execution(run_component, tmp_path):
    output_sd_path = tmp_path / "sd"

    run_component(
        [
            "--input",
            meta["resources_dir"] + "/xenium_tiny",
            "--output",
            output_sd_path,
        ]
    )

    assert os.path.exists(output_sd_path), "Output zarr folder was not created"

    sdata = sd.read_zarr(output_sd_path)
    assert isinstance(sdata, sd.SpatialData), (
        "the generated output is not a SpatialData object"
    )

    assert os.path.exists(output_sd_path / "images"), "images folder was not created"
    assert os.path.exists(output_sd_path / "labels"), "labels folder was not created"
    assert os.path.exists(output_sd_path / "points"), "images folder was not created"
    assert os.path.exists(output_sd_path / "shapes"), "shapes folder was not created"
    assert os.path.exists(output_sd_path / "tables"), "tables folder was not created"
    assert (output_sd_path / "zarr.json").is_file(), "zarr metadata file was not created"


def test_compressed_input(run_component, tmp_path):
    output_sd_path = tmp_path / "sd"
    zipped_input = tmp_path / "xenium_tiny.zip"

    subprocess.run(
        ["zip", "-r", str(zipped_input), "xenium_tiny"],
        cwd=meta["resources_dir"],
        check=True,
    )
    run_component(
        [
            "--input",
            zipped_input,
            "--output",
            output_sd_path,
        ]
    )

    assert os.path.exists(output_sd_path), "Output zarr folder was not created"

    sdata = sd.read_zarr(output_sd_path)
    assert isinstance(sdata, sd.SpatialData), (
        "the generated output is not a SpatialData object"
    )

    assert os.path.exists(output_sd_path / "images"), "images folder was not created"
    assert os.path.exists(output_sd_path / "labels"), "labels folder was not created"
    assert os.path.exists(output_sd_path / "points"), "images folder was not created"
    assert os.path.exists(output_sd_path / "shapes"), "shapes folder was not created"
    assert os.path.exists(output_sd_path / "tables"), "tables folder was not created"
    assert (output_sd_path / "zarr.json").is_file(), "zarr metadata file was not created"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
