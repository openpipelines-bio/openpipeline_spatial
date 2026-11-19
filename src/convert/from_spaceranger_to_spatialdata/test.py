import sys

import numpy as np
import pytest
import spatialdata as sd

## VIASH START
meta = {
    "executable": "./target/executable/convert/from_spaceranger_to_spatialdata/from_spaceranger_to_spatialdata",
    "resources_dir": "resources_test/",
    "config": "src/convert/from_spaceranger_to_spatialdata/config.vsh.yaml",
}
## VIASH END

input = f"{meta['resources_dir']}/Visium_FFPE_Human_Ovarian_Cancer_tiny_spaceranger"


def test_simple_execution(run_component, tmp_path):
    output = tmp_path / "output.zarr"
    run_component(["--input", input, "--output", str(output)])

    assert output.is_dir(), "output zarr store was not created"

    sdata = sd.read_zarr(output)

    assert "table" in sdata.tables, (
        f"expected a table named 'table', got {list(sdata.tables)}"
    )
    table = sdata.tables["table"]
    assert table.n_obs > 0, "table has no spots"
    assert table.n_vars > 0, "table has no genes"
    assert table.X is not None, "table .X has no counts"


def test_output_type_raw(run_component, tmp_path):
    filtered = tmp_path / "filtered.zarr"
    run_component(
        ["--input", input, "--output", str(filtered), "--output_type", "filtered"]
    )

    raw = tmp_path / "raw.zarr"
    run_component(["--input", input, "--output", str(raw), "--output_type", "raw"])

    n_filtered = sd.read_zarr(filtered).tables["table"].n_obs
    n_raw = sd.read_zarr(raw).tables["table"].n_obs

    assert n_raw > n_filtered, (
        f"raw matrix ({n_raw} barcodes) should have more spots than "
        f"filtered matrix ({n_filtered} barcodes)"
    )


def test_output_layer(run_component, tmp_path):
    output = tmp_path / "output_layer.zarr"
    run_component(
        ["--input", input, "--output", str(output), "--output_layer", "counts"]
    )

    table = sd.read_zarr(output).tables["table"]

    assert "counts" in table.layers, "layer 'counts' not found"
    assert np.allclose(table.layers["counts"].toarray(), table.X.toarray()), (
        "layer 'counts' is not the same as .X"
    )


def test_output_table(run_component, tmp_path):
    output = tmp_path / "output_table.zarr"
    run_component(
        [
            "--input",
            input,
            "--output",
            str(output),
            "--output_table",
            "spots",
            "--output_layer",
            "counts",
        ]
    )

    sdata = sd.read_zarr(output)

    assert list(sdata.tables) == ["spots"], (
        f"expected only table 'spots', got {list(sdata.tables)}"
    )
    table = sdata.tables["spots"]
    assert table.n_obs > 0, "table has no spots"
    assert "counts" in table.layers, "layer 'counts' not found in table 'spots'"


def test_dataset_id(run_component, tmp_path):
    output = tmp_path / "output_dataset_id.zarr"
    run_component(
        ["--input", input, "--output", str(output), "--dataset_id", "my_sample"]
    )

    sdata = sd.read_zarr(output)

    assert "my_sample" in sdata.shapes, (
        f"expected shapes 'my_sample', got {list(sdata.shapes)}"
    )
    for image in ["my_sample_hires_image", "my_sample_lowres_image"]:
        assert image in sdata.images, (
            f"expected image '{image}', got {list(sdata.images)}"
        )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
