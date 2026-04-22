import pytest
import sys
import spatialdata as sd


def test_move_x_to_layer(run_component, tmp_path):
    """Test moving X matrix to layer."""
    input_path = meta["resources_dir"] + "/xenium_tiny.zarr"
    output_path = tmp_path / "output.zarr"

    run_component(
        [
            "--input",
            input_path,
            "--output",
            output_path,
            "--output_layer",
            "layer_1",
        ]
    )

    # Verify output
    assert output_path.exists(), "output zarr was not created"
    result = sd.read_zarr(str(output_path))

    assert result["table"].X is None, "X matrix should be None after moving"
    assert "layer_1" in result["table"].layers, "layer_1 should exist"
    assert result["table"].layers["layer_1"] is not None, "layer_1 should have data"


def test_move_named_layer_to_x(run_component, tmp_path):
    """Test moving a named layer to X matrix."""
    input_path = meta["resources_dir"] + "/xenium_tiny.zarr"
    output_path = tmp_path / "output.zarr"

    sdata = sd.read_zarr(input_path)
    sdata["table"].layers["test_layer"] = sdata["table"].X.copy()
    sdata["table"].X = None
    tmp_input = tmp_path / "input_with_layer.zarr"
    sdata.write(str(tmp_input))

    run_component(
        [
            "--input",
            str(tmp_input),
            "--output",
            str(output_path),
            "--input_layer",
            "test_layer",
        ]
    )

    # Verify output
    assert output_path.exists(), "output zarr was not created"
    result = sd.read_zarr(str(output_path))

    assert "test_layer" not in result["table"].layers, "test_layer should be removed"
    assert result["table"].X is not None, "X matrix should have data"


def test_move_layer_to_layer(run_component, tmp_path):
    """Test moving a named layer to another named layer."""
    input_path = meta["resources_dir"] + "/xenium_tiny.zarr"
    output_path = tmp_path / "output.zarr"

    sdata = sd.read_zarr(input_path)
    sdata["table"].layers["test_layer"] = sdata["table"].X.copy()
    tmp_input = tmp_path / "input_with_layer.zarr"
    sdata.write(str(tmp_input))

    run_component(
        [
            "--input",
            str(tmp_input),
            "--output",
            str(output_path),
            "--input_layer",
            "test_layer",
            "--output_layer",
            "layer_2",
        ]
    )

    # Verify output
    assert output_path.exists(), "output zarr was not created"
    result = sd.read_zarr(str(output_path))

    assert "test_layer" not in result["table"].layers, "test_layer should be removed"
    assert "layer_2" in result["table"].layers, "layer_2 should exist"
    assert result["table"].layers["layer_2"] is not None, "layer_2 should have data"
