import sys

import pytest
import spatialdata as sd


def test_two_inputs(run_component, tmp_path):
    zarr_path = meta["resources_dir"] + "/xenium_tiny.zarr"
    output = tmp_path / "output.zarr"

    run_component(
        [
            "--inputs",
            zarr_path,
            "--inputs",
            zarr_path,
            "--output",
            str(output),
        ]
    )

    assert output.is_dir(), "Output Zarr store was not created"

    result = sd.read_zarr(output)
    source = sd.read_zarr(zarr_path)

    # Each spatial element should appear twice, once per dataset name
    # (xenium_tiny and xenium_tiny_2)
    for element_type, elements in [
        ("images", source.images),
        ("labels", source.labels),
        ("points", source.points),
        ("shapes", source.shapes),
    ]:
        result_keys = set(getattr(result, element_type).keys())
        for key in elements.keys():
            assert f"{key}-xenium_tiny" in result_keys, (
                f"Expected '{key}-xenium_tiny' in {element_type}, got {result_keys}"
            )
            assert f"{key}-xenium_tiny_2" in result_keys, (
                f"Expected '{key}-xenium_tiny_2' in {element_type}, got {result_keys}"
            )

    # Table should have 2× the original obs count
    assert "table" in result.tables, "Main table not found in output"
    assert result["table"].n_obs == source["table"].n_obs * 2, (
        f"Expected {source['table'].n_obs * 2} observations, got {result['table'].n_obs}"
    )
    assert result["table"].n_vars == source["table"].n_vars, (
        f"Expected {source['table'].n_vars} variables, got {result['table'].n_vars}"
    )


def test_anndata_label(run_component, tmp_path):
    """--anndata_label should add a batch column with one value per input."""
    zarr_path = meta["resources_dir"] + "/xenium_tiny.zarr"
    output = tmp_path / "output.zarr"

    run_component(
        [
            "--inputs",
            zarr_path,
            "--inputs",
            zarr_path,
            "--output",
            str(output),
            "--anndata_label",
            "batch",
        ]
    )

    assert output.is_dir(), "Output Zarr store was not created"

    result = sd.read_zarr(output)
    table = result["table"]

    assert "batch" in table.obs.columns, "Batch label column not found in obs"
    assert len(table.obs["batch"].unique()) == 2, (
        f"Expected 2 batch values, got {table.obs['batch'].unique()}"
    )


def test_single_input_passthrough(run_component, tmp_path):
    """A single input should be written through unchanged."""
    zarr_path = meta["resources_dir"] + "/xenium_tiny.zarr"
    output = tmp_path / "output.zarr"

    run_component(
        [
            "--inputs",
            zarr_path,
            "--output",
            str(output),
        ]
    )

    assert output.is_dir(), "Output Zarr store was not created"

    source = sd.read_zarr(zarr_path)
    result = sd.read_zarr(output)

    assert result["table"].n_obs == source["table"].n_obs, (
        "Single-input passthrough changed the table obs count"
    )
    assert result["table"].n_vars == source["table"].n_vars, (
        "Single-input passthrough changed the table var count"
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
