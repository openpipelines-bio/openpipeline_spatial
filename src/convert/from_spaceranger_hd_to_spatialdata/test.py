import sys
import pytest
import spatialdata as sd

## VIASH START
meta = {
    "executable": "./target/executable/convert/from_spaceranger_hd_to_spatialdata/from_spaceranger_hd_to_spatialdata",
    "resources_dir": "resources_test/",
    "config": "src/convert/from_spaceranger_hd_to_spatialdata/config.vsh.yaml",
}
## VIASH END

input = f"{meta['resources_dir']}/visium_hd/Visium_HD_Mouse_Brain_tiny_spaceranger"


def test_segmented_cells(run_component, tmp_path):
    output = tmp_path / "output.zarr"
    run_component(
        [
            "--input",
            input,
            "--output",
            str(output),
            "--dataset_id",
            "Visium_HD_Mouse_Brain",
        ]
    )

    assert output.is_dir(), "output zarr store was not created"

    sdata = sd.read_zarr(output)

    assert "table" in sdata.tables, (
        f"expected a 'table' with segmented cells, got {list(sdata.tables)}"
    )
    table = sdata.tables["table"]
    assert table.n_obs > 0, "cell table has no cells"

    # Cell boundary polygons are stored as a shapes element.
    assert any("cell_segmentations" in key for key in sdata.shapes), (
        f"expected cell segmentation shapes, got {list(sdata.shapes)}"
    )


def test_bins(run_component, tmp_path):
    output = tmp_path / "bins.zarr"
    run_component(
        [
            "--input",
            input,
            "--output",
            str(output),
            "--mode",
            "bins",
            "--bin_size",
            "8",
            "--dataset_id",
            "Visium_HD_Mouse_Brain",
        ]
    )

    assert output.is_dir(), "output zarr store was not created"

    sdata = sd.read_zarr(output)

    assert "table" in sdata.tables, (
        f"expected a 'table' with binned counts, got {list(sdata.tables)}"
    )
    assert sdata.tables["table"].n_obs > 0, "binned table has no bins"

    # Binned data is stored as square-bin shapes named after the bin size.
    assert any("square_008um" in key for key in sdata.shapes), (
        f"expected 8um bin shapes, got {list(sdata.shapes)}"
    )


def test_output_layer(run_component, tmp_path):
    output = tmp_path / "output_layer.zarr"
    run_component(
        [
            "--input",
            input,
            "--output",
            str(output),
            "--dataset_id",
            "Visium_HD_Mouse_Brain",
            "--output_layer",
            "counts",
        ]
    )

    sdata = sd.read_zarr(output)
    assert "counts" in sdata.tables["table"].layers, (
        "expected counts to be stored as a named layer"
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
