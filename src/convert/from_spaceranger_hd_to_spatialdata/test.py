import re
import subprocess
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

    # Space Ranger cell-type annotations (segmented_outputs/cell_types/) are
    # merged into .obs.
    assert "fine_cell_type" in table.obs, (
        f"expected cell-type annotations in .obs, got {list(table.obs.columns)}"
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

    assert list(sdata.tables) == ["table"], (
        f"expected only a 'table' with binned counts, got {list(sdata.tables)}"
    )
    assert sdata.tables["table"].n_obs > 0, "binned table has no bins"

    # Binned data is stored as square-bin shapes named after the bin size.
    # The cell segmentation is not kept in bins mode.
    assert list(sdata.shapes) == ["Visium_HD_Mouse_Brain_square_008um"], (
        f"expected only 8um bin shapes, got {list(sdata.shapes)}"
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
            "--output_table",
            "cells",
        ]
    )

    sdata = sd.read_zarr(output)
    assert list(sdata.tables) == ["cells"], f"got tables {list(sdata.tables)}"
    assert "counts" in sdata.tables["cells"].layers, (
        "expected counts to be stored as a named layer"
    )


def test_bins_multiple_default_names(run_component, tmp_path):
    output = tmp_path / "bins_multiple.zarr"
    run_component(
        [
            "--input",
            input,
            "--output",
            str(output),
            "--mode",
            "bins",
            "--bin_size",
            "2",
            "--bin_size",
            "8",
            "--bin_size",
            "16",
            "--dataset_id",
            "Visium_HD_Mouse_Brain",
        ]
    )

    sdata = sd.read_zarr(output)

    # Several bin tables keep the spatialdata-io names.
    expected = ["square_002um", "square_008um", "square_016um"]
    assert sorted(sdata.tables) == expected, f"got tables {list(sdata.tables)}"
    for name in expected:
        region = sdata.tables[name].uns["spatialdata_attrs"]["region"]
        assert region == f"Visium_HD_Mouse_Brain_{name}"
        assert region in sdata.shapes, f"shapes '{region}' missing"


def test_bins_multiple_renamed(run_component, tmp_path):
    output = tmp_path / "bins_renamed.zarr"
    run_component(
        [
            "--input",
            input,
            "--output",
            str(output),
            "--mode",
            "bins",
            "--bin_size",
            "16",
            "--bin_size",
            "8",
            "--output_table",
            "bins_16",
            "--output_table",
            "bins_8",
            "--output_layer",
            "counts",
            "--dataset_id",
            "Visium_HD_Mouse_Brain",
        ]
    )

    sdata = sd.read_zarr(output)

    assert sorted(sdata.tables) == ["bins_16", "bins_8"], (
        f"got tables {list(sdata.tables)}"
    )
    # Names are matched to bin sizes by position: 8um has more bins than 16um.
    assert sdata.tables["bins_8"].n_obs > sdata.tables["bins_16"].n_obs
    assert (
        sdata.tables["bins_8"].uns["spatialdata_attrs"]["region"]
        == "Visium_HD_Mouse_Brain_square_008um"
    )
    for name in ["bins_16", "bins_8"]:
        assert "counts" in sdata.tables[name].layers, (
            f"expected counts layer in table '{name}'"
        )


def test_fail_output_table_count(run_component, tmp_path):
    output = tmp_path / "should_not_exist.zarr"

    with pytest.raises(subprocess.CalledProcessError) as err:
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
                "--bin_size",
                "16",
                "--output_table",
                "bins",
                "--dataset_id",
                "Visium_HD_Mouse_Brain",
            ]
        )
    assert re.search(
        r"'--output_table' has 1 names, but 2 tables will be loaded",
        err.value.stdout.decode("utf-8"),
    )
    assert not output.exists()


def test_fail_missing_bin_size(run_component, tmp_path):
    output = tmp_path / "should_not_exist.zarr"

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input,
                "--output",
                str(output),
                "--mode",
                "bins",
                "--bin_size",
                "4",
                "--dataset_id",
                "Visium_HD_Mouse_Brain",
            ]
        )
    assert re.search(
        r"Bin sizes \[4\] not found .* Available bin sizes: \[2, 8, 16\]",
        err.value.stdout.decode("utf-8"),
    )
    assert not output.exists()


def test_fail_nucleus_missing_files(run_component, tmp_path):
    output = tmp_path / "should_not_exist.zarr"

    # The test data has no nucleus segmentation files.
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input,
                "--output",
                str(output),
                "--load_nucleus_segmentations",
                "true",
                "--dataset_id",
                "Visium_HD_Mouse_Brain",
            ]
        )
    assert re.search(
        r"'--load_nucleus_segmentations' is set, but no nucleus segmentation was loaded",
        err.value.stdout.decode("utf-8"),
    )
    assert not output.exists()


def test_fail_nucleus_bins_mode(run_component, tmp_path):
    output = tmp_path / "should_not_exist.zarr"

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input,
                "--output",
                str(output),
                "--mode",
                "bins",
                "--load_nucleus_segmentations",
                "true",
                "--dataset_id",
                "Visium_HD_Mouse_Brain",
            ]
        )
    assert re.search(
        r"'--load_nucleus_segmentations' can only be used with '--mode segmented_cells'",
        err.value.stdout.decode("utf-8"),
    )
    assert not output.exists()


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
