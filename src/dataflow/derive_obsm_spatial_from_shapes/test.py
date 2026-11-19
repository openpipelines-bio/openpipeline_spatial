import subprocess
import sys

import geopandas as gpd
import numpy as np
import pytest
import spatialdata as sd
from spatialdata.models import ShapesModel, TableModel

## VIASH START
meta = {
    "executable": "./target/executable/dataflow/derive_obsm_spatial_from_shapes/derive_obsm_spatial_from_shapes",
    "resources_dir": "resources_test/",
    "config": "src/dataflow/derive_obsm_spatial_from_shapes/config.vsh.yaml",
}
## VIASH END

input = f"{meta['resources_dir']}/xenium_tiny.zarr"


def _annotate_shapes(sdata, region):
    """Make the table annotate the shapes element `region` by `cell_id` and
    remove .obsm['spatial'], simulating a reader that never populates it."""
    adata = sdata.tables["table"].copy()
    del adata.obsm["spatial"]
    adata.obs["region"] = region
    adata.obs["region"] = adata.obs["region"].astype("category")
    sdata.tables["table"] = TableModel.parse(
        adata,
        region=region,
        region_key="region",
        instance_key="cell_id",
        overwrite_metadata=True,
    )
    return sdata


def _cell_boundaries_input(tmp_path):
    """Write an input where the table annotates the polygon cell_boundaries."""
    sdata = sd.read_zarr(input, selection=("tables", "shapes"))
    assert (sdata.shapes["cell_boundaries"].geom_type == "Polygon").all()
    sdata = _annotate_shapes(sdata, "cell_boundaries")

    input_zarr = tmp_path / "input.zarr"
    sdata.write(input_zarr)
    return input_zarr


def test_derives_spatial_from_polygon_centroids(run_component, tmp_path):
    input_zarr = _cell_boundaries_input(tmp_path)
    output = tmp_path / "output.zarr"
    run_component(["--input", str(input_zarr), "--output", str(output)])

    sdata = sd.read_zarr(output)
    table = sdata.tables["table"]
    assert "spatial" in table.obsm, ".obsm['spatial'] was not added"

    coords = table.obsm["spatial"]
    assert coords.shape == (table.n_obs, 2)
    assert not np.isnan(coords).any(), ".obsm['spatial'] has missing values"

    # The derived coordinates are the shapes' centroids, in the row order of
    # the table's instance key
    centroids = sdata.shapes["cell_boundaries"].loc[table.obs["cell_id"]].centroid
    assert np.allclose(coords, np.column_stack([centroids.x, centroids.y]))


def test_derived_values_match_original_obsm_spatial(run_component, tmp_path):
    sdata = sd.read_zarr(input, selection=("tables", "shapes"))
    table = sdata.tables["table"]
    original_spatial = table.obsm["spatial"].copy()

    # Make circle shapes centred on the original coordinates, as
    # spatialdata_io.xenium() does for cell_circles
    circles = gpd.GeoDataFrame(
        {"radius": np.ones(table.n_obs)},
        geometry=gpd.points_from_xy(original_spatial[:, 0], original_spatial[:, 1]),
        index=table.obs["cell_id"].to_numpy(),
    )
    sdata.shapes["cell_circles"] = ShapesModel.parse(circles)
    sdata = _annotate_shapes(sdata, "cell_circles")

    input_zarr = tmp_path / "input.zarr"
    sdata.write(input_zarr)

    output = tmp_path / "output.zarr"
    run_component(["--input", str(input_zarr), "--output", str(output)])

    derived = sd.read_zarr(output).tables["table"].obsm["spatial"]
    assert np.allclose(derived, original_spatial), (
        "derived centroids are not the same as the original .obsm['spatial']"
    )


def test_output_same_as_input(run_component, tmp_path):
    input_zarr = _cell_boundaries_input(tmp_path)
    run_component(["--input", str(input_zarr), "--output", str(input_zarr)])

    sdata = sd.read_zarr(input_zarr)
    assert "spatial" in sdata.tables["table"].obsm, ".obsm['spatial'] was not added"
    assert "cell_boundaries" in sdata.shapes, (
        f"expected shapes 'cell_boundaries' to be kept, got {list(sdata.shapes)}"
    )


def test_output_obsm_spatial(run_component, tmp_path):
    input_zarr = _cell_boundaries_input(tmp_path)
    output = tmp_path / "output.zarr"
    run_component(
        [
            "--input",
            str(input_zarr),
            "--output",
            str(output),
            "--output_obsm_spatial",
            "cell_centroid",
        ]
    )

    table = sd.read_zarr(output).tables["table"]
    assert "cell_centroid" in table.obsm, ".obsm['cell_centroid'] was not added"
    assert "spatial" not in table.obsm, ".obsm['spatial'] should not be added"


def test_input_table(run_component, tmp_path):
    sdata = sd.read_zarr(input, selection=("tables", "shapes"))
    sdata = _annotate_shapes(sdata, "cell_boundaries")
    sdata.tables["cells"] = sdata.tables["table"]
    del sdata.tables["table"]

    input_zarr = tmp_path / "input.zarr"
    sdata.write(input_zarr)

    output = tmp_path / "output.zarr"
    run_component(
        [
            "--input",
            str(input_zarr),
            "--output",
            str(output),
            "--input_table",
            "cells",
        ]
    )

    sdata = sd.read_zarr(output)
    assert list(sdata.tables) == ["cells"], (
        f"expected only table 'cells', got {list(sdata.tables)}"
    )
    assert "spatial" in sdata.tables["cells"].obsm, ".obsm['spatial'] was not added"


def test_missing_input_table_raises(run_component, tmp_path):
    output = tmp_path / "output.zarr"

    with pytest.raises(subprocess.CalledProcessError) as exc_info:
        run_component(
            [
                "--input",
                input,
                "--output",
                str(output),
                "--input_table",
                "missing",
            ]
        )

    assert b"Table 'missing' not found" in exc_info.value.output
    assert b"Available tables: ['table']" in exc_info.value.output


def test_region_not_shapes_raises(run_component, tmp_path):
    # The table in the test data annotates the labels element 'cell_labels'
    output = tmp_path / "output.zarr"

    with pytest.raises(subprocess.CalledProcessError) as exc_info:
        run_component(["--input", input, "--output", str(output)])

    assert b"Region 'cell_labels'" in exc_info.value.output
    assert b"is not a shapes element" in exc_info.value.output


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
