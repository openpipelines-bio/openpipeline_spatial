import spatialdata as sd
import sys
import pytest

## VIASH START
par = {"input": "input.zarr"}
meta = {"resources_dir": "resources_test"}
## VIASH END


def test_run():
    sdata = sd.read_zarr(par["input"])

    assert "table" in sdata.tables, "Output should contain a 'table' element."
    table = sdata.tables["table"]
    assert table.n_obs > 0, "Table should contain cells."

    # Cell boundary polygons from the segmentation.
    assert any("cell_segmentations" in key for key in sdata.shapes), (
        "Output should contain cell segmentation shapes (polygons)."
    )

    # Space Ranger cell-type annotations merged into the table's .obs.
    assert "fine_cell_type" in table.obs, (
        "Output table .obs should contain Space Ranger cell-type annotations."
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "--import-mode=importlib"]))
