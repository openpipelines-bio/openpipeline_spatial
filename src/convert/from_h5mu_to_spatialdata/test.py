import pytest
import sys
import spatialdata as sd
import mudata as mu
import os

def test_simple_execution(run_component, tmp_path):
    input = meta["resources_dir"] + "/xenium_tiny.h5mu"
    input_spatialdata = meta["resources_dir"] + "/xenium_tiny.zarr"
    output = tmp_path / "output.zarr"

    run_component(
        [
            "--input",
            input,
            "--input_spatialdata",
            input_spatialdata,
            "--output",
            output,
        ]
    )
    assert os.path.exists(output), "output zarr was not created"
    
    mdata = mu.read_h5mu(input)
    mod = mdata.mod["rna"]

    sdata = sd.read_zarr(output)
    table = sdata["table"]

    # Check that the main table in the SpatialData object matches the selected modality in the MuData object
    assert table.n_obs == mod.n_obs, "The number of observations in the SpatialData table does not match the selected modality in the MuData object."
    assert table.obs_names.equals(mod.obs_names), "The observation names in the SpatialData table do not match the selected modality in the MuData object."
    assert table.n_vars == mod.n_vars, "The number of variables in the SpatialData table does not match the selected modality in the MuData object."
    assert table.var_names.equals(mod.var_names), "The variable names in the SpatialData table do not match the selected modality in the MuData object."
    assert table.obs.keys().equals(mod.obs.keys()), "The observation metadata columns in the SpatialData table do not match those in the selected modality of the MuData object."
    assert table.var.keys().equals(mod.var.keys()), "The variable metadata columns in the SpatialData table do not match those in the selected modality of the MuData object."
    assert set(table.uns.keys()) == set(mod.uns.keys()), "The unstructured metadata keys in the SpatialData table do not match those in the selected modality of the MuData object."
    assert set(table.obsm.keys()) == set(mod.obsm.keys()), "The obsm keys in the SpatialData table do not match those in the selected modality of the MuData object."
    assert set(table.varm.keys()) == set(mod.varm.keys()), "The varm keys in the SpatialData table do not match those in the selected modality of the MuData object."

    sdata_existing = sd.read_zarr(input_spatialdata)

    # Check that remaining slots in the SpatialData object are filled from the existing SpatialData when provided
    assert sdata.images.keys() == sdata_existing.images.keys(), "The image keys in the output SpatialData do not match those in the existing SpatialData."
    assert sdata.labels.keys() == sdata_existing.labels.keys(), "The label keys in the output SpatialData do not match those in the existing SpatialData."
    assert sdata.shapes.keys() == sdata_existing.shapes.keys(), "The shape keys in the output SpatialData do not match those in the existing SpatialData."
    assert sdata.points.keys() == sdata_existing.points.keys(), "The point keys in the output SpatialData do not match those in the existing SpatialData."


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
