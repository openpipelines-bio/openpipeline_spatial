import pytest
import sys
import mudata as mu
import subprocess

## VIASH START
meta = {
    "executable": "./target/executable/convert/from_cells2stats_to_h5mu/from_cells2stats_to_h5mu",
    "resources_dir": "resources_test/aviti/",
}
## VIASH END

input = f"{meta['resources_dir']}/aviti/teton_cells2stats_tiny/"


def test_simple_execution(run_component, tmp_path):
    output = tmp_path / "aviti.h5mu"

    # run component
    run_component(
        ["--input", input, "--output", str(output), "--output_compression", "gzip"]
    )

    assert output.is_file(), "output file was not created"

    mdata = mu.read_h5mu(output)
    assert list(mdata.mod.keys()) == ["rna"], "Expected modality rna"
    adata = mdata.mod["rna"]

    assert adata.X.dtype.kind == "f"
    expected_obs_keys = [
        "AreaUm",
        "Area",
        "Tile",
        "WellLabel",
        "Well",
        "Cell",
        "NuclearAreaUm",
        "NuclearArea",
    ]
    assert all([obs in expected_obs_keys for obs in adata.obs.columns])
    obs_counts = ["Area", "Cell", "NuclearArea"]
    assert all([adata.obs[obs].dtype.kind == "u" for obs in obs_counts])
    obs_areas = ["AreaUm", "NuclearAreaUm"]
    assert all([adata.obs[obs].dtype.kind == "f" for obs in obs_areas])
    obs_categories = ["Tile", "WellLabel", "Well"]
    assert all([adata.obs[obs].dtype.kind == "O" for obs in obs_categories])

    expected_obsm_keys = ["spatial", "spatial_um"]
    assert list(adata.obsm.keys()) == expected_obsm_keys
    assert list(adata.uns.keys()) == expected_obsm_keys
    assert all(adata.obsm[obsm].dtype.kind == "f" for obsm in expected_obsm_keys)


def test_compressed_input(run_component, tmp_path):
    output = tmp_path / "aviti.h5mu"
    zipped_input = tmp_path / "aviti.zip"

    subprocess.run(["zip", "-r", str(zipped_input), "."], cwd=input, check=True)

    # run component
    run_component(
        [
            "--input",
            zipped_input,
            "--output",
            str(output),
            "--output_compression",
            "gzip",
        ]
    )

    assert output.is_file(), "output file was not created"

    mdata = mu.read_h5mu(output)
    assert list(mdata.mod.keys()) == ["rna"], "Expected modality rna"
    adata = mdata.mod["rna"]

    assert adata.X.dtype.kind == "f"
    expected_obs_keys = [
        "AreaUm",
        "Area",
        "Tile",
        "WellLabel",
        "Well",
        "Cell",
        "NuclearAreaUm",
        "NuclearArea",
    ]
    assert all([obs in expected_obs_keys for obs in adata.obs.columns])
    obs_counts = ["Area", "Cell", "NuclearArea"]
    assert all([adata.obs[obs].dtype.kind == "u" for obs in obs_counts])
    obs_areas = ["AreaUm", "NuclearAreaUm"]
    assert all([adata.obs[obs].dtype.kind == "f" for obs in obs_areas])
    obs_categories = ["Tile", "WellLabel", "Well"]
    assert all([adata.obs[obs].dtype.kind == "O" for obs in obs_categories])

    expected_obsm_keys = ["spatial", "spatial_um"]
    assert list(adata.obsm.keys()) == expected_obsm_keys
    assert list(adata.uns.keys()) == expected_obsm_keys
    assert all(adata.obsm[obsm].dtype.kind == "f" for obsm in expected_obsm_keys)


def test_extended_parameters(run_component, tmp_path):
    output = tmp_path / "aviti_ext.h5mu"

    # run component
    run_component(
        [
            "--input",
            input,
            "--modality",
            "mod1",
            "--output",
            str(output),
            "--layer_nuclear_counts",
            "nuclear_counts",
            "--obsm_coordinates",
            "coords",
            "--obsm_cell_paint",
            "cell_paint",
            "--obsm_cell_paint_nuclear",
            "cell_paint_nuclear",
            "--obsm_cell_profiler",
            "cell_profiler",
            "--obsm_unassigned_targets",
            "unassigned_targets",
            "--output_compression",
            "gzip",
        ]
    )

    assert output.is_file(), "output file was not created"

    mdata = mu.read_h5mu(output)
    assert list(mdata.mod.keys()) == ["mod1"]
    adata = mdata.mod["mod1"]

    assert list(adata.layers) == ["nuclear_counts"]
    assert adata.layers["nuclear_counts"].dtype.kind == "f"

    expected_obsm_keys = [
        "cell_paint",
        "cell_paint_nuclear",
        "cell_profiler",
        "coords",
        "coords_um",
        "unassigned_targets",
    ]

    assert list(adata.uns.keys()) == expected_obsm_keys
    assert list(adata.obsm.keys()) == expected_obsm_keys


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
