import pytest
import sys
import mudata as mu
import subprocess


def test_simple_execution(run_component, tmp_path):
    output = tmp_path / "cosmx_tiny.h5mu"

    run_component(
        [
            "--input",
            meta["resources_dir"] + "/Lung5_Rep2_tiny",
            "--output",
            output,
        ]
    )
    assert output.is_file(), "output file was not created"

    mdata = mu.read_h5mu(output)
    assert list(mdata.mod.keys()) == ["rna"], "Expected modality rna"

    adata = mdata.mod["rna"]

    assert list(adata.obs.keys()) == [
        "fov",
        "Area",
        "AspectRatio",
        "CenterX_global_px",
        "CenterY_global_px",
        "Width",
        "Height",
        "Mean.MembraneStain",
        "Max.MembraneStain",
        "Mean.PanCK",
        "Max.PanCK",
        "Mean.CD45",
        "Max.CD45",
        "Mean.CD3",
        "Max.CD3",
        "Mean.DAPI",
        "Max.DAPI",
        "cell_ID",
    ]

    assert list(adata.uns.keys()) == ["spatial"]
    assert list(adata.obsm.keys()) == ["spatial", "spatial_fov"]

    assert adata.obsm["spatial"].dtype == "int"
    assert adata.obsm["spatial_fov"].dtype == "float"


def test_compressed_input(run_component, tmp_path):
    output = tmp_path / "cosmx_tiny.h5mu"
    zipped_input = tmp_path / "Lung5_Rep2_tiny.zip"

    subprocess.run(
        ["zip", "-r", str(zipped_input), "Lung5_Rep2_tiny"],
        cwd=meta["resources_dir"],
        check=True,
    )

    run_component(
        [
            "--input",
            zipped_input,
            "--output",
            output,
        ]
    )
    assert output.is_file(), "output file was not created"

    mdata = mu.read_h5mu(output)
    assert list(mdata.mod.keys()) == ["rna"], "Expected modality rna"

    adata = mdata.mod["rna"]

    assert list(adata.obs.keys()) == [
        "fov",
        "Area",
        "AspectRatio",
        "CenterX_global_px",
        "CenterY_global_px",
        "Width",
        "Height",
        "Mean.MembraneStain",
        "Max.MembraneStain",
        "Mean.PanCK",
        "Max.PanCK",
        "Mean.CD45",
        "Max.CD45",
        "Mean.CD3",
        "Max.CD3",
        "Mean.DAPI",
        "Max.DAPI",
        "cell_ID",
    ]

    assert list(adata.uns.keys()) == ["spatial"]
    assert list(adata.obsm.keys()) == ["spatial", "spatial_fov"]

    assert adata.obsm["spatial"].dtype == "int"
    assert adata.obsm["spatial_fov"].dtype == "float"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
