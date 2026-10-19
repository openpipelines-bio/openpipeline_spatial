import pytest
import mudata as mu
import numpy as np
import re
import subprocess
import sys

## VIASH START
meta = {
    "executable": "./target/executable/neighbors/spatial_neighborhood_graph/spatial_neighborhood_graph",
}
## VIASH END

input_xenium = f"{meta['resources_dir']}/xenium_tiny.h5mu"
input_cosmx = f"{meta['resources_dir']}/Lung5_Rep2_tiny.h5mu"


def test_simple_execution_xenium(run_component, tmp_path):
    output = tmp_path / "nc_xenium.h5mu"

    # run component
    run_component(
        [
            "--input",
            input_xenium,
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

    expected_obsp_keys = ["spatial_connectivities", "spatial_distances"]
    assert all([obsp in adata.obsp.keys() for obsp in expected_obsp_keys]), (
        "Not all expected obsp keys found"
    )
    assert all(adata.obsp[obsp].dtype.kind == "f" for obsp in expected_obsp_keys), (
        "Expected obsp matrices to be float type"
    )


def test_multiple_libraries(run_component, tmp_path):
    # split the observations into two libraries, stored as a non-categorical column
    mdata = mu.read_h5mu(input_xenium)
    n_obs = mdata.mod["rna"].n_obs
    libraries = np.where(np.arange(n_obs) % 2 == 0, "sample_a", "sample_b")
    mdata.mod["rna"].obs["sample_id"] = libraries.astype(object)
    input_multi = tmp_path / "xenium_multi_library.h5mu"
    mdata.write_h5mu(input_multi)

    output = tmp_path / "nc_xenium_multi_library.h5mu"

    run_component(
        [
            "--input",
            str(input_multi),
            "--input_obs_library_key",
            "sample_id",
            "--output",
            str(output),
        ]
    )

    assert output.is_file(), "output file was not created"
    adata = mu.read_h5mu(output).mod["rna"]
    assert adata.obs["sample_id"].dtype.name == "category", (
        "Expected library key to be stored as categorical"
    )

    connectivities = adata.obsp["spatial_connectivities"]
    is_a = libraries == "sample_a"
    is_b = libraries == "sample_b"
    assert connectivities[is_a][:, is_a].nnz > 0, "Expected edges within sample_a"
    assert connectivities[is_b][:, is_b].nnz > 0, "Expected edges within sample_b"
    assert connectivities[is_a][:, is_b].nnz == 0, (
        "Expected no edges between observations of different libraries"
    )


def test_missing_library_key(run_component, tmp_path):
    output = tmp_path / "nc_xenium_missing_key.h5mu"

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_xenium,
                "--input_obs_library_key",
                "does_not_exist",
                "--output",
                str(output),
            ]
        )
    assert re.search(
        r"--input_obs_library_key 'does_not_exist' not found in .obs of modality 'rna'.",
        err.value.stdout.decode("utf-8"),
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
