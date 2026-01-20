import pytest
import mudata as mu
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


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
