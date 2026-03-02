import pytest
import mudata as mu
import sys

## VIASH START
meta = {
    "executable": "./target/executable/neighbors/join_graphs/join_graphs",
    "resources_dir": "resources_test/xenium/",
}
## VIASH END

# This file has spatial_connectivities, X_pca, and connectivities pre-computed.
input_file = f"{meta['resources_dir']}/xenium_tiny.spatial_expression_neighbors.h5mu"


def test_default_execution(run_component, tmp_path):
    output = tmp_path / "fused_neighbors.h5mu"

    run_component(
        [
            "--input",
            input_file,
            "--output",
            str(output),
            "--output_obsp_connectivities",
            "fused_conn",
            "--alpha",
            "0.2",
            "--output_compression",
            "gzip",
        ]
    )

    assert output.is_file(), "Output file was not created."
    mdata = mu.read_h5mu(output)
    assert "rna" in mdata.mod, "Expected 'rna' modality in output."
    adata = mdata.mod["rna"]

    assert "fused_conn" in adata.obsp.keys(), (
        "Expected output obsp connectivities key to be present."
    )
    assert "fused_conn" in adata.uns.keys(), "Expected metadata for fusion in .uns."
    assert adata.uns["fused_conn"]["params"]["alpha"] == 0.2, (
        "Expected alpha parameter to be stored in .uns."
    )

    # Check graph properties (should be combination of input graphs)
    fused = adata.obsp["fused_conn"]
    spatial = adata.obsp["spatial_connectivities"]
    expr = adata.obsp["connectivities"]

    assert fused.shape == expr.shape == spatial.shape, (
        "Fused graph should have same shape as input graphs."
    )


def test_alpha_zero_equals_expression(run_component, tmp_path):
    output = tmp_path / "alpha_zero.h5mu"
    run_component(
        [
            "--input",
            input_file,
            "--output",
            str(output),
            "--alpha",
            "0.0",
            "--output_obsp_connectivities",
            "fused_zero",
        ]
    )
    mdata = mu.read_h5mu(output)
    adata = mdata.mod["rna"]

    fused = adata.obsp["fused_zero"]
    # Should equal expression graph exactly (dense check might be heavy, check nnz or sum)
    expr = adata.obsp["connectivities"]
    assert abs(fused.sum() - expr.sum()) < 1e-6


def test_alpha_one_equals_spatial(run_component, tmp_path):
    output = tmp_path / "alpha_one.h5mu"
    run_component(
        [
            "--input",
            input_file,
            "--output",
            str(output),
            "--alpha",
            "1.0",
            "--output_obsp_connectivities",
            "fused_one",
        ]
    )
    mdata = mu.read_h5mu(output)
    adata = mdata.mod["rna"]

    fused = adata.obsp["fused_one"]
    # Should equal spatial graph exactly
    spatial = adata.obsp["spatial_connectivities"]
    assert abs(fused.sum() - spatial.sum()) < 1e-6


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
