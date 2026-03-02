import pytest
import mudata as mu
import sys

## VIASH START
meta = {
    "executable": "./target/executable/cluster/leiden_spatial/leiden_spatial",
    "resources_dir": "resources_test/xenium/",
}
## VIASH END

# This file has spatial_connectivities, X_pca, and connectivities pre-computed.
input_file = f"{meta['resources_dir']}/xenium_tiny.spatial_expression_neighbors.h5mu"


def test_default_execution(run_component, tmp_path):
    output = tmp_path / "leiden_spatial_output.h5mu"

    run_component(
        [
            "--input", input_file,
            "--output", str(output),
            "--output_compression", "gzip",
        ]
    )

    assert output.is_file(), "Output file was not created."
    mdata = mu.read_h5mu(output)
    assert "rna" in mdata.mod, "Expected 'rna' modality in output."
    adata = mdata.mod["rna"]

    assert "leiden_spatial" in adata.obs.columns, (
        "Expected 'leiden_spatial' column in .obs."
    )
    assert str(adata.obs["leiden_spatial"].dtype) == "category", (
        "Expected 'leiden_spatial' column to be categorical."
    )
    n_domains = adata.obs["leiden_spatial"].nunique()
    assert n_domains >= 1, "Expected at least one spatial domain."


def test_custom_obs_label(run_component, tmp_path):
    output = tmp_path / "custom_label_output.h5mu"

    run_component(
        [
            "--input", input_file,
            "--output", str(output),
            "--obs_label", "my_spatial_domains",
        ]
    )

    assert output.is_file(), "Output file was not created."
    mdata = mu.read_h5mu(output)
    adata = mdata.mod["rna"]
    assert "my_spatial_domains" in adata.obs.columns, (
        "Expected custom obs label 'my_spatial_domains' in .obs."
    )


def test_spatial_only_alpha(run_component, tmp_path):
    """Test that alpha=1.0 (purely spatial clustering) runs without errors."""
    output = tmp_path / "alpha_one_output.h5mu"

    run_component(
        [
            "--input", input_file,
            "--output", str(output),
            "--alpha", "1.0",
        ]
    )

    assert output.is_file(), "Output file was not created."
    mdata = mu.read_h5mu(output)
    assert "leiden_spatial" in mdata.mod["rna"].obs.columns


def test_expression_only_alpha(run_component, tmp_path):
    """Test that alpha=0.0 (purely expression-based clustering) runs without errors."""
    output = tmp_path / "alpha_zero_output.h5mu"

    run_component(
        [
            "--input", input_file,
            "--output", str(output),
            "--alpha", "0.0",
        ]
    )

    assert output.is_file(), "Output file was not created."
    mdata = mu.read_h5mu(output)
    assert "leiden_spatial" in mdata.mod["rna"].obs.columns


def test_spatial_connectivities_key_not_found(run_component, tmp_path):
    """Test that a clear error is raised when spatial connectivities are missing."""
    output = tmp_path / "error_output.h5mu"

    with pytest.raises(Exception):
        run_component(
            [
                "--input", input_file,
                "--output", str(output),
                "--input_obsp_spatial_connectivities", "nonexistent_key",
            ]
        )


def test_expression_connectivities_key_not_found(run_component, tmp_path):
    """Test that a clear error is raised when the expression connectivities key is missing."""
    # Raw data has no expression connectivities (.obsp only has spatial_connectivities)
    raw_mdata = mu.read_h5mu(input_file)
    no_expr = tmp_path / "no_expr_neighbors.h5mu"
    adata = raw_mdata.mod["rna"]
    del adata.obsp["connectivities"]
    del adata.obsp["distances"]
    mu.MuData({"rna": adata}).write_h5mu(no_expr)

    output = tmp_path / "error_output.h5mu"
    with pytest.raises(Exception):
        run_component(
            [
                "--input", str(no_expr),
                "--output", str(output),
                "--input_obsp_expression_connectivities", "connectivities",
            ]
        )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
