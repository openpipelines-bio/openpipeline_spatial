import mudata as mu
import sys
import pytest
import pandas as pd

## VIASH START
meta = {
    "resources_dir": "resources_test/xenium",
    "executable": "./calculate_spatial_autocorrelation",
}
## VIASH END


def test_calculate_spatial_autocorrelation_moran(run_component, tmp_path):
    input_path = meta["resources_dir"] + "/xenium_tiny_neighbors.h5mu"
    output_path = tmp_path / "output_moran.h5mu"

    print(f"Running component with Moran's I on {input_path}")

    run_component(
        [
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--modality",
            "rna",
            "--graph_key",
            "spatial_connectivities",
            "--mode",
            "moran",
            "--n_perms",
            "10",  # Reduce permutations for speed
        ]
    )

    assert output_path.exists(), "Output file not created"

    mdata = mu.read_h5mu(output_path)
    adata = mdata.mod["rna"]

    assert "moranI" in adata.uns, "moranI key missing from .uns"
    df = adata.uns["moranI"]
    assert isinstance(df, pd.DataFrame), "moranI should be a DataFrame"
    assert not df.empty, "moranI DataFrame is empty"

    # Check essential columns
    expected_cols = ["I", "pval_norm", "var_norm"]
    for col in expected_cols:
        assert col in df.columns, f"Missing column {col} in moranI dataframe"

    # Check values range
    assert df["I"].max() <= 1.0, "Moran's I > 1 (theoretical max)"
    assert df["I"].min() >= -1.0, "Moran's I < -1 (theoretical min)"


def test_calculate_spatial_autocorrelation_geary(run_component, tmp_path):
    input_path = meta["resources_dir"] + "/xenium_tiny_neighbors.h5mu"
    output_path = tmp_path / "output_geary.h5mu"

    mdata = mu.read_h5mu(input_path)
    genes = list(mdata.mod["rna"].var_names[:5])  # Test with first 5 genes
    genes_str = ",".join(genes)

    cmd = [
        "--input",
        str(input_path),
        "--output",
        str(output_path),
        "--modality",
        "rna",
        "--graph_key",
        "spatial_connectivities",
        "--mode",
        "geary",
        "--n_perms",
        "10",
        "--genes",
        genes_str,
    ]

    print(
        f"Running component with Geary's C on {input_path} for subset of genes: {genes_str}"
    )
    run_component(cmd)

    assert output_path.exists()

    mdata = mu.read_h5mu(output_path)
    adata = mdata.mod["rna"]

    assert "gearyC" in adata.uns
    df = adata.uns["gearyC"]

    assert len(df) == 5, f"Expected results for 5 genes, got {len(df)}"
    assert all(g in df.index for g in genes), (
        "Not all requested genes are in output index"
    )

    assert "C" in df.columns


def test_calculate_spatial_autocorrelation_obs(run_component, tmp_path):
    input_path = meta["resources_dir"] + "/xenium_tiny_neighbors.h5mu"
    output_path = tmp_path / "output_obs.h5mu"

    features = "total_counts,cell_area,nucleus_area"

    print(
        f"Running component with Moran's I on {input_path} for obs features: {features}"
    )

    run_component(
        [
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--modality",
            "rna",
            "--graph_key",
            "spatial_connectivities",
            "--mode",
            "moran",
            "--attr",
            "obs",
            "--genes",
            features,
            "--n_perms",
            "10",
        ]
    )

    assert output_path.exists()

    mdata = mu.read_h5mu(output_path)
    adata = mdata.mod["rna"]

    assert "moranI" in adata.uns
    df = adata.uns["moranI"]

    assert len(df) == 3, f"Expected results for 3 features, got {len(df)}"
    # When attr='obs', the index contains the column names
    for feature in features.split(","):
        assert feature in df.index, f"Feature {feature} missing from index: {df.index}"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
