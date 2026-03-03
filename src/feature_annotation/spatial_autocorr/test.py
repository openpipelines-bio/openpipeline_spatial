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
            "--obsp_neighborhood_graph",
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
        "--obsp_neighborhood_graph",
        "spatial_connectivities",
        "--mode",
        "geary",
        "--n_perms",
        "10",
        "--input_genes",
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
            "--obsp_neighborhood_graph",
            "spatial_connectivities",
            "--mode",
            "moran",
            "--attr",
            "obs",
            "--input_genes",
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


def test_calculate_spatial_autocorrelation_use_all_genes(run_component, tmp_path):
    input_path = meta["resources_dir"] + "/xenium_tiny_neighbors.h5mu"
    temp_input_path = tmp_path / "xenium_tiny_hvg.h5mu"
    output_path_hvg = tmp_path / "output_hvg.h5mu"
    output_path_all = tmp_path / "output_all.h5mu"

    # Create dataset with highly variable genes
    mdata = mu.read_h5mu(input_path)
    adata = mdata.mod["rna"]
    
    # Mark first 5 genes as highly variable
    adata.var["highly_variable"] = False
    genes = list(adata.var_names)
    for g in genes[:5]:
         adata.var.loc[g, "highly_variable"] = True
    
    mdata.write_h5mu(temp_input_path)
    
    # Test 1: Default behavior (should use only HVG)
    print("Running component with default settings (expecting HVG usage)")
    run_component([
        "--input", str(temp_input_path),
        "--output", str(output_path_hvg),
        "--modality", "rna",
        "--obsp_neighborhood_graph", "spatial_connectivities",
        "--mode", "moran",
        "--n_perms", "10"
    ])
    
    assert output_path_hvg.exists()
    mdata_hvg = mu.read_h5mu(output_path_hvg)
    df_hvg = mdata_hvg.mod["rna"].uns["moranI"]
    assert len(df_hvg) == 5, f"Expected 5 genes (HVG), got {len(df_hvg)}"
    
    # Test 2: Override to use all genes
    print("Running component with --use_all_genes (expecting all genes)")
    run_component([
        "--input", str(temp_input_path),
        "--output", str(output_path_all),
        "--modality", "rna",
        "--obsp_neighborhood_graph", "spatial_connectivities",
        "--mode", "moran",
        "--n_perms", "10",
        "--use_all_genes", "true"
    ])
    
    assert output_path_all.exists()
    mdata_all = mu.read_h5mu(output_path_all)
    df_all = mdata_all.mod["rna"].uns["moranI"]
    assert len(df_all) == len(genes), f"Expected all {len(genes)} genes, got {len(df_all)}"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
