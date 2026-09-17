import sys
import numpy as np
import mudata as mu
import scanpy as sc
import pytest

## VIASH START
meta = {
    "executable": "./target/executable/banksy/banksy/banksy",
}
## VIASH END

input_path_xenium = f"{meta['resources_dir']}/xenium_tiny.h5mu"


@pytest.fixture
def input_normalized(tmp_path):
    mdata = mu.read_h5mu(input_path_xenium)
    adata = mdata.mod["rna"]
    sc.pp.filter_genes(adata, min_cells=1)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    output = tmp_path / "xenium_tiny_normalized.h5mu"
    mdata.write_h5mu(output)
    return str(output)


def test_simple_execution(run_component, tmp_path, input_normalized):
    output = tmp_path / "banksy_xenium.h5mu"

    run_component(
        [
            "--input",
            input_normalized,
            "--k_geom",
            "5",
            "--pca_dims",
            "10",
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

    assert "banksy_embedding" in adata.obsm, "Expected BANKSY embedding in .obsm"
    assert adata.obsm["banksy_embedding"].shape[0] == adata.n_obs, (
        "Embedding should have one row per cell"
    )
    assert adata.obsm["banksy_embedding"].dtype.kind == "f", (
        "Expected embedding to be float type"
    )

    assert "banksy_cluster" in adata.obs, "Expected cluster labels in .obs"
    assert adata.obs["banksy_cluster"].notna().all(), (
        "Expected no missing cluster labels"
    )


def test_custom_output_keys(run_component, tmp_path, input_normalized):
    output = tmp_path / "banksy_xenium_custom_keys.h5mu"

    run_component(
        [
            "--input",
            input_normalized,
            "--k_geom",
            "5",
            "--pca_dims",
            "10",
            "--output_obsm_embedding",
            "my_embedding",
            "--output_obs_cluster",
            "my_cluster",
            "--output",
            str(output),
        ]
    )

    assert output.is_file(), "output file was not created"
    mdata = mu.read_h5mu(output)
    adata = mdata.mod["rna"]
    assert "my_embedding" in adata.obsm, "Expected embedding under custom obsm key"
    assert "my_cluster" in adata.obs, "Expected cluster labels under custom obs key"


def test_lambda_changes_embedding(run_component, tmp_path, input_normalized):
    # lambda=0 ignores the spatial neighbourhood entirely (plain cell typing);
    # lambda=0.8 (domain segmentation) should give a materially different
    # embedding, since it's dominated by the neighbourhood expression blocks.
    embeddings = {}
    for lam in ["0.0", "0.8"]:
        output = tmp_path / f"banksy_lambda_{lam}.h5mu"
        run_component(
            [
                "--input",
                input_normalized,
                "--k_geom",
                "5",
                "--pca_dims",
                "10",
                "--lambda_param",
                lam,
                "--output",
                str(output),
            ]
        )
        mdata = mu.read_h5mu(output)
        embeddings[lam] = mdata.mod["rna"].obsm["banksy_embedding"]

    assert embeddings["0.0"].shape == embeddings["0.8"].shape
    assert not np.allclose(embeddings["0.0"], embeddings["0.8"]), (
        "Expected lambda=0.0 and lambda=0.8 to produce different embeddings"
    )


def test_var_input(run_component, tmp_path, input_normalized):
    mdata = mu.read_h5mu(input_normalized)
    adata = mdata.mod["rna"]
    np.random.seed(42)
    n_select = min(500, adata.n_vars)
    selected = np.random.choice(adata.n_vars, size=n_select, replace=False)
    mask = np.zeros(adata.n_vars, dtype=bool)
    mask[selected] = True
    adata.var["filter_with_hvg"] = mask
    input_with_hvg = tmp_path / "xenium_with_hvg_filter.h5mu"
    mdata.write_h5mu(input_with_hvg)

    output = tmp_path / "banksy_xenium_var_input.h5mu"
    run_component(
        [
            "--input",
            str(input_with_hvg),
            "--var_input",
            "filter_with_hvg",
            "--k_geom",
            "5",
            "--pca_dims",
            "10",
            "--output",
            str(output),
        ]
    )

    assert output.is_file(), "output file was not created"
    mdata_out = mu.read_h5mu(output)
    adata_out = mdata_out.mod["rna"]
    # output retains the full (unfiltered) gene space
    assert adata_out.n_vars == adata.n_vars, (
        "var_input should only affect the genes used for computation, not the output gene space"
    )
    assert "banksy_embedding" in adata_out.obsm


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
