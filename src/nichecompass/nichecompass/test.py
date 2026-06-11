import pytest
import mudata as mu
import numpy as np
import sys

## VIASH START
meta = {
    "executable": "./target/executable/nichecompass/nichecompass/nichecompass",
}
## VIASH END

input_path_xenium = f"{meta['resources_dir']}/xenium_tiny.h5mu"
input_path_cosmx = f"{meta['resources_dir']}/Lung5_Rep2_tiny.h5mu"
gp_mask = f"{meta['resources_dir']}/prior_knowledge_gp_mask.json"


@pytest.fixture
def input_xenium_with_hvg_filter(tmp_path):
    mdata = mu.read_h5mu(input_path_xenium)
    adata = mdata.mod["rna"]
    np.random.seed(42)
    n_select = min(1000, adata.n_vars)
    selected = np.random.choice(adata.n_vars, size=n_select, replace=False)
    mask = np.zeros(adata.n_vars, dtype=bool)
    mask[selected] = True
    adata.var["filter_with_hvg"] = mask
    output = tmp_path / "xenium_with_hvg_filter.h5mu"
    mdata.write_h5mu(output)
    return str(output)


def test_simple_execution_xenium(run_component, tmp_path, input_xenium_with_hvg_filter):
    output = tmp_path / "nc_xenium.h5mu"

    # run component
    run_component(
        [
            "--input",
            input_xenium_with_hvg_filter,
            "--input_gp_mask",
            gp_mask,
            "--var_input",
            "filter_with_hvg",
            "--n_epochs",
            "1",
            "--n_epochs_all_gps",
            "0",
            "--n_epochs_no_edge_recon",
            "0",
            "--n_epochs_no_cat_covariates_contrastive",
            "0",
            "--output",
            str(output),
            "--output_model",
            "test_model",
            "--output_compression",
            "gzip",
        ]
    )

    assert output.is_file(), "output file was not created"
    mdata = mu.read_h5mu(output)
    assert list(mdata.mod.keys()) == ["rna"], "Expected modality rna"
    adata = mdata.mod["rna"]

    expected_uns_keys = [
        "nichecompass_sources_categories_label_encoder",
        "nichecompass_targets_categories_label_encoder",
        "nichecompass_source_genes_idx",
        "nichecompass_target_genes_idx",
        "nichecompass_genes_idx",
        "nichecompass_gp_names",
        "nichecompass_active_gp_names",
    ]
    assert all([uns in adata.uns.keys() for uns in expected_uns_keys]), (
        f"Expected uns keys: {expected_uns_keys}, found: {list(adata.uns.keys())}"
    )
    assert len(adata.uns["nichecompass_gp_names"]) > len(
        adata.uns["nichecompass_active_gp_names"]
    ), "Expected less active GP names than total GP names"

    expected_obsm_keys = ["nichecompass_latent"]
    assert all([obsm in adata.obsm.keys() for obsm in expected_obsm_keys]), (
        "Not all expected obsm keys found"
    )
    assert all(adata.obsm[obsm].dtype.kind == "f" for obsm in expected_obsm_keys), (
        "Expected obsm matrices to be float type"
    )

    expected_varm_keys = [
        "nichecompass_gp_sources",
        "nichecompass_gp_targets",
        "nichecompass_gp_sources_categories",
        "nichecompass_gp_targets_categories",
    ]
    assert all([varm in adata.varm.keys() for varm in expected_varm_keys]), (
        "Not all expected varm keys found"
    )
    assert (
        adata.varm["nichecompass_gp_targets"].shape
        == adata.varm["nichecompass_gp_sources"].shape
    ), "Expected GP targets and sources varm to have same shape"


def test_no_var_input(run_component, tmp_path, input_xenium_with_hvg_filter):
    output = tmp_path / "nc_xenium_no_var_input.h5mu"

    # run component without var_input (all genes are used)
    run_component(
        [
            "--input",
            input_xenium_with_hvg_filter,
            "--input_gp_mask",
            gp_mask,
            "--var_input",
            "",
            "--n_epochs",
            "1",
            "--n_epochs_all_gps",
            "0",
            "--n_epochs_no_edge_recon",
            "0",
            "--n_epochs_no_cat_covariates_contrastive",
            "0",
            "--output",
            str(output),
            "--output_model",
            "test_model_no_var_input",
            "--output_compression",
            "gzip",
        ]
    )

    assert output.is_file(), "output file was not created"
    mdata = mu.read_h5mu(output)
    assert list(mdata.mod.keys()) == ["rna"], "Expected modality rna"
    adata = mdata.mod["rna"]

    expected_uns_keys = [
        "nichecompass_gp_names",
        "nichecompass_active_gp_names",
    ]
    assert all([uns in adata.uns.keys() for uns in expected_uns_keys]), (
        f"Expected uns keys: {expected_uns_keys}, found: {list(adata.uns.keys())}"
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
