import sys
import numpy as np
import pytest
from mudata import read_h5mu

##VIASH START
par = {
    "input": "scvi_cellcharter/output.h5mu",
    "input_obs_library_key": "sample_id",
    "n_dimensions_latent_space": 10,
    "n_aggregated_blocks": 4,
    "expected_n_clusters": None,
    "expect_stability": True,
}

meta = {"resources_dir": "resources_test"}
##VIASH END


@pytest.fixture(scope="module")
def adata():
    mdata = read_h5mu(par["input"])
    assert "rna" in mdata.mod, f"Expected modality 'rna', found: {list(mdata.mod)}."
    return mdata.mod["rna"]


def test_embeddings(adata):
    assert "X_scvi" in adata.obsm, "Expected the scVI embedding in .obsm['X_scvi']."
    assert adata.obsm["X_scvi"].shape == (
        adata.n_obs,
        par["n_dimensions_latent_space"],
    )

    assert "X_cellcharter" in adata.obsm, (
        "Expected the aggregated features in .obsm['X_cellcharter']."
    )
    aggregated = np.asarray(adata.obsm["X_cellcharter"])
    assert aggregated.shape == (
        adata.n_obs,
        par["n_dimensions_latent_space"] * par["n_aggregated_blocks"],
    )
    np.testing.assert_allclose(
        aggregated[:, : par["n_dimensions_latent_space"]],
        adata.obsm["X_scvi"],
        rtol=1e-5,
        atol=1e-5,
        err_msg="Expected the first block of the aggregated features to be the scVI embedding.",
    )


def test_spatial_graph(adata):
    assert "spatial_connectivities" in adata.obsp, (
        "Expected the spatial neighborhood graph in .obsp['spatial_connectivities']."
    )
    connectivities = adata.obsp["spatial_connectivities"]
    assert connectivities.nnz > 0, "Expected edges in the spatial neighborhood graph."

    libraries = adata.obs[par["input_obs_library_key"]].to_numpy()
    rows, cols = connectivities.nonzero()
    assert (libraries[rows] == libraries[cols]).all(), (
        "Expected no spatial edges between observations of different libraries."
    )


def test_clusters(adata):
    assert "cellcharter_cluster" in adata.obs, (
        "Expected the spatial domains in .obs['cellcharter_cluster']."
    )
    labels = adata.obs["cellcharter_cluster"]
    assert labels.dtype.name == "category", "Expected categorical cluster labels."
    assert not labels.isna().any(), "Expected every observation to have a cluster."
    if par["expected_n_clusters"] is not None:
        assert len(labels.cat.categories) == par["expected_n_clusters"]


def test_stability(adata):
    if not par["expect_stability"]:
        assert "cellcharter_autok" not in adata.uns
        return

    assert "cellcharter_autok" in adata.uns, (
        "Expected the stability analysis results in .uns['cellcharter_autok']."
    )
    stability = adata.uns["cellcharter_autok"]
    for key in ["n_clusters", "stability", "stability_mean", "best_k", "peaks", "k"]:
        assert key in stability, f"Expected '{key}' in the stability analysis results."
    assert stability["k"] in stability["n_clusters"]
    assert len(adata.obs["cellcharter_cluster"].cat.categories) == stability["k"]


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "--import-mode=importlib"]))
