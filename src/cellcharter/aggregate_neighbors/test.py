import re
import subprocess
import sys
import mudata as mu
import numpy as np
import pytest
import scipy.sparse as sp

## VIASH START
meta = {
    "executable": "./target/executable/cellcharter/aggregate_neighbors/aggregate_neighbors",
    "resources_dir": "resources_test/xenium/",
}
## VIASH END

input_xenium = f"{meta['resources_dir']}/xenium_tiny.spatial_expression_neighbors.h5mu"


def to_dense(x):
    return x.toarray() if sp.issparse(x) else np.asarray(x)


@pytest.fixture
def input_multi_sample(tmp_path):
    # split the observations into two samples and add a layer to aggregate
    mdata = mu.read_h5mu(input_xenium)
    adata = mdata.mod["rna"]
    samples = np.where(np.arange(adata.n_obs) % 2 == 0, "sample_a", "sample_b")
    adata.obs["sample_id"] = samples.astype(object)
    adata.layers["log_normalized"] = adata.X.copy()
    output = tmp_path / "xenium_multi_sample.h5mu"
    mdata.write_h5mu(output)
    return output


def test_simple_execution(run_component, tmp_path):
    output = tmp_path / "output.h5mu"

    run_component(
        [
            "--input",
            input_xenium,
            "--input_obsm_features",
            "X_pca",
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

    assert "X_cellcharter" in adata.obsm, "Expected X_cellcharter in .obsm"
    aggregated = to_dense(adata.obsm["X_cellcharter"])
    x_pca = adata.obsm["X_pca"]
    n_pcs = x_pca.shape[1]

    # default n_layers=3 with mean aggregation: cell | 1-hop | 2-hop | 3-hop
    assert aggregated.shape == (adata.n_obs, n_pcs * 4)
    np.testing.assert_allclose(aggregated[:, :n_pcs], x_pca, rtol=1e-5, atol=1e-5)
    assert np.isfinite(aggregated).all(), "Expected no NaN or infinite values"


def test_multiple_layers_and_aggregations(run_component, tmp_path):
    output = tmp_path / "output.h5mu"

    run_component(
        [
            "--input",
            input_xenium,
            "--input_obsm_features",
            "X_pca",
            "--n_layers",
            "0;1;3",
            "--aggregations",
            "mean;var",
            "--output",
            str(output),
            "--output_obsm_aggregated",
            "X_aggregated",
        ]
    )

    adata = mu.read_h5mu(output).mod["rna"]
    assert "X_aggregated" in adata.obsm, "Expected X_aggregated in .obsm"
    assert "X_cellcharter" not in adata.obsm
    n_pcs = adata.obsm["X_pca"].shape[1]

    # cell | 1-hop mean | 1-hop var | 3-hop mean | 3-hop var
    aggregated = to_dense(adata.obsm["X_aggregated"])
    assert aggregated.shape == (adata.n_obs, n_pcs * 5)
    variances = np.hstack(
        [aggregated[:, 2 * n_pcs : 3 * n_pcs], aggregated[:, 4 * n_pcs :]]
    )
    assert (variances > -1e-4).all(), "Expected non-negative variances"


def test_sample_key_and_layer(run_component, tmp_path, input_multi_sample):
    output = tmp_path / "output.h5mu"

    run_component(
        [
            "--input",
            str(input_multi_sample),
            "--layer",
            "log_normalized",
            "--input_obs_sample_key",
            "sample_id",
            "--n_layers",
            "1",
            "--output",
            str(output),
        ]
    )

    adata = mu.read_h5mu(output).mod["rna"]
    features = to_dense(adata.layers["log_normalized"])
    aggregated = to_dense(adata.obsm["X_cellcharter"])
    assert aggregated.shape == (adata.n_obs, adata.n_vars * 2)

    # The 1-hop mean only takes the neighbors from the same sample into account
    connectivities = adata.obsp["spatial_connectivities"]
    samples = adata.obs["sample_id"].to_numpy()
    assert any(
        samples[i] != samples[j] for i, j in zip(*connectivities.nonzero()) if i != j
    ), "Expected the spatial graph to contain edges between samples"

    expected = np.zeros_like(features)
    for i in range(adata.n_obs):
        neighbors = [
            j
            for j in connectivities[i].nonzero()[1]
            if j != i and samples[j] == samples[i]
        ]
        if neighbors:
            expected[i] = features[neighbors].mean(axis=0)
    np.testing.assert_allclose(
        aggregated[:, adata.n_vars :], expected, rtol=1e-5, atol=1e-5
    )


def test_missing_connectivities(run_component, tmp_path):
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_xenium,
                "--input_obsp_spatial_connectivities",
                "does_not_exist",
                "--output",
                str(tmp_path / "output.h5mu"),
            ]
        )
    assert re.search(
        r"--input_obsp_spatial_connectivities 'does_not_exist' not found in .obsp of modality 'rna'",
        err.value.stdout.decode("utf-8"),
    )


def test_missing_sample_key(run_component, tmp_path):
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_xenium,
                "--input_obs_sample_key",
                "does_not_exist",
                "--output",
                str(tmp_path / "output.h5mu"),
            ]
        )
    assert re.search(
        r"--input_obs_sample_key 'does_not_exist' not found in .obs of modality 'rna'",
        err.value.stdout.decode("utf-8"),
    )


def test_obsm_and_layer_mutually_exclusive(run_component, tmp_path):
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_xenium,
                "--input_obsm_features",
                "X_pca",
                "--layer",
                "log_normalized",
                "--output",
                str(tmp_path / "output.h5mu"),
            ]
        )
    assert re.search(
        r"Only one of --input_obsm_features and --layer can be provided",
        err.value.stdout.decode("utf-8"),
    )


def test_n_layers_without_zero(run_component, tmp_path):
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_xenium,
                "--input_obsm_features",
                "X_pca",
                "--n_layers",
                "1;2",
                "--output",
                str(tmp_path / "output.h5mu"),
            ]
        )
    assert re.search(
        r"When providing multiple values for --n_layers, 0 \(the cells' own features\) must be included",
        err.value.stdout.decode("utf-8"),
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
