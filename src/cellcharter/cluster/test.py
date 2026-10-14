import re
import subprocess
import sys
import anndata as ad
import mudata as mu
import numpy as np
import pytest
from sklearn.metrics import adjusted_rand_score

## VIASH START
meta = {
    "executable": "./target/executable/cellcharter/cluster/cluster",
    "resources_dir": "resources_test/xenium/",
}
## VIASH END

input_xenium = f"{meta['resources_dir']}/xenium_tiny.qc.all_neighbors.pca.h5mu"

N_BLOBS = 3


@pytest.fixture
def input_blobs(tmp_path):
    # Three well-separated Gaussian blobs as neighborhood-aggregated features
    rng = np.random.default_rng(0)
    n_per_blob, n_features = 200, 8
    centers = rng.normal(scale=10, size=(N_BLOBS, n_features))
    features = np.vstack(
        [center + rng.normal(size=(n_per_blob, n_features)) for center in centers]
    ).astype(np.float32)
    truth = np.repeat(np.arange(N_BLOBS), n_per_blob)

    adata = ad.AnnData(
        X=np.zeros((features.shape[0], 1), dtype=np.float32),
        obsm={"X_cellcharter": features},
    )
    adata.obs["truth"] = truth.astype(str)
    output = tmp_path / "blobs.h5mu"
    mu.MuData({"rna": adata}).write_h5mu(output)
    return output


def test_fixed_n_clusters(run_component, tmp_path, input_blobs):
    output = tmp_path / "output.h5mu"

    run_component(
        [
            "--input",
            str(input_blobs),
            "--n_clusters",
            str(N_BLOBS),
            "--output",
            str(output),
            "--output_compression",
            "gzip",
        ]
    )

    assert output.is_file(), "output file was not created"
    adata = mu.read_h5mu(output).mod["rna"]
    assert "cellcharter_cluster" in adata.obs, "Expected cluster labels in .obs"
    labels = adata.obs["cellcharter_cluster"]
    assert labels.dtype.name == "category", "Expected categorical cluster labels"
    assert list(labels.cat.categories) == [str(i) for i in range(N_BLOBS)]
    assert adjusted_rand_score(adata.obs["truth"], labels) > 0.99, (
        "Expected the clusters to recover the blobs"
    )
    assert "cellcharter_autok" not in adata.uns, (
        "Expected no stability analysis without --n_clusters_range"
    )


def test_autok(run_component, tmp_path, input_blobs):
    output = tmp_path / "output.h5mu"

    run_component(
        [
            "--input",
            str(input_blobs),
            "--n_clusters_range",
            "2;5",
            "--max_runs",
            "3",
            "--covariance_type",
            "diag",
            "--output",
            str(output),
            "--output_obs_cluster",
            "domains",
            "--output_uns_stability",
            "autok",
        ]
    )

    adata = mu.read_h5mu(output).mod["rna"]
    stability = adata.uns["autok"]
    np.testing.assert_array_equal(stability["n_clusters"], [2, 3, 4, 5])
    assert stability["stability"].shape[0] == 4
    np.testing.assert_allclose(
        stability["stability_mean"], stability["stability"].mean(axis=1)
    )
    assert stability["best_k"] == N_BLOBS, "Expected the most stable K to be 3"
    assert stability["k"] == N_BLOBS

    labels = adata.obs["domains"]
    assert labels.nunique() == N_BLOBS
    assert adjusted_rand_score(adata.obs["truth"], labels) > 0.99


def test_autok_with_n_clusters(run_component, tmp_path, input_blobs):
    output = tmp_path / "output.h5mu"

    run_component(
        [
            "--input",
            str(input_blobs),
            "--n_clusters_range",
            "2;5",
            "--n_clusters",
            "4",
            "--max_runs",
            "2",
            "--covariance_type",
            "diag",
            "--output",
            str(output),
        ]
    )

    adata = mu.read_h5mu(output).mod["rna"]
    stability = adata.uns["cellcharter_autok"]
    assert stability["k"] == 4, "Expected the provided --n_clusters to be used"
    assert list(adata.obs["cellcharter_cluster"].cat.categories) == [
        "0",
        "1",
        "2",
        "3",
    ]


def test_reproducibility(run_component, tmp_path):
    def run(output):
        run_component(
            [
                "--input",
                input_xenium,
                "--input_obsm_features",
                "X_pca",
                "--n_clusters",
                "3",
                "--covariance_type",
                "diag",
                "--init_strategy",
                "kmeans++",
                "--random_state",
                "42",
                "--output",
                str(output),
            ]
        )
        return mu.read_h5mu(output).mod["rna"].obs["cellcharter_cluster"]

    labels_1 = run(tmp_path / "output_1.h5mu")
    labels_2 = run(tmp_path / "output_2.h5mu")
    assert len(labels_1) == mu.read_h5mu(input_xenium).mod["rna"].n_obs
    assert labels_1.equals(labels_2), "Expected identical labels with the same seed"


def test_missing_n_clusters(run_component, tmp_path, input_blobs):
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            ["--input", str(input_blobs), "--output", str(tmp_path / "output.h5mu")]
        )
    assert re.search(
        r"At least one of --n_clusters or --n_clusters_range must be provided",
        err.value.stdout.decode("utf-8"),
    )


@pytest.mark.parametrize(
    "extra_args, message",
    [
        (["--n_clusters_range", "2;5;8"], r"must contain exactly two values"),
        (["--n_clusters_range", "5;2"], r"must contain exactly two values"),
        (
            ["--n_clusters_range", "2;5", "--n_clusters", "6"],
            r"--n_clusters \(6\) must be within --n_clusters_range \[2, 5\]",
        ),
        (
            ["--n_clusters", "3", "--input_obsm_features", "does_not_exist"],
            r"--input_obsm_features 'does_not_exist' not found in .obsm of modality 'rna'",
        ),
    ],
)
def test_invalid_arguments(run_component, tmp_path, input_blobs, extra_args, message):
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                str(input_blobs),
                "--output",
                str(tmp_path / "output.h5mu"),
            ]
            + extra_args
        )
    assert re.search(message, err.value.stdout.decode("utf-8"))


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
