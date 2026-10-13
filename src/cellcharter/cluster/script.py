import sys
import anndata as ad
import cellcharter as cc
import mudata as mu
import numpy as np
import scipy.sparse as sp
import torch
from pytorch_lightning import seed_everything

## VIASH START
par = {
    # Inputs
    "input": "resources_test/xenium/xenium_tiny.qc.all_neighbors.pca.h5mu",
    "modality": "rna",
    "input_obsm_features": "X_pca",
    # Number of clusters
    "n_clusters": None,
    "n_clusters_range": [2, 4],
    "max_runs": 3,
    "stability_convergence_tolerance": 0.01,
    # Gaussian Mixture Model
    "covariance_type": "diag",
    "init_strategy": "kmeans",
    "convergence_tolerance": 0.001,
    "covariance_regularization": 1e-6,
    "max_epochs": 100,
    "batch_size": None,
    "accelerator": "auto",
    "random_state": 0,
    # Outputs
    "output": "output.h5mu",
    "output_obs_cluster": "cellcharter_cluster",
    "output_uns_stability": "cellcharter_autok",
    "output_compression": None,
}
meta = {"resources_dir": "src/utils/", "cpus": None}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

## Validate the number of clusters
n_clusters = par["n_clusters"]
n_clusters_range = par["n_clusters_range"]
if n_clusters is None and not n_clusters_range:
    raise ValueError(
        "At least one of --n_clusters or --n_clusters_range must be provided."
    )
if n_clusters_range:
    if len(n_clusters_range) != 2 or n_clusters_range[0] > n_clusters_range[1]:
        raise ValueError(
            "--n_clusters_range must contain exactly two values: the minimum and maximum number of clusters."
        )
    if n_clusters is not None and not (
        n_clusters_range[0] <= n_clusters <= n_clusters_range[1]
    ):
        raise ValueError(
            f"--n_clusters ({n_clusters}) must be within --n_clusters_range {n_clusters_range}."
        )

## Read in data
logger.info("Reading input data from %s...", par["input"])
mdata = mu.read_h5mu(par["input"])
if par["modality"] not in mdata.mod:
    raise ValueError(
        f"Modality '{par['modality']}' not found in input. Available modalities: {list(mdata.mod.keys())}."
    )
adata = mdata.mod[par["modality"]]

features_key = par["input_obsm_features"]
if features_key not in adata.obsm:
    raise ValueError(
        f"--input_obsm_features '{features_key}' not found in .obsm of modality '{par['modality']}'. "
        "Make sure the neighborhood-aggregated features were computed, e.g. with cellcharter/aggregate_neighbors."
    )
features = adata.obsm[features_key]
if sp.issparse(features):
    raise ValueError(
        f".obsm['{features_key}'] is sparse. Clustering requires a dense, low-dimensional representation."
    )
features = np.asarray(features, dtype=np.float32)
logger.info("Clustering %s cells using %s features.", *features.shape)

## Set up reproducibility and compute resources
# cellcharter does not use its random_state argument to seed the model initialization,
# so seed all random number generators explicitly.
seed_everything(par["random_state"], verbose=False)
if meta["cpus"]:
    torch.set_num_threads(meta["cpus"])

# Use a single device: when using multiple devices, predictions are only returned for a subset of the cells.
trainer_params = {
    "accelerator": par["accelerator"],
    "devices": 1,
    "max_epochs": par["max_epochs"],
    "enable_progress_bar": False,
}
model_params = {
    "covariance_type": par["covariance_type"],
    "init_strategy": par["init_strategy"],
    "convergence_tolerance": par["convergence_tolerance"],
    "covariance_regularization": par["covariance_regularization"],
    "batch_size": par["batch_size"],
    "trainer_params": trainer_params,
    "random_state": par["random_state"],
}

# Fit on a minimal AnnData, so that cellcharter does not store its parameters in the output .uns
adata_features = ad.AnnData(obs=adata.obs[[]], obsm={"features": features})

## Cluster
if n_clusters_range:
    logger.info(
        "Running the stability analysis for %s to %s clusters with at most %s runs...",
        *n_clusters_range,
        par["max_runs"],
    )
    autok = cc.tl.ClusterAutoK(
        n_clusters=tuple(n_clusters_range),
        max_runs=par["max_runs"],
        convergence_tol=par["stability_convergence_tolerance"],
        model_params=model_params,
    )
    autok.fit(adata_features, use_rep="features")

    k = autok.best_k if n_clusters is None else n_clusters
    logger.info(
        "Number of clusters with the highest stability: %s. Peaks of the stability curve: %s.",
        autok.best_k,
        list(autok.peaks),
    )
    logger.info("Clustering cells into %s clusters...", k)
    labels = autok.predict(adata_features, use_rep="features", k=k)

    # ClusterAutoK also fits the numbers of clusters adjacent to the range to compute the stability.
    # The stability is only available for the numbers of clusters within the range.
    evaluated_n_clusters = np.array(autok.n_clusters[1:-1])
    stability = np.asarray(autok.stability)
    adata.uns[par["output_uns_stability"]] = {
        "n_clusters": evaluated_n_clusters,
        "stability": stability,
        "stability_mean": stability.mean(axis=1),
        "best_k": int(autok.best_k),
        "peaks": np.asarray(autok.peaks, dtype=int),
        "k": int(k),
    }
else:
    logger.info("Clustering cells into %s clusters...", n_clusters)
    model = cc.tl.Cluster(n_clusters=n_clusters, **model_params)
    model.fit(adata_features, use_rep="features")
    labels = model.predict(adata_features, use_rep="features")

adata.obs[par["output_obs_cluster"]] = labels
logger.info(
    "Cluster sizes: %s", adata.obs[par["output_obs_cluster"]].value_counts().to_dict()
)

## Save output
logger.info("Writing output to %s...", par["output"])
mdata.write_h5mu(par["output"], compression=par["output_compression"])
