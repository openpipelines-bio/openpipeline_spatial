import sys
import mudata as mu
import scanpy as sc

from banksy.initialize_banksy import initialize_banksy
from banksy.embed_banksy import generate_banksy_matrix

## VIASH START
par = {
    "input": "resources_test/xenium/xenium_tiny.h5mu",
    "modality": "rna",
    "layer": None,
    "input_obsm_spatial_coords": "spatial",
    "var_input": None,
    "k_geom": 15,
    "nbr_weight_decay": "scaled_gaussian",
    "max_m": 1,
    "lambda_param": 0.8,
    "pca_dims": 20,
    "resolution": 1.0,
    "random_state": 0,
    "output": "foo.h5mu",
    "output_obsm_embedding": "banksy_embedding",
    "output_obs_cluster": "banksy_cluster",
    "output_compression": None,
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from subset_vars import subset_vars

logger = setup_logger()

## Read in data
logger.info("Reading input data...")
adata = mu.read_h5ad(par["input"], mod=par["modality"])

# Work on a (possibly gene-subsetted) copy; the original adata (full gene
# space) is what the embedding/cluster labels get written back onto below.
if par["var_input"]:
    logger.info(f"Subsetting features on .var['{par['var_input']}']...")
    model_adata = subset_vars(adata, subset_col=par["var_input"])
else:
    model_adata = adata.copy()

# BANKSY always reads expression from .X, so swap in the requested layer.
if par["layer"]:
    logger.info(f"Using .layers['{par['layer']}'] as expression input...")
    model_adata.X = model_adata.layers[par["layer"]]

## Build the spatial neighbourhood graph
# coord_keys is a (obs_x_col, obs_y_col, obsm_key) triple that BANKSY's API
# expects; the first two are only used by BANKSY's own diagnostic plotting
# helpers (all disabled below via plt_*=False), so dummy names are fine here.
logger.info(
    f"Building spatial neighbourhood graph (k_geom={par['k_geom']}, "
    f"nbr_weight_decay={par['nbr_weight_decay']}, max_m={par['max_m']})..."
)
coord_keys = ("x", "y", par["input_obsm_spatial_coords"])
banksy_dict = initialize_banksy(
    model_adata,
    coord_keys,
    num_neighbours=par["k_geom"],
    nbr_weight_decay=par["nbr_weight_decay"],
    max_m=par["max_m"],
    plt_edge_hist=False,
    plt_nbr_weights=False,
    plt_agf_angles=False,
    plt_theta=False,
)

## Build the BANKSY-augmented (own expression + neighbourhood expression [+ AGF]) matrix
logger.info(f"Computing BANKSY matrix (lambda={par['lambda_param']})...")
banksy_dict, _ = generate_banksy_matrix(
    model_adata, banksy_dict, [par["lambda_param"]], par["max_m"], verbose=False
)
bm_adata = banksy_dict[par["nbr_weight_decay"]][par["lambda_param"]]["adata"]

## Downstream clustering: standard scanpy PCA -> neighbours -> Leiden.
# generate_banksy_matrix does not copy .obsm/.uns/.layers onto bm_adata, but
# none of that is needed here: PCA/neighbours/Leiden only touch .X.
n_comps = min(par["pca_dims"], min(bm_adata.shape) - 1)
if n_comps < par["pca_dims"]:
    logger.warning(
        f"Requested --pca_dims {par['pca_dims']} exceeds what the BANKSY "
        f"matrix shape {bm_adata.shape} supports; using {n_comps} instead."
    )

logger.info(f"Running PCA ({n_comps} components)...")
sc.pp.pca(bm_adata, n_comps=n_comps, random_state=par["random_state"])

logger.info("Computing neighbour graph on the BANKSY PCA embedding...")
sc.pp.neighbors(bm_adata, random_state=par["random_state"])

logger.info(f"Running Leiden clustering (resolution={par['resolution']})...")
sc.tl.leiden(
    bm_adata,
    resolution=par["resolution"],
    key_added="cluster",
    flavor="igraph",
    n_iterations=2,
    random_state=par["random_state"],
)

## Store results back onto the original (full gene space) adata.
# Both outputs are indexed on .obs, which subset_vars never touches, so no
# reindexing is needed regardless of whether --var_input was set.
logger.info("Storing results...")
adata.obsm[par["output_obsm_embedding"]] = bm_adata.obsm["X_pca"]
adata.obs[par["output_obs_cluster"]] = bm_adata.obs["cluster"].to_numpy()

logger.info(f"Writing output to '{par['output']}'...")
mdata = mu.MuData({par["modality"]: adata})
mdata.write_h5mu(par["output"], compression=par["output_compression"])
