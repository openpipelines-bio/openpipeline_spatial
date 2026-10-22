import sys
import anndata as ad
import mudata as mu
import cellcharter as cc
import scipy.sparse as sp

## VIASH START
par = {
    # Inputs
    "input": "resources_test/xenium/xenium_tiny.qc.all_neighbors.pca.h5mu",
    "modality": "rna",
    "input_obsm_features": "X_pca",
    "layer": None,
    "input_obsp_spatial_connectivities": "spatial_connectivities",
    "input_obs_sample_key": None,
    # Neighborhood aggregation
    "n_layers": [3],
    "aggregations": ["mean"],
    # Outputs
    "output": "output.h5mu",
    "output_obsm_aggregated": "X_cellcharter",
    "output_compression": None,
}
meta = {"resources_dir": "src/utils/"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

## Read in data
logger.info("Reading input data from %s...", par["input"])
mdata = mu.read_h5mu(par["input"])
if par["modality"] not in mdata.mod:
    raise ValueError(
        f"Modality '{par['modality']}' not found in input. Available modalities: {list(mdata.mod.keys())}."
    )
adata = mdata.mod[par["modality"]]

## Select the features to aggregate
if par["input_obsm_features"] and par["layer"]:
    raise ValueError("Only one of --input_obsm_features and --layer can be provided.")
if par["input_obsm_features"]:
    if par["input_obsm_features"] not in adata.obsm:
        raise ValueError(
            f"--input_obsm_features '{par['input_obsm_features']}' not found in .obsm of modality '{par['modality']}'."
        )
    logger.info("Using .obsm['%s'] as features.", par["input_obsm_features"])
    features = adata.obsm[par["input_obsm_features"]]
elif par["layer"]:
    if par["layer"] not in adata.layers:
        raise ValueError(
            f"--layer '{par['layer']}' not found in .layers of modality '{par['modality']}'."
        )
    logger.info("Using .layers['%s'] as features.", par["layer"])
    features = adata.layers[par["layer"]]
else:
    logger.info("Using .X as features.")
    features = adata.X

# cellcharter does not support sparse features in combination with a sample key or the 'var' aggregation.
# Densify the input.
if sp.issparse(features):
    logger.info("Converting sparse features to a dense matrix...")
    features = features.toarray()

## Validate the spatial graph and sample key
connectivity_key = par["input_obsp_spatial_connectivities"]
if connectivity_key not in adata.obsp:
    raise ValueError(
        f"--input_obsp_spatial_connectivities '{connectivity_key}' not found in .obsp of modality '{par['modality']}'. "
        "Make sure a spatial neighborhood graph was computed, e.g. with neighbors/spatial_neighborhood_graph."
    )

sample_key = par["input_obs_sample_key"]
if sample_key and sample_key not in adata.obs.columns:
    raise ValueError(
        f"--input_obs_sample_key '{sample_key}' not found in .obs of modality '{par['modality']}'."
    )

## Validate the neighborhood layers
# A single value n means "all layers up to n"; multiple values are explicit hop distances
# for which the cells' own features (layer 0) need to be included.
n_layers = par["n_layers"]
if len(n_layers) == 1:
    n_layers = n_layers[0]
else:
    n_layers = sorted(set(n_layers))
    if 0 not in n_layers:
        raise ValueError(
            "When providing multiple values for --n_layers, 0 (the cells' own features) must be included."
        )

## Aggregate neighbors
# Build a minimal AnnData holding only what cellcharter needs. This allows using .layers as input,
# and avoids cellcharter implicitly using an .obs column named 'sample' when no sample key was provided.
adata_features = ad.AnnData(
    obs=adata.obs[[sample_key]] if sample_key else adata.obs[[]],
    obsm={"features": features},
    obsp={connectivity_key: adata.obsp[connectivity_key]},
)

logger.info(
    "Aggregating neighbors with n_layers=%s and aggregations=%s...",
    n_layers,
    par["aggregations"],
)
aggregated = cc.gr.aggregate_neighbors(
    adata_features,
    n_layers=n_layers,
    aggregations=par["aggregations"],
    connectivity_key=connectivity_key,
    use_rep="features",
    sample_key=sample_key,
    copy=True,
)
logger.info("Aggregated features have shape %s.", aggregated.shape)
adata.obsm[par["output_obsm_aggregated"]] = aggregated

## Save output
logger.info("Writing output to %s...", par["output"])
mdata.write_h5mu(par["output"], compression=par["output_compression"])
