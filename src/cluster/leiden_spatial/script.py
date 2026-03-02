import sys
import scanpy as sc
import mudata as mu

## VIASH START
par = {
    # Inputs
    "input": "resources_test/xenium/xenium_tiny.qc.neighbors.h5mu",
    "modality": "rna",
    "input_obsp_expression_connectivities": "connectivities",
    "input_obsp_spatial_connectivities": "spatial_connectivities",
    # Clustering options
    "alpha": 0.2,
    "resolution": 1.0,
    "obs_label": "leiden_spatial",
    # Outputs
    "output": "foo.h5mu",
    "output_compression": None,
}
meta = {"resources_dir": "src/utils/"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

## Read data
logger.info("Reading input data...")
adata = mu.read_h5ad(par["input"], mod=par["modality"])

## Validate spatial graph
spatial_key = par["input_obsp_spatial_connectivities"]
if spatial_key not in adata.obsp:
    raise ValueError(
        f"Spatial connectivities key '{spatial_key}' not found in .obsp. "
        "Run the neighbors/spatial_neighborhood_graph component first."
    )

## Validate and read expression connectivity graph
expr_key = par["input_obsp_expression_connectivities"]
if expr_key not in adata.obsp:
    raise ValueError(
        f"Expression connectivities key '{expr_key}' not found in .obsp. "
        "Run a neighbors component before this component."
    )
logger.info(f"Using expression graph from .obsp['{expr_key}']...")
nn_graph_genes = adata.obsp[expr_key]
nn_graph_space = adata.obsp[spatial_key]

## Combine graphs
alpha = par["alpha"]
logger.info(
    f"Combining graphs (alpha={alpha}: spatial weight, "
    f"{1 - alpha:.2f}: expression weight)..."
)
joint_graph = (1 - alpha) * nn_graph_genes + alpha * nn_graph_space

## Run Leiden
logger.info(
    f"Running Leiden clustering (resolution={par['resolution']})..."
)
sc.tl.leiden(
    adata,
    adjacency=joint_graph,
    resolution=par["resolution"],
    key_added=par["obs_label"],
    flavor="igraph",
    n_iterations=2,
)

n_domains = adata.obs[par["obs_label"]].nunique()
logger.info(f"Identified {n_domains} spatial domain(s), stored in .obs['{par['obs_label']}'].")

## Write output
logger.info("Saving output data...")
mdata = mu.MuData({par["modality"]: adata})
mdata.write_h5mu(par["output"], compression=par["output_compression"])
