import sys
import scanpy as sc
import mudata as mu

## VIASH START
par = {
    # Inputs
    "input": "resources_test/xenium/xenium_tiny.spatial_expression_neighbors.h5mu",
    "modality": "rna",
    "input_obsp_expression_connectivities": "connectivities",
    "input_obsp_spatial_connectivities": "spatial_connectivities",
    # Fusion options
    "alpha": 0.2,
    # Outputs
    "output": "foo.h5mu",
    "output_compression": None,
    "output_obsp_connectivities": "spatial_expression_connectivities",
}
meta = {"resources_dir": "src/utils/"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

## Read data
logger.info("Reading input data...")
adata = mu.read_h5ad(par["input"], mod=par["modality"])

## Validate inputs
spatial_key = par["input_obsp_spatial_connectivities"]
if spatial_key not in adata.obsp:
    raise ValueError(
        f"Spatial connectivities key '{spatial_key}' not found in .obsp."
    )

expr_key = par["input_obsp_expression_connectivities"]
if expr_key not in adata.obsp:
    raise ValueError(
        f"Expression connectivities key '{expr_key}' not found in .obsp."
    )

nn_graph_genes = adata.obsp[expr_key]
nn_graph_space = adata.obsp[spatial_key]

## Combine graphs
alpha = par["alpha"]
logger.info(
    f"Combining graphs (alpha={alpha}: spatial weight, "
    f"{1 - alpha:.2f}: expression weight)..."
)
joint_graph = (1 - alpha) * nn_graph_genes + alpha * nn_graph_space
out_key = par["output_obsp_connectivities"]
logger.info(f"Storing result in .obsp['{out_key}']...")
adata.obsp[out_key] = joint_graph
adata.uns[out_key] = {
    "params": {"alpha": alpha},
    "inputs": {"expression_connectivities": expr_key, "spatial_connectivities": spatial_key},
}

## Write output
logger.info("Saving output data...")
mdata = mu.MuData({par["modality"]: adata})
mdata.write_h5mu(par["output"], compression=par["output_compression"])
