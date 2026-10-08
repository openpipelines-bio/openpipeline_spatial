import sys

import h5py
import mudata as md
import numpy as np

from scipy.spatial import ConvexHull

## VIASH START
par = {
    "input": "resources_test/visium/visium_tiny_neighbours.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "obsm_spatial_coordinates": "spatial",
    "obsp_spatial_graph": "spatial_connectivities",
    "obs_total_counts": "total_counts",
    "output_prefix": "spatial_",
    "uns_spatial_stats": "spatial_stats",
    "output_compression": None,
    "tissue_edge_max_neighbors": 6,
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()


def calculate_neighbors_metrics(
    adata, conn, prefix, obs_total_counts, tissue_edge_max_neighbors
):
    """Add per-spot neighbour metrics to adata.obs."""
    adata.obs[f"{prefix}n_neighbors"] = conn.count_nonzero(axis=1)
    adata.obs[f"{prefix}tissue_edge"] = (
        adata.obs[f"{prefix}n_neighbors"] < tissue_edge_max_neighbors
    )

    if obs_total_counts in adata.obs:
        counts = adata.obs[obs_total_counts].values
        # Binarize the graph so the counts of neighbouring spots are summed
        # as-is, rather than scaled by edge weights or distances
        adjacency = (conn != 0).astype(float)
        adata.obs[f"{prefix}local_expression_density"] = np.asarray(
            adjacency @ counts
        ).flatten()
    else:
        logger.warning(
            f"'{obs_total_counts}' not found in .obs; skipping local_expression_density"
        )


def calculate_position_features(adata, spatial_coords, prefix):
    """Calculate position-based features."""
    # Tissue centroid
    centroid = spatial_coords.mean(axis=0)

    # Distance to centroid
    distances_to_centroid = np.linalg.norm(spatial_coords - centroid, axis=1)
    adata.obs[f"{prefix}distance_to_centroid"] = distances_to_centroid

    # Normalized coordinates (0-1 scale)
    min_coords = spatial_coords.min(axis=0)
    max_coords = spatial_coords.max(axis=0)
    coord_range = max_coords - min_coords

    normalized_coords = (spatial_coords - min_coords) / coord_range
    adata.obs[f"{prefix}norm_x"] = normalized_coords[:, 0]
    adata.obs[f"{prefix}norm_y"] = normalized_coords[:, 1]

    # Distance to convex hull boundary
    hull = ConvexHull(spatial_coords)

    # Each row of hull.equations is an outward edge normal and offset, so
    # -(normal . x + offset) is the distance from an interior point to that
    # edge; the nearest edge gives the distance to the boundary.
    homogeneous_coords = np.c_[spatial_coords, np.ones(len(spatial_coords))]
    distances_to_boundary = np.clip(
        -(hull.equations @ homogeneous_coords.T).max(axis=0), 0, None
    )
    adata.obs[f"{prefix}distance_to_boundary"] = distances_to_boundary
    logger.info("Added: distance_to_centroid, norm_x, norm_y, distance_to_boundary")


def calculate_global_statistics(spatial_coords):
    """Calculate global spatial statistics."""
    stats = {}

    # Global density
    hull = ConvexHull(spatial_coords)
    area = hull.volume  # In 2D, volume of convex hull is the area
    stats["area_calculation_method"] = "convex_hull"
    stats["spot_density"] = len(spatial_coords) / area
    stats["total_area"] = area
    stats["n_spots"] = len(spatial_coords)

    # Spatial extent
    min_coords = spatial_coords.min(axis=0)
    max_coords = spatial_coords.max(axis=0)
    stats["spatial_extent_x"] = float(max_coords[0] - min_coords[0])
    stats["spatial_extent_y"] = float(max_coords[1] - min_coords[1])
    stats["centroid_x"] = float(spatial_coords[:, 0].mean())
    stats["centroid_y"] = float(spatial_coords[:, 1].mean())

    logger.info(f"Global stats: spot_density={stats['spot_density']:.4f}")

    return stats


def main(par):
    with h5py.File(par["input"], "r") as h5mu:
        available_modalities = list(h5mu["mod"].keys())
    if par["modality"] not in available_modalities:
        raise KeyError(
            f"Modality '{par['modality']}' not found in MuData. "
            f"Available modalities: {available_modalities}"
        )

    logger.info(f"Reading modality '{par['modality']}' from '{par['input']}'...")
    adata = md.read_h5ad(par["input"], mod=par["modality"])
    logger.info(adata)

    logger.info(
        f"Extracting spatial coordinates from .obsm['{par['obsm_spatial_coordinates']}']..."
    )
    if par["obsm_spatial_coordinates"] not in adata.obsm:
        raise KeyError(
            f"Spatial key '{par['obsm_spatial_coordinates']}' not found in .obsm. "
            f"Available keys: {list(adata.obsm.keys())}"
        )

    spatial_coords = adata.obsm[par["obsm_spatial_coordinates"]]
    if spatial_coords.shape[1] != 2:
        raise ValueError(
            f"Expected 2D spatial coordinates, got shape {spatial_coords.shape}"
        )
    logger.info(f"Shape: {spatial_coords.shape} (n_spots x 2)")

    prefix = par["output_prefix"]

    logger.info(
        f"Extracting spatial graph from .obsp['{par['obsp_spatial_graph']}']..."
    )
    if par["obsp_spatial_graph"] not in adata.obsp:
        raise KeyError(
            f"Spatial graph key '{par['obsp_spatial_graph']}' not found in .obsp. "
            f"Available keys: {list(adata.obsp.keys())}"
        )
    conn = adata.obsp[par["obsp_spatial_graph"]]
    logger.info(f"Shape: {conn.shape} (n_spots x n_spots)")

    logger.info("Calculating neighbor metrics...")
    calculate_neighbors_metrics(
        adata,
        conn,
        prefix,
        par["obs_total_counts"],
        par["tissue_edge_max_neighbors"],
    )

    logger.info("Calculating position-based features...")
    calculate_position_features(adata, spatial_coords, prefix)

    logger.info("Calculating global spatial statistics...")
    global_stats = calculate_global_statistics(spatial_coords)

    # Store in uns
    uns_key = par["uns_spatial_stats"]
    if uns_key not in adata.uns:
        adata.uns[uns_key] = {}
    adata.uns[uns_key].update(global_stats)

    logger.info(f"Writing output to '{par['output']}'...")
    write_h5ad_to_h5mu_with_compression(
        par["output"],
        par["input"],
        par["modality"],
        adata,
        par["output_compression"],
    )

    logger.info("Done!")


if __name__ == "__main__":
    sys.exit(main(par))
