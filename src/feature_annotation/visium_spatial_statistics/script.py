import sys
import warnings

import mudata as md
import numpy as np

from scipy.spatial import ConvexHull, distance_matrix

## VIASH START
par = {
    "input": "resources_test/visium/visium_tiny_neighbours.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "obsm_spatial_coordinates": "spatial",
    "obsp_spatial_graph": "spatial_connectivities",
    "obs_total_counts": "total_counts",
    "output_prefix": "spatial_",
    "tissue_edge_max_neighbors": 6,
}
## VIASH END


def calculate_neighbors_metrics(
    adata, conn, prefix, obs_total_counts, tissue_edge_max_neighbors
):
    """Add per-spot neighbour metrics to adata.obs."""
    adata.obs[f"{prefix}n_neighbors"] = np.asarray(conn.sum(axis=1)).flatten()
    adata.obs[f"{prefix}tissue_edge"] = (
        adata.obs[f"{prefix}n_neighbors"] < tissue_edge_max_neighbors
    )

    if obs_total_counts in adata.obs:
        counts = adata.obs[obs_total_counts].values
        adata.obs[f"{prefix}local_expression_density"] = np.asarray(
            conn @ counts
        ).flatten()
    else:
        warnings.warn(
            f"'{obs_total_counts}' not found in .obs; skipping local_expression_density"
        )


def calculate_position_features(adata, spatial_coords, prefix):
    """Calculate position-based features."""
    print("  Calculating position-based features...", flush=True)

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
    try:
        hull = ConvexHull(spatial_coords)
        hull_points = spatial_coords[hull.vertices]

        # For each point, find distance to nearest hull vertex (approximation)
        distances_to_boundary = np.min(
            distance_matrix(spatial_coords, hull_points), axis=1
        )
        adata.obs[f"{prefix}distance_to_boundary"] = distances_to_boundary
        print(
            "    Added: distance_to_centroid, norm_x, norm_y, distance_to_boundary",
            flush=True,
        )
    except Exception as e:
        warnings.warn(f"Could not calculate convex hull: {e}")
        print("    Added: distance_to_centroid, norm_x, norm_y", flush=True)


def calculate_global_statistics(spatial_coords):
    """Calculate global spatial statistics."""
    print("  Calculating global spatial statistics...", flush=True)

    stats = {}

    # Global density
    try:
        hull = ConvexHull(spatial_coords)
        area = hull.volume  # In 2D, volume of convex hull is the area
        stats["area_calculation_method"] = "convex_hull"
        stats["spot_density"] = len(spatial_coords) / area
        stats["total_area"] = area
        stats["n_spots"] = len(spatial_coords)
    except Exception as e:
        print(f"    Error: Could not calculate convex hull area ({e})", flush=True)

    # Spatial extent
    min_coords = spatial_coords.min(axis=0)
    max_coords = spatial_coords.max(axis=0)
    stats["spatial_extent_x"] = float(max_coords[0] - min_coords[0])
    stats["spatial_extent_y"] = float(max_coords[1] - min_coords[1])
    stats["centroid_x"] = float(spatial_coords[:, 0].mean())
    stats["centroid_y"] = float(spatial_coords[:, 1].mean())

    if "spot_density" in stats:
        print(
            f"    Global stats: spot_density={stats['spot_density']:.4f}",
            flush=True,
        )

    return stats


def main(par):
    print("====== Visium spatial statistics ======", flush=True)

    print(f"\n>>> Reading MuData from '{par['input']}'...", flush=True)
    mdata = md.read_h5mu(par["input"])
    print(mdata, flush=True)

    print(f"\n>>> Extracting modality '{par['modality']}'...", flush=True)
    if par["modality"] not in mdata.mod:
        raise KeyError(
            f"Modality '{par['modality']}' not found in MuData. "
            f"Available modalities: {list(mdata.mod.keys())}"
        )
    adata = mdata[par["modality"]]
    print(adata, flush=True)

    print(
        f"\n>>> Extracting spatial coordinates from .obsm['{par['obsm_spatial_coordinates']}']...",
        flush=True,
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
    print(f"  Shape: {spatial_coords.shape} (n_spots x 2)", flush=True)

    prefix = par["output_prefix"]

    print(
        f"\n>>> Extracting spatial graph from .obsp['{par['obsp_spatial_graph']}']...",
        flush=True,
    )
    if par["obsp_spatial_graph"] not in adata.obsp:
        raise KeyError(
            f"Spatial graph key '{par['obsp_spatial_graph']}' not found in .obsp. "
            f"Available keys: {list(adata.obsp.keys())}"
        )
    conn = adata.obsp[par["obsp_spatial_graph"]]
    print(f"  Shape: {conn.shape} (n_spots x n_spots)", flush=True)

    print("\n>>> Calculating neighbor metrics...", flush=True)
    calculate_neighbors_metrics(
        adata,
        conn,
        prefix,
        par["obs_total_counts"],
        par["tissue_edge_max_neighbors"],
    )

    print("\n>>> Calculating position-based features...", flush=True)
    calculate_position_features(adata, spatial_coords, prefix)

    print("\n>>> Calculating global spatial statistics...", flush=True)
    global_stats = calculate_global_statistics(spatial_coords)

    # Store in uns
    if "spatial_stats" not in adata.uns:
        adata.uns["spatial_stats"] = {}
    adata.uns["spatial_stats"].update(global_stats)

    mdata.mod[par["modality"]] = adata

    print(f"\n>>> Writing output to '{par['output']}'...", flush=True)
    mdata.write_h5mu(par["output"])

    print("\n>>> Done!\n", flush=True)


if __name__ == "__main__":
    sys.exit(main(par))
