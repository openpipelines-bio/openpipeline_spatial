import sys
import warnings

import mudata as mu
import numpy as np
import squidpy as sq
from scipy.spatial import ConvexHull, Voronoi, distance_matrix
from sklearn.neighbors import KernelDensity

## VIASH START
par = {
    "input": "resources_test/xenium/xenium_tiny_qc.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "spatial_coords_key": "spatial",
    "density_bandwidth": 50.0,
    "calculate_ripley_l": False,
    "n_subsample_ripley": -1,
    "output_prefix": "spatial_",
}
## VIASH END


def calculate_morphology_metrics(adata, prefix):
    """Calculate cell morphology ratios and features."""
    print("  Calculating morphology metrics...", flush=True)

    if "cell_area" not in adata.obs or "nucleus_area" not in adata.obs:
        warnings.warn(
            "cell_area and/or nucleus_area not found in .obs. "
            "Skipping morphology metrics."
        )
        return

    # Nucleus to cell area ratio - can be removed once https://github.com/openpipelines-bio/openpipeline_qc/issues/18 is fixed.
    adata.obs[f"{prefix}nucleus_cell_ratio"] = (
        adata.obs["nucleus_area"] / adata.obs["cell_area"]
    )

    # Cell area percentile (relative size)
    adata.obs[f"{prefix}cell_area_percentile"] = (
        adata.obs["cell_area"].rank(pct=True) * 100
    )

    print(
        f"    Added: {prefix}nucleus_cell_ratio, {prefix}cell_area_percentile",
        flush=True,
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


def calculate_density_metrics(adata, spatial_coords, bandwidth, prefix):
    """Calculate local density metrics."""
    print("  Calculating density metrics...", flush=True)

    # Kernel density estimation
    kde = KernelDensity(bandwidth=bandwidth, kernel="gaussian")
    kde.fit(spatial_coords)
    log_density = kde.score_samples(spatial_coords)
    adata.obs[f"{prefix}kernel_density"] = np.exp(log_density)

    # Calculate degree centrality (as a proxy for local density / number of neighbors)
    # This assumes spatial_connectivities is already present in .obsp (from build_spatial_graph)
    if "spatial_connectivities" in adata.obsp:
        spatial_connectivities = adata.obsp["spatial_connectivities"]
        degrees = spatial_connectivities.sum(axis=1)
        # If sparse matrix, it returns matrix object, need to convert to array
        if hasattr(degrees, "A1"):
            degrees = degrees.A1
        adata.obs[f"{prefix}graph_degree"] = degrees

        print("    Added: kernel_density, graph_degree", flush=True)
    else:
        print(
            "    Warning: 'spatial_connectivities' not found in .obsp. Skipping graph_degree.",
            flush=True,
        )
        print("    Added: kernel_density", flush=True)


def calculate_voronoi_metrics(adata, spatial_coords, prefix):
    """Calculate Voronoi tessellation statistics."""
    print("  Calculating Voronoi tessellation...", flush=True)

    try:
        vor = Voronoi(spatial_coords)

        polygon_areas = []
        neighbor_counts = []

        for point_idx in range(len(spatial_coords)):
            region_idx = vor.point_region[point_idx]
            region = vor.regions[region_idx]

            # Skip infinite regions
            if -1 in region or len(region) == 0:
                polygon_areas.append(np.nan)
                neighbor_counts.append(np.nan)
                continue

            # Calculate polygon area using shoelace formula
            vertices = vor.vertices[region]
            x = vertices[:, 0]
            y = vertices[:, 1]
            area = 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))
            polygon_areas.append(area)

            # Count neighbors (shared vertices)
            neighbor_counts.append(len(region))

        adata.obs[f"{prefix}voronoi_area"] = polygon_areas
        adata.obs[f"{prefix}voronoi_neighbors"] = neighbor_counts

        print("    Added: voronoi_area, voronoi_neighbors", flush=True)
    except Exception as e:
        warnings.warn(f"Could not calculate Voronoi tessellation: {e}")


def calculate_global_statistics(adata, spatial_coords, par):
    """Calculate global spatial pattern statistics."""
    print("  Calculating global spatial statistics...", flush=True)

    stats = {}

    # Global density
    try:
        hull = ConvexHull(spatial_coords)
        area = hull.volume  # In 2D, volume of convex hull is the area
        stats["area_calculation_method"] = "convex_hull"
        stats["cell_density"] = len(spatial_coords) / area
        stats["total_area"] = area
        stats["n_cells"] = len(spatial_coords)
    except Exception as e:
        print(f"    Error: Could not calculate convex hull area ({e})", flush=True)

    # Ripley's L (using squidpy)
    if par["calculate_ripley_l"]:
        print("    Computing Ripley's L function (Squidpy)...", flush=True)

        n_subsample = par["n_subsample_ripley"]

        try:
            # Create a working AnnData for Ripley's
            if n_subsample > 0 and len(adata) > n_subsample:
                print(f"      Subsampling to {n_subsample} random cells...", flush=True)

                # Use numpy to generate random indices
                indices = np.random.choice(len(adata), n_subsample, replace=False)
                adata_ripley = adata[indices].copy()
                stats["ripley_l_subsampled"] = True
            else:
                adata_ripley = adata.copy()
                stats["ripley_l_subsampled"] = False

            # Squidpy requires a cluster key. We create a dummy one for "global" context.
            adata_ripley.obs["_temp_global"] = "all"
            adata_ripley.obs["_temp_global"] = adata_ripley.obs["_temp_global"].astype(
                "category"
            )

            # Calculate Ripley's L statistic
            # This stores the result in adata.uns['all_L']
            sq.gr.ripley(
                adata_ripley,
                cluster_key="_temp_global",
                mode="L",
                spatial_key=par["spatial_coords_key"],
                n_simulations=50,
            )

            # Extract basic stats from the results
            if "all_L" in adata_ripley.uns and "L_stat" in adata_ripley.uns["all_L"]:
                l_stat_df = adata_ripley.uns["all_L"]["L_stat"]
                stats["ripley_l_max"] = float(l_stat_df.max().max())
                stats["ripley_l_mean"] = float(l_stat_df.mean().mean())

            # Clean up
            if "all_L" in adata_ripley.uns:
                del adata_ripley.uns["all_L"]

        except Exception as e:
            print(f"    Error calculating Ripley's L: {e}", flush=True)
            import traceback

            traceback.print_exc()
    else:
        print("    Skipping Ripley's L (disabled in config)...", flush=True)

    # Spatial extent
    min_coords = spatial_coords.min(axis=0)
    max_coords = spatial_coords.max(axis=0)

    stats["spatial_extent_x"] = float(max_coords[0] - min_coords[0])
    stats["spatial_extent_y"] = float(max_coords[1] - min_coords[1])
    stats["centroid_x"] = float(spatial_coords[:, 0].mean())
    stats["centroid_y"] = float(spatial_coords[:, 1].mean())

    if "cell_density" in stats:
        print(
            f"    Global stats: cell_density={stats['cell_density']:.4f}",
            flush=True,
        )

    return stats


def main(par):
    print(f"\n>>> Reading MuData from '{par['input']}'...", flush=True)
    mdata = mu.read_h5mu(par["input"])
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
        f"\n>>> Extracting spatial coordinates from .obsm['{par['spatial_coords_key']}']...",
        flush=True,
    )
    if par["spatial_coords_key"] not in adata.obsm:
        raise KeyError(
            f"Spatial key '{par['spatial_coords_key']}' not found in .obsm. "
            f"Available keys: {list(adata.obsm.keys())}"
        )

    spatial_coords = adata.obsm[par["spatial_coords_key"]]
    if spatial_coords.shape[1] != 2:
        raise ValueError(
            f"Expected 2D spatial coordinates, got shape {spatial_coords.shape}"
        )
    print(f"  Shape: {spatial_coords.shape} (n_cells × 2)", flush=True)

    prefix = par["output_prefix"]

    # Calculate morphology metrics
    print("\n>>> Calculating morphology metrics...", flush=True)
    calculate_morphology_metrics(adata, prefix)

    # Calculate position features
    print("\n>>> Calculating position-based features...", flush=True)
    calculate_position_features(adata, spatial_coords, prefix)

    # Calculate density metrics
    print("\n>>> Calculating density metrics...", flush=True)
    calculate_density_metrics(
        adata,
        spatial_coords,
        par["density_bandwidth"],
        prefix,
    )

    # Calculate Voronoi tessellation
    print("\n>>> Calculating Voronoi tessellation...", flush=True)
    calculate_voronoi_metrics(adata, spatial_coords, prefix)

    # Calculate global statistics
    print("\n>>> Calculating global spatial statistics...", flush=True)
    global_stats = calculate_global_statistics(adata, spatial_coords, par)

    # Store in uns
    if "spatial_stats" not in adata.uns:
        adata.uns["spatial_stats"] = {}
    adata.uns["spatial_stats"].update(global_stats)

    print(f"\n>>> Writing output to '{par['output']}'...", flush=True)
    mdata.write_h5mu(par["output"])

    print("\n>>> Done!\n", flush=True)


if __name__ == "__main__":
    sys.exit(main(par))
