"""Crop a converted Xenium SpatialData object down to one dense patch of cells.

Unlike ``subset_xenium.py`` (which crops the raw 10x ``outs/`` files *before*
conversion), this operates on an *already-converted* SpatialData Zarr store, using
SpatialData's own ``bounding_box_query()``. That sidesteps the biggest problem with
cropping raw Xenium output: ``spatialdata_io.xenium()`` reads a proprietary,
undocumented, version-branched file (``cells.zarr.zip``) to build the pixel-space
segmentation label rasters and cross-check the cell table, and cropping that file
correctly would mean reverse-engineering its internal indexing invariants. Converting
first means ``cells.zarr.zip`` is parsed exactly once, by ``spatialdata_io`` itself,
into SpatialData's own standard ``Labels2D``/``AnnData`` representation -- and
*that* is generically and correctly croppable via the public API, no custom format
handling needed. The result keeps everything: images, raster labels, boundary
shapes, transcripts, and the annotation table, all consistently subset together.

Steps:

1. Cluster cell centroids (``table.obsm["spatial"]``, in microns) with ``DBSCAN`` to
   find spatially separated patches -- useful for source datasets (like 10x's XOA
   v4.0 example data) that were themselves "artificially subset to N square patches"
   within one larger image, so a plain bounding box over *all* cells would span the
   (mostly empty) gaps between patches instead of just one dense region.
2. Pick one patch (largest, by default) and compute its micron bounding box with a
   margin.
3. Convert that box to pixel coordinates and run ``bounding_box_query()``.

Coordinate systems gotcha: elements from ``spatialdata_io.xenium()`` don't all use
the same units under the shared ``"global"`` coordinate system. Images and labels use
an ``Identity`` transform (so ``"global"`` is pixel space for them), while shapes and
points use a ``Scale(1 / pixel_size)`` transform (their values are stored in microns,
scaled *into* that same pixel-space ``"global"``). Querying with micron-valued bounds
silently returns ``None`` for shapes/points (correctly read as "no matches") while
silently cropping the *wrong* region for images/labels (no error -- your micron
numbers just get reinterpreted as pixel indices). This script extracts ``pixel_size``
straight from a shape element's own transform and converts before querying, so this
mismatch never has to be hand-tracked by the caller.

Usage:
    crop_xenium_to_patch.py --input FULL_SDATA.zarr --output PATCH.zarr [--patch-rank 0]
"""

import argparse

import numpy as np
import spatialdata as sd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from spatialdata import bounding_box_query
from spatialdata.transformations import get_transformation


def cluster_by_distance(coords, eps):
    """Group points into clusters by single-linkage distance threshold.

    Equivalent in spirit to DBSCAN's core operation (points within `eps` of each
    other join the same cluster) via connected components on a radius graph, using
    only scipy -- scikit-learn isn't part of this project's pinned Python
    environments, and dragging it in as a one-off dependency for a resource script
    isn't worth it when scipy alone covers this case (patches are dense and well
    separated, so DBSCAN's noise/core-point distinction isn't needed; isolated
    points just end up as their own singleton cluster and get filtered out by
    `--min-samples` at the ranking step instead).
    """
    n = len(coords)
    pairs = cKDTree(coords).query_pairs(r=eps, output_type="ndarray")
    if len(pairs) == 0:
        return np.arange(n)
    graph = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    _, labels = connected_components(graph, directed=False)
    return labels


def get_pixel_size(sdata):
    """Recover microns-per-pixel from a shape element's Scale transform to 'global'.

    Shapes/points elements are stored in microns and transformed into the shared
    pixel-space 'global' coordinate system via Scale(1 / pixel_size); the x-axis
    scale factor is the reciprocal of that.
    """
    shapes_name = next(iter(sdata.shapes))
    transform = get_transformation(sdata.shapes[shapes_name], get_all=True)["global"]
    return 1.0 / transform.scale[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", required=True, help="converted Xenium SpatialData .zarr"
    )
    parser.add_argument("--output", required=True, help="output cropped .zarr")
    parser.add_argument(
        "--margin-um",
        type=float,
        default=10.0,
        help="padding (in microns) around the selected patch's cell bounding box",
    )
    parser.add_argument(
        "--eps-um",
        type=float,
        default=30.0,
        help="DBSCAN neighborhood radius (microns) for patch clustering; cells "
        "within this distance of each other are considered the same patch",
    )
    parser.add_argument(
        "--min-samples",
        type=int,
        default=5,
        help="DBSCAN min_samples: minimum cells to form a patch (smaller clusters "
        "and isolated cells are treated as noise and excluded from selection)",
    )
    parser.add_argument(
        "--patch-rank",
        type=int,
        default=0,
        help="which patch to keep, ranked by cell count descending (0 = largest)",
    )
    args = parser.parse_args()

    print(f">>> reading {args.input}")
    sdata = sd.read_zarr(args.input)
    table = sdata.tables["table"]
    pixel_size = get_pixel_size(sdata)
    print(f"  {table.n_obs} cells total, pixel_size={pixel_size}")

    print(">>> clustering cell centroids into patches")
    coords_um = table.obsm["spatial"]
    labels = cluster_by_distance(coords_um, eps=args.eps_um)
    all_ids, all_counts = np.unique(labels, return_counts=True)
    keep = all_counts >= args.min_samples
    patch_ids, counts = all_ids[keep], all_counts[keep]
    order = np.argsort(-counts)
    patch_ids, counts = patch_ids[order], counts[order]
    for patch_id, count in zip(patch_ids, counts):
        print(f"  patch {patch_id}: {count} cells")
    n_noise = int(all_counts[~keep].sum())
    if n_noise:
        print(f"  {n_noise} cells in clusters smaller than --min-samples (excluded)")

    if args.patch_rank >= len(patch_ids):
        raise ValueError(
            f"--patch-rank {args.patch_rank} out of range: only {len(patch_ids)} "
            "patches found"
        )
    chosen_patch = patch_ids[args.patch_rank]
    patch_coords = coords_um[labels == chosen_patch]
    print(f">>> selected patch {chosen_patch} ({len(patch_coords)} cells)")

    x_min = patch_coords[:, 0].min() - args.margin_um
    x_max = patch_coords[:, 0].max() + args.margin_um
    y_min = patch_coords[:, 1].min() - args.margin_um
    y_max = patch_coords[:, 1].max() + args.margin_um
    print(
        f"  bbox (um): x=[{x_min:.1f},{x_max:.1f}] y=[{y_min:.1f},{y_max:.1f}] "
        f"-> ({(x_max - x_min) / pixel_size:.0f}x{(y_max - y_min) / pixel_size:.0f} px)"
    )

    print(">>> cropping with bounding_box_query")
    cropped = bounding_box_query(
        sdata,
        axes=("x", "y"),
        min_coordinate=[x_min / pixel_size, y_min / pixel_size],
        max_coordinate=[x_max / pixel_size, y_max / pixel_size],
        target_coordinate_system="global",
    )
    print(
        f"  {cropped.tables['table'].n_obs} cells, {len(cropped.points['transcripts'])} transcripts kept"
    )

    print(f">>> writing {args.output}")
    cropped.write(args.output, overwrite=True)


if __name__ == "__main__":
    main()
