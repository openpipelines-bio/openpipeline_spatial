import json
import re
import shutil
import sys
import tempfile
import zipfile
from contextlib import contextmanager
from pathlib import Path

import h5py
import numcodecs
import numpy as np
import pandas as pd
import tifffile
import zarr
from scipy.sparse import coo_matrix, csc_matrix, csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

## VIASH START
par = {
    "input": "resources_test/xenium/xenium_multicellseg_tiny",
    "output": "xenium_multicellseg_tiny_cropped",
    "margin_um": 10.0,
    "eps_um": 30.0,
    "min_cells_per_patch": 5,
    "patch_rank": 0,
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

# Xenium's *.zarr.zip stores are zarr v2 (.zgroup/.zarray) with Blosc/zstd
# compression. xeniumranger can only read v2, so every store is written back in
# that format, rather than zarr-python 3's default v3 layout.
ZARR_COMPRESSOR = numcodecs.Blosc(
    cname="zstd", clevel=5, shuffle=numcodecs.Blosc.SHUFFLE
)

# Transcripts with a Phred-scaled quality value at or above this threshold are
# "high quality" in XOA's outputs (transcripts.zarr.zip tiles, density grids).
HIGH_QV_THRESHOLD = 20


def cluster_by_distance(coords, eps):
    """Group points into clusters by single-linkage distance threshold.

    Equivalent in spirit to DBSCAN's core operation (points within `eps` of
    each other join the same cluster) via connected components on a radius
    graph, using only scipy: patches are dense and well separated, so DBSCAN's
    noise/core-point distinction isn't needed, isolated points just end up as
    their own singleton cluster and get filtered out by min_cells_per_patch at
    the ranking step instead.
    """
    n = len(coords)
    pairs = cKDTree(coords).query_pairs(r=eps, output_type="ndarray")
    if len(pairs) == 0:
        return np.arange(n)
    graph = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    _, labels = connected_components(graph, directed=False)
    return labels


def cell_id_str_from_prefix_suffix(prefix, suffix):
    """Decode a Xenium (cell_id_prefix, dataset_suffix) uint32 pair into its string cell ID.

    See https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/analysis/xoa-output-zarr#cellID
    """
    hex_shift = {str(i): chr(ord("a") + i) for i in range(10)} | {
        chr(ord("a") + i): chr(ord("a") + 10 + i) for i in range(6)
    }
    prefix_hex = [hex(int(x))[2:] for x in prefix]
    prefix_shifted = ["".join(hex_shift[c] for c in x) for x in prefix_hex]
    return np.array([p.rjust(8, "a") + f"-{s}" for p, s in zip(prefix_shifted, suffix)])


@contextmanager
def read_zarr_zip(path):
    with tempfile.TemporaryDirectory() as extract_dir:
        with zipfile.ZipFile(path) as zf:
            zf.extractall(extract_dir)
        yield zarr.open_group(extract_dir, mode="r")


@contextmanager
def write_zarr_zip(path):
    """Yield a fresh zarr v2 group, zipped (uncompressed, like XOA) to `path` on exit."""
    with tempfile.TemporaryDirectory() as store_dir:
        yield zarr.open_group(store_dir, mode="w", zarr_format=2)
        path = Path(path)
        if path.exists():
            path.unlink()
        with zipfile.ZipFile(path, "w", zipfile.ZIP_STORED) as zf:
            for file_path in sorted(Path(store_dir).rglob("*")):
                if file_path.is_file():
                    zf.write(file_path, arcname=file_path.relative_to(store_dir))


def _create_array(group, name, data, attrs=None):
    """Create a zarr array and populate it, without relying on create_array's
    data= kwarg: not present in zarr 3.0.x (our pinned version), only added in
    later 3.x releases.
    """
    data = np.asarray(data)
    arr = group.create_array(
        name, shape=data.shape, dtype=data.dtype, compressors=ZARR_COMPRESSOR
    )
    if data.size:
        arr[...] = data
    if attrs:
        arr.attrs.update(dict(attrs))
    return arr


def _label_pixel_counts(mask_full, mask_crop, n_labels):
    """Per label row (label L -> row L-1): its pixel count in total and inside the crop."""
    total = np.bincount(mask_full.ravel(), minlength=n_labels + 1)[1 : n_labels + 1]
    inside = np.bincount(mask_crop.ravel(), minlength=n_labels + 1)[1 : n_labels + 1]
    return total, inside


def _relabel(mask, keep_rows, n_labels):
    """Renumber a label raster to 1..N following `keep_rows` order; other labels become 0."""
    lut = np.zeros(n_labels + 1, dtype=mask.dtype)
    lut[keep_rows + 1] = np.arange(1, len(keep_rows) + 1, dtype=mask.dtype)
    return lut[mask]


def _shift_vertices(vertices, origin_x_um, origin_y_um):
    # Rows are flattened (x0, y0, x1, y1, ...) pairs, padded by repeating the
    # last vertex, so every even/odd column is a valid x/y coordinate.
    shifted = vertices.copy()
    shifted[:, 0::2] -= origin_x_um
    shifted[:, 1::2] -= origin_y_um
    return shifted


def crop_cells_zarr(src_zip, dst_zip, window_px, origin_x_um, origin_y_um):
    """Crop cells.zarr.zip and return the (sorted) row indices of the kept cells.

    Only cells lying *entirely* inside the crop window, with all of their
    nuclei, are kept: cells cut by the window edge would otherwise end up with
    truncated masks and polygons, and transcripts/nuclei outside the bundle.
    Both label rasters are renumbered to a contiguous 1..N range matching the
    row order of the cropped tables (label L -> row L-1), and the polygon sets
    are subset to the same cells/nuclei, since xeniumranger derives its imported
    cell/nucleus counts from them and requires them to match the rasters.
    """
    px_x0, px_x1, px_y0, px_y1 = window_px
    with read_zarr_zip(src_zip) as src:
        # polygon_set_names is ["nucleus", "cell"]: index 0 holds the nuclei
        # (masks/0, polygon_sets/0), index 1 the cells (masks/1, polygon_sets/1).
        nucleus_mask = src["masks"]["0"][...]
        cell_mask = src["masks"]["1"][...]
        nucleus_mask_crop = nucleus_mask[px_y0:px_y1, px_x0:px_x1]
        cell_mask_crop = cell_mask[px_y0:px_y1, px_x0:px_x1]

        n_cells = src["cell_id"].shape[0]
        nucleus_cell_index = src["polygon_sets"]["0"]["cell_index"][...]
        n_nuclei = len(nucleus_cell_index)

        cell_total, cell_in_crop = _label_pixel_counts(
            cell_mask, cell_mask_crop, n_cells
        )
        nucleus_total, nucleus_in_crop = _label_pixel_counts(
            nucleus_mask, nucleus_mask_crop, n_nuclei
        )
        cell_inside = (cell_total > 0) & (cell_in_crop == cell_total)
        cells_with_cut_nucleus = np.unique(
            nucleus_cell_index[nucleus_in_crop < nucleus_total]
        )
        cell_inside[cells_with_cut_nucleus] = False

        kept_cells = np.nonzero(cell_inside)[0]
        kept_nuclei = np.nonzero(np.isin(nucleus_cell_index, kept_cells))[0]
        new_cell_index = np.full(n_cells, -1, dtype=np.int64)
        new_cell_index[kept_cells] = np.arange(len(kept_cells))

        cell_summary = src["cell_summary"][...][kept_cells].copy()
        column_names = list(src["cell_summary"].attrs["column_names"])
        # rebase centroid columns to the crop's own origin, same as every other
        # micron-valued coordinate (see the note in main() for why).
        for col, offset in [
            ("cell_centroid_x", origin_x_um),
            ("nucleus_centroid_x", origin_x_um),
            ("cell_centroid_y", origin_y_um),
            ("nucleus_centroid_y", origin_y_um),
        ]:
            cell_summary[:, column_names.index(col)] -= offset

        with write_zarr_zip(dst_zip) as dst:
            dst.attrs.update(dict(src.attrs))
            dst.attrs["number_cells"] = int(len(kept_cells))

            masks_grp = dst.create_group("masks")
            _create_array(
                masks_grp, "0", _relabel(nucleus_mask_crop, kept_nuclei, n_nuclei)
            )
            _create_array(masks_grp, "1", _relabel(cell_mask_crop, kept_cells, n_cells))
            _create_array(
                masks_grp,
                "homogeneous_transform",
                src["masks"]["homogeneous_transform"][...],
            )

            _create_array(dst, "cell_id", src["cell_id"][...][kept_cells])
            _create_array(dst, "cell_summary", cell_summary, src["cell_summary"].attrs)

            dst_polygon_sets = dst.create_group("polygon_sets")
            for key, keep_rows in [("0", kept_nuclei), ("1", kept_cells)]:
                src_set = src["polygon_sets"][key]
                dst_set = dst_polygon_sets.create_group(key)
                dst_set.attrs.update(dict(src_set.attrs))
                for arr_name, arr in src_set.arrays():
                    data = arr[...][keep_rows]
                    if arr_name == "cell_index":
                        data = new_cell_index[data].astype(arr.dtype)
                    elif arr_name == "vertices":
                        data = _shift_vertices(data, origin_x_um, origin_y_um)
                    _create_array(dst_set, arr_name, data, arr.attrs)

        logger.info(
            f"cells.zarr.zip: {len(kept_cells)} of {n_cells} cells and "
            f"{len(kept_nuclei)} of {n_nuclei} nuclei lie entirely inside the crop"
        )
        cell_id_str = cell_id_str_from_prefix_suffix(
            src["cell_id"][...][kept_cells, 0], src["cell_id"][...][kept_cells, 1]
        )
    return kept_cells, cell_id_str


def crop_cell_feature_matrix(src_path, dst_path, keep_cell_ids):
    """Reslice a CellRanger-format HDF5 feature-barcode matrix to the kept barcodes."""
    with h5py.File(src_path, "r") as src:
        barcodes = src["matrix/barcodes"][...]
        barcodes_str = np.array([b.decode() for b in barcodes])
        keep_idx = np.nonzero(np.isin(barcodes_str, list(keep_cell_ids)))[0]

        shape = src["matrix/shape"][...]
        n_features = int(shape[0])
        n_barcodes = int(shape[1])
        data = src["matrix/data"][...]
        indices = src["matrix/indices"][...]
        indptr = src["matrix/indptr"][...]

        mat = csc_matrix((data, indices, indptr), shape=(n_features, n_barcodes))
        mat_kept = mat[:, keep_idx]
        mat_kept.sort_indices()

        with h5py.File(dst_path, "w") as dst:
            for key, val in src.attrs.items():
                dst.attrs[key] = val
            grp = dst.create_group("matrix")
            for key, val in src["matrix"].attrs.items():
                grp.attrs[key] = val
            grp.create_dataset("data", data=mat_kept.data.astype(data.dtype))
            grp.create_dataset("indices", data=mat_kept.indices.astype(indices.dtype))
            grp.create_dataset("indptr", data=mat_kept.indptr.astype(indptr.dtype))
            grp.create_dataset(
                "shape", data=np.array([n_features, len(keep_idx)], dtype=shape.dtype)
            )
            grp.create_dataset("barcodes", data=barcodes[keep_idx])

            feat_grp = grp.create_group("features")
            src_feat = src["matrix/features"]
            for key in src_feat.keys():
                feat_grp.create_dataset(key, data=src_feat[key][...])
            for key, val in src_feat.attrs.items():
                feat_grp.attrs[key] = val


def crop_cell_feature_matrix_zarr(src_zip, dst_zip, kept_cells):
    """Subset cell_feature_matrix.zarr.zip (Xenium Explorer's copy of the matrix)
    to the kept cells: stored both feature-major (CSR, top level) and cell-major
    (`csc/`), with cells in the same row order as cells.zarr.zip.
    """
    with read_zarr_zip(src_zip) as src, write_zarr_zip(dst_zip) as dst:
        src_cf = src["cell_features"]
        n_features = int(src_cf.attrs["number_features"])
        n_cells = int(src_cf.attrs["number_cells"])
        mat = csr_matrix(
            (src_cf["data"][...], src_cf["indices"][...], src_cf["indptr"][...]),
            shape=(n_features, n_cells),
        )[:, kept_cells]
        by_feature = mat.tocsr()
        by_feature.sort_indices()
        by_cell = mat.tocsc()
        by_cell.sort_indices()

        dst_cf = dst.create_group("cell_features")
        dst_cf.attrs.update(dict(src_cf.attrs))
        dst_cf.attrs["number_cells"] = int(len(kept_cells))
        _create_array(dst_cf, "cell_id", src_cf["cell_id"][...][kept_cells])
        for dst_grp, src_grp, sparse in [
            (dst_cf, src_cf, by_feature),
            (dst_cf.create_group("csc"), src_cf["csc"], by_cell),
        ]:
            for name in ["data", "indices", "indptr"]:
                _create_array(
                    dst_grp, name, getattr(sparse, name).astype(src_grp[name].dtype)
                )


def crop_analysis_zarr(src_zip, dst_zip, kept_cells, n_cells):
    """Subset analysis.zarr.zip's cell groupings (clusterings: per grouping, a
    CSR of group -> cell row indices) to the kept cells, renumbering the cell
    indices to the cropped row order. Groups may end up empty; they're kept so
    group_names stay aligned.
    """
    new_cell_index = np.full(n_cells, -1, dtype=np.int64)
    new_cell_index[kept_cells] = np.arange(len(kept_cells))
    with read_zarr_zip(src_zip) as src, write_zarr_zip(dst_zip) as dst:
        dst.attrs.update(dict(src.attrs))
        src_groups = src["cell_groups"]
        dst_groups = dst.create_group("cell_groups")
        dst_groups.attrs.update(dict(src_groups.attrs))
        for key, grouping in src_groups.groups():
            indices = grouping["indices"][...]
            indptr = grouping["indptr"][...]
            new_indices, new_indptr = [], [0]
            for start, end in zip(indptr[:-1], indptr[1:]):
                members = new_cell_index[indices[start:end]]
                members = members[members >= 0]
                new_indices.append(members)
                new_indptr.append(new_indptr[-1] + len(members))
            dst_grouping = dst_groups.create_group(key)
            _create_array(
                dst_grouping,
                "indices",
                np.concatenate(new_indices).astype(indices.dtype),
            )
            _create_array(
                dst_grouping, "indptr", np.array(new_indptr, dtype=indptr.dtype)
            )


def _tile_layout(location, gene_identity, is_high_quality, grid_size, n_genes):
    """Split transcripts (or clusters) over the square tiles of one grid level.

    Returns a list of (tile key, row order, gene_offset) tuples. Within a tile,
    XOA stores all high quality rows first, then all low quality rows, each
    sorted by gene; gene_offset gives, per gene, the [start, end) of its low
    quality rows followed by those of its high quality rows (0, 0 when empty).
    """
    tile_x = np.floor(location[:, 0] / grid_size).astype(np.int64)
    tile_y = np.floor(location[:, 1] / grid_size).astype(np.int64)
    tiles = []
    for tx, ty in sorted(set(zip(tile_x.tolist(), tile_y.tolist()))):
        rows = np.nonzero((tile_x == tx) & (tile_y == ty))[0]
        rows = rows[np.lexsort((gene_identity[rows], ~is_high_quality[rows]))]
        gene_offset = np.zeros((n_genes, 4), dtype=np.uint32)
        n_high = int(is_high_quality[rows].sum())
        for block_start, block, cols in [
            (0, rows[:n_high], [2, 3]),
            (n_high, rows[n_high:], [0, 1]),
        ]:
            counts = np.bincount(gene_identity[block], minlength=n_genes)
            ends = block_start + np.cumsum(counts)
            starts = ends - counts
            present = counts > 0
            gene_offset[present, cols[0]] = starts[present]
            gene_offset[present, cols[1]] = ends[present]
        tiles.append((f"{tx},{ty}", rows, gene_offset))
    return tiles


def _density_csr(location, identity, n_ids, origin, grid_size, rows, cols):
    """Count transcripts per (identity, grid row, grid col), as a CSR with
    one row per (identity, grid row) pair, the layout of transcripts.zarr.zip's
    density/gene and density/codeword groups.
    """
    grid_col = ((location[:, 0] - origin[0]) // grid_size[0]).astype(np.int64)
    grid_row = ((location[:, 1] - origin[1]) // grid_size[1]).astype(np.int64)
    mat = coo_matrix(
        (
            np.ones(len(location), dtype=np.int64),
            (identity.astype(np.int64) * rows + grid_row, grid_col),
        ),
        shape=(n_ids * rows, cols),
    ).tocsr()
    mat.sum_duplicates()
    mat.sort_indices()
    return mat


def crop_transcripts_zarr(src_zip, dst_zip, keep_transcript_ids, origin_um, window_um):
    """Rebuild transcripts.zarr.zip from the transcripts kept in transcripts.parquet.

    The store holds a multi-resolution grid of square tiles: level 0 has one
    row per transcript (every field kept verbatim, bar the rebased location),
    coarser levels hold per-gene clusters of nearby transcripts for zoomed-out
    rendering in Xenium Explorer. Since rebasing the coordinates moves
    transcripts between tiles, every level is recomputed from the kept level 0
    rows, as are the high quality transcript density grids. The coarse
    metrics_density grid (derived QC metrics, not recomputable from the
    transcripts alone) is sliced to the bins overlapping the crop instead.
    """
    origin_x_um, origin_y_um = origin_um
    with read_zarr_zip(src_zip) as src:
        src_grids = src["grids"]
        grid_size = float(src_grids.attrs["grid_size"][0])
        n_levels = int(src_grids.attrs["number_levels"])
        n_genes = int(src.attrs["number_genes"])
        n_codewords = int(src.attrs["codeword_count"])

        level0_keys = src_grids.attrs["grid_keys"][0]
        fields = [
            name
            for name, _ in src_grids["0"][level0_keys[0]].arrays()
            if name != "gene_offset"
        ]
        level0 = {
            name: np.concatenate(
                [src_grids["0"][key][name][...] for key in level0_keys]
            )
            for name in fields
        }
        uuid = level0["uuid"].astype(np.uint64)
        transcript_id = (uuid[:, 1] << np.uint64(32)) | uuid[:, 0]
        keep = np.isin(transcript_id, keep_transcript_ids)
        level0 = {name: values[keep] for name, values in level0.items()}
        location = level0["location"].copy()
        location[:, 0] -= origin_x_um
        location[:, 1] -= origin_y_um
        level0["location"] = location
        gene_identity = level0["gene_identity"][:, 0]
        is_high_quality = level0["quality_score"][:, 0] >= HIGH_QV_THRESHOLD
        codeword = level0["codeword_identity"][:, 0]
        logger.info(f"transcripts.zarr.zip: {int(keep.sum())} of {len(keep)} kept")

        with write_zarr_zip(dst_zip) as dst:
            dst.attrs.update(dict(src.attrs))
            dst.attrs["number_rnas"] = int(keep.sum())

            grids_attrs = dict(src_grids.attrs)
            grids_attrs["codeword_to_transcript_counts"] = np.bincount(
                codeword, minlength=n_codewords
            ).tolist()
            grid_keys, grid_number_objects, objects_per_tile_per_gene = [], [], []
            dst_grids = dst.create_group("grids")
            for level in range(n_levels):
                if level == 0:
                    level_data = level0
                    level_gene = gene_identity
                    level_high_quality = is_high_quality
                else:
                    # Merge transcripts of the same gene and quality class
                    # within a bin that doubles in size with every level, as
                    # an approximation of XOA's own (undocumented) clustering
                    # (bin sizes chosen to give similar cluster counts).
                    bin_um = 4.0 * 2**level
                    group_keys = np.column_stack(
                        [
                            gene_identity,
                            is_high_quality,
                            np.floor(location[:, 0] / bin_um),
                            np.floor(location[:, 1] / bin_um),
                        ]
                    )
                    cluster_keys, cluster, cluster_count = np.unique(
                        group_keys, axis=0, return_inverse=True, return_counts=True
                    )
                    cluster = cluster.ravel()
                    n_clusters = len(cluster_count)
                    cluster_location = (
                        np.column_stack(
                            [
                                np.bincount(cluster, location[:, dim], n_clusters)
                                for dim in range(3)
                            ]
                        )
                        / cluster_count[:, None]
                    )
                    level_gene = cluster_keys[:, 0].astype(gene_identity.dtype)
                    level_high_quality = cluster_keys[:, 1].astype(bool)
                    level_data = {
                        "location": cluster_location.astype(np.float32),
                        "cluster_count": cluster_count.astype(np.uint32)[:, None],
                        "gene_identity": level_gene.astype(np.uint16)[:, None],
                    }

                tiles = _tile_layout(
                    level_data["location"],
                    level_gene,
                    level_high_quality,
                    grid_size * 2**level,
                    n_genes,
                )
                # every tile array carries column_names/column_descriptions
                # attributes, which xeniumranger reads the columns by: reuse
                # those of the source's first tile at the same level.
                src_tile = src_grids[str(level)][src_grids.attrs["grid_keys"][level][0]]
                dst_level = dst_grids.create_group(str(level))
                high_counts, low_counts = [], []
                for key, rows, gene_offset in tiles:
                    dst_tile = dst_level.create_group(key)
                    for name, values in level_data.items():
                        _create_array(
                            dst_tile, name, values[rows], src_tile[name].attrs
                        )
                    _create_array(
                        dst_tile,
                        "gene_offset",
                        gene_offset,
                        src_tile["gene_offset"].attrs,
                    )
                    high_counts.append(gene_offset[:, 3] - gene_offset[:, 2])
                    low_counts.append(gene_offset[:, 1] - gene_offset[:, 0])
                grid_keys.append([key for key, _, _ in tiles])
                grid_number_objects.append([len(rows) for _, rows, _ in tiles])
                objects_per_tile_per_gene.append(
                    {
                        "high_qscore": np.max(high_counts, axis=0).tolist(),
                        "low_qscore": np.max(low_counts, axis=0).tolist(),
                    }
                )
            grids_attrs["grid_keys"] = grid_keys
            grids_attrs["grid_number_objects"] = grid_number_objects
            grids_attrs["grid_array_shapes"] = [
                [{} for _ in keys] for keys in grid_keys
            ]
            grids_attrs["number_objects_per_tile_per_gene"] = objects_per_tile_per_gene
            dst_grids.attrs.update(grids_attrs)

            dst_density = dst.create_group("density")
            for name, identity, n_ids in [
                ("gene", gene_identity, n_genes),
                ("codeword", codeword, n_codewords),
            ]:
                src_density = src["density"][name]
                density_attrs = dict(src_density.attrs)
                density_grid = [float(g) for g in density_attrs["grid_size"]]
                origin = [
                    float(np.floor(location[:, dim].min() / density_grid[dim]))
                    * density_grid[dim]
                    for dim in range(2)
                ]
                cols, rows = [
                    int((location[:, dim].max() - origin[dim]) // density_grid[dim]) + 1
                    for dim in range(2)
                ]
                # Only high quality transcripts are counted, and (for
                # codewords) only a subset of identities is covered, by an
                # undocumented rule: keep the source's own selection, i.e.
                # every identity with a non-empty density row in the source.
                src_rows_per_id = int(src_density.attrs["rows"])
                covered = (
                    np.diff(src_density["indptr"][...])
                    .reshape(n_ids, src_rows_per_id)
                    .sum(axis=1)
                    > 0
                )
                counted = is_high_quality & covered[identity]
                density = _density_csr(
                    location[counted],
                    identity[counted],
                    n_ids,
                    origin,
                    density_grid,
                    rows,
                    cols,
                )
                density_attrs.update(
                    {
                        "origin": {"x": origin[0], "y": origin[1]},
                        "rows": rows,
                        "cols": cols,
                    }
                )
                dst_grp = dst_density.create_group(name)
                dst_grp.attrs.update(density_attrs)
                for arr_name in ["data", "indices", "indptr"]:
                    _create_array(
                        dst_grp,
                        arr_name,
                        getattr(density, arr_name).astype(src_density[arr_name].dtype),
                    )

            metrics_density = src["metrics_density"][...]
            slices = []
            for axis, dim, (window_start, window_end), crop_origin in [
                (1, "x", window_um[0], origin_x_um),
                (0, "y", window_um[1], origin_y_um),
            ]:
                grid_origin = src.attrs[f"metrics_density_{dim}_origin"]
                spacing = src.attrs[f"metrics_density_{dim}_spacing"]
                count = metrics_density.shape[axis]
                start = int(np.clip((window_start - grid_origin) // spacing, 0, count))
                end = int(
                    np.clip(np.ceil((window_end - grid_origin) / spacing), start, count)
                )
                slices.append(slice(start, end))
                dst.attrs[f"metrics_density_{dim}_origin"] = float(
                    grid_origin + start * spacing - crop_origin
                )
                dst.attrs[f"metrics_density_{dim}_count"] = end - start
            _create_array(
                dst,
                "metrics_density",
                metrics_density[slices[1], slices[0]],
                src["metrics_density"].attrs,
            )

            for name, arr in src.arrays():
                if name != "metrics_density":
                    _create_array(dst, name, arr[...], arr.attrs)


def crop_ome_tiff(src_path, dst_path, window_px):
    """Crop every plane of a (possibly multi-plane, pyramidal) XOA OME-TIFF.

    The layout of the source is preserved: tiled, same compression, and the
    same number of 2x-downsampled sub-resolutions per plane, stored as SubIFDs.
    The reader (for multi-channel datasets) requires a valid OME-XML
    ImageDescription with named channels, and reconstructs the multi-file
    series from it; a plain tifffile.imwrite(..., ome=True) auto-generates
    unnamed channels and breaks that. So the file's own original OME-XML
    (identical across the channel files bar its own document UUID) is reused
    verbatim, with only SizeX/SizeY patched to the cropped extent.
    """
    px_x0, px_x1, px_y0, px_y1 = window_px
    ome_xml = tifffile.tiffcomment(src_path)
    ome_xml = re.sub(r'SizeX="\d+"', f'SizeX="{px_x1 - px_x0}"', ome_xml)
    ome_xml = re.sub(r'SizeY="\d+"', f'SizeY="{px_y1 - px_y0}"', ome_xml)
    with tifffile.TiffFile(src_path) as tif, tifffile.TiffWriter(dst_path) as writer:
        first_page = tif.pages[0]
        n_sublevels = len(first_page.subifds or ())
        options = {
            "photometric": "minisblack",
            "compression": first_page.compression,
            "tile": (
                (first_page.tilelength, first_page.tilewidth)
                if first_page.is_tiled
                else None
            ),
            "metadata": None,
        }
        # Iterate the file's own top-level pages (one per plane), not its OME
        # series: the OME metadata references sibling channel files, which
        # would otherwise make tifffile stitch them into one (C, Y, X) array.
        for page_idx, page in enumerate(tif.pages):
            plane = page.asarray()[px_y0:px_y1, px_x0:px_x1]
            writer.write(
                plane,
                subifds=n_sublevels,
                description=ome_xml.encode("utf-8") if page_idx == 0 else None,
                **options,
            )
            for _ in range(n_sublevels):
                plane = plane[::2, ::2]
                writer.write(plane, subfiletype=1, **options)
        logger.info(
            f"{Path(src_path).name}: {len(tif.pages)} x {first_page.shape} -> "
            f"{(px_y1 - px_y0, px_x1 - px_x0)}"
        )


def main(par):
    input_dir = Path(par["input"])
    output_dir = Path(par["output"])
    if output_dir.exists():
        # Clear contents rather than removing and recreating the directory
        # itself: output_dir may be a Docker bind-mount point (e.g. when
        # re-running against an existing, possibly empty, output path), which
        # can't be rmdir'd from inside the container.
        for child in output_dir.iterdir():
            if child.is_dir():
                shutil.rmtree(child)
            else:
                child.unlink()
    else:
        output_dir.mkdir(parents=True)

    logger.info("Reading cells.parquet")
    cells = pd.read_parquet(input_dir / "cells.parquet")

    logger.info("Clustering cell centroids into patches")
    coords_um = cells[["x_centroid", "y_centroid"]].to_numpy()
    labels = cluster_by_distance(coords_um, eps=par["eps_um"])
    all_ids, all_counts = np.unique(labels, return_counts=True)
    keep = all_counts >= par["min_cells_per_patch"]
    patch_ids, counts = all_ids[keep], all_counts[keep]
    order = np.argsort(-counts, kind="stable")
    patch_ids, counts = patch_ids[order], counts[order]
    for patch_id, count in zip(patch_ids, counts):
        logger.info(f"patch {patch_id}: {count} cells")
    n_noise = int(all_counts[~keep].sum())
    if n_noise:
        logger.info(
            f"{n_noise} cells in clusters smaller than min_cells_per_patch (excluded)"
        )

    if par["patch_rank"] >= len(patch_ids):
        raise ValueError(
            f"patch_rank {par['patch_rank']} out of range: only "
            f"{len(patch_ids)} patches found"
        )
    chosen_patch = patch_ids[par["patch_rank"]]
    patch_coords = coords_um[labels == chosen_patch]
    logger.info(f"Selected patch {chosen_patch} ({len(patch_coords)} cells)")

    margin_um = par["margin_um"]
    x_min = patch_coords[:, 0].min() - margin_um
    x_max = patch_coords[:, 0].max() + margin_um
    y_min = patch_coords[:, 1].min() - margin_um
    y_max = patch_coords[:, 1].max() + margin_um
    logger.info(f"bbox (um): x=[{x_min:.1f},{x_max:.1f}] y=[{y_min:.1f},{y_max:.1f}]")

    with open(input_dir / "experiment.xenium") as f:
        experiment = json.load(f)
    pixel_size = experiment["pixel_size"]

    focus_paths = sorted((input_dir / "morphology_focus").glob("*.ome.tif"))
    with tifffile.TiffFile(focus_paths[0]) as tif:
        image_height, image_width = tif.pages[0].shape
    px_x0 = max(int(x_min / pixel_size), 0)
    px_x1 = min(int(np.ceil(x_max / pixel_size)), image_width)
    px_y0 = max(int(y_min / pixel_size), 0)
    px_y1 = min(int(np.ceil(y_max / pixel_size)), image_height)
    logger.info(f"bbox (px): x=[{px_x0},{px_x1}] y=[{px_y0},{px_y1}]")
    window_px = (px_x0, px_x1, px_y0, px_y1)
    # The crop window in (source) microns, matching the cropped pixel extent.
    window_um = (
        (px_x0 * pixel_size, px_x1 * pixel_size),
        (px_y0 * pixel_size, px_y1 * pixel_size),
    )

    # The cropped images/masks are new arrays starting at pixel (0, 0); every
    # micron-valued coordinate (cell/nucleus centroids, boundary vertices,
    # transcript locations) must be rebased by the same origin so they stay
    # aligned with the cropped pixel data, since the raw bundle format has no
    # separate offset/translation field to record a crop's origin.
    origin_x_um = window_um[0][0]
    origin_y_um = window_um[1][0]

    logger.info("Cropping cells.zarr.zip")
    kept_cells, kept_cell_ids = crop_cells_zarr(
        input_dir / "cells.zarr.zip",
        output_dir / "cells.zarr.zip",
        window_px,
        origin_x_um,
        origin_y_um,
    )
    keep_cell_ids = set(kept_cell_ids.tolist())

    tables = {}
    cells_kept = cells[cells["cell_id"].isin(keep_cell_ids)].copy()
    cells_kept["x_centroid"] -= origin_x_um
    cells_kept["y_centroid"] -= origin_y_um
    tables["cells"] = cells_kept
    logger.info(f"cells.parquet: {len(cells_kept)} rows kept (of {len(cells)})")

    for name in ["cell_boundaries", "nucleus_boundaries"]:
        df = pd.read_parquet(input_dir / f"{name}.parquet")
        df_kept = df[df["cell_id"].isin(keep_cell_ids)].copy()
        df_kept["vertex_x"] -= origin_x_um
        df_kept["vertex_y"] -= origin_y_um
        tables[name] = df_kept
        logger.info(f"{name}.parquet: {len(df_kept)} rows kept (of {len(df)})")

    logger.info("Cropping transcripts.parquet")
    transcripts = pd.read_parquet(input_dir / "transcripts.parquet")
    tmask = (
        (transcripts["x_location"] >= window_um[0][0])
        & (transcripts["x_location"] < window_um[0][1])
        & (transcripts["y_location"] >= window_um[1][0])
        & (transcripts["y_location"] < window_um[1][1])
    )
    transcripts_kept = transcripts[tmask].copy()
    # transcripts inside the window assigned to a cell cut by the window edge
    # (and so dropped) become unassigned.
    dropped_cell = ~transcripts_kept["cell_id"].isin(keep_cell_ids) & (
        transcripts_kept["cell_id"] != "UNASSIGNED"
    )
    transcripts_kept.loc[dropped_cell, "cell_id"] = "UNASSIGNED"
    transcripts_kept.loc[dropped_cell, "overlaps_nucleus"] = 0
    keep_transcript_ids = transcripts_kept["transcript_id"].to_numpy(dtype=np.uint64)
    transcripts_kept["x_location"] -= origin_x_um
    transcripts_kept["y_location"] -= origin_y_um
    tables["transcripts"] = transcripts_kept
    logger.info(
        f"{int(tmask.sum())} transcripts kept (of {len(transcripts)}), "
        f"{int(dropped_cell.sum())} unassigned from dropped edge cells"
    )

    # The parquet tables are what the converters read; the legacy *.csv.gz
    # duplicates are only written when the source bundle has them, so the
    # cropped bundle keeps the same file set as the source (xeniumranger copies
    # or regenerates each of them in its output bundle).
    for name, df in tables.items():
        df.to_parquet(output_dir / f"{name}.parquet", index=False)
        if (input_dir / f"{name}.csv.gz").exists():
            df.to_csv(output_dir / f"{name}.csv.gz", index=False, compression="gzip")

    logger.info("Cropping cell_feature_matrix.h5")
    crop_cell_feature_matrix(
        input_dir / "cell_feature_matrix.h5",
        output_dir / "cell_feature_matrix.h5",
        keep_cell_ids,
    )

    # Xenium Explorer companions, required by xeniumranger. Optional here, so
    # bundles without them (e.g. from older versions of this component) can
    # still be cropped for the converters.
    optional_zarr_stores = {
        "transcripts.zarr.zip": lambda src, dst: crop_transcripts_zarr(
            src, dst, keep_transcript_ids, (origin_x_um, origin_y_um), window_um
        ),
        "cell_feature_matrix.zarr.zip": lambda src, dst: crop_cell_feature_matrix_zarr(
            src, dst, kept_cells
        ),
        "analysis.zarr.zip": lambda src, dst: crop_analysis_zarr(
            src, dst, kept_cells, len(cells)
        ),
    }
    for fname, crop in optional_zarr_stores.items():
        if (input_dir / fname).exists():
            logger.info(f"Cropping {fname}")
            crop(input_dir / fname, output_dir / fname)
        else:
            logger.warning(f"{fname} not found in input bundle, skipping")

    logger.info("Cropping morphology images")
    (output_dir / "morphology_focus").mkdir()
    tif_paths = list(focus_paths)
    if (input_dir / "morphology.ome.tif").exists():
        tif_paths.append(input_dir / "morphology.ome.tif")
    for tif_path in tif_paths:
        crop_ome_tiff(tif_path, output_dir / tif_path.relative_to(input_dir), window_px)

    experiment["num_cells"] = int(len(kept_cells))
    with open(output_dir / "experiment.xenium", "w") as f:
        json.dump(experiment, f, indent=4)

    # Copied through unchanged: xeniumranger requires these to be present, and
    # regenerates the metrics and summary for its own output bundle. The
    # *.tar.gz archives (Xenium Explorer / aux outputs) are dropped.
    for fname in ["metrics_summary.csv", "gene_panel.json", "analysis_summary.html"]:
        if (input_dir / fname).exists():
            shutil.copy2(input_dir / fname, output_dir / fname)

    logger.info(f"Wrote cropped bundle to {output_dir}")


if __name__ == "__main__":
    main(par)
