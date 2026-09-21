import json
import re
import shutil
import sys
import tempfile
import zipfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import tifffile
import zarr
from scipy.sparse import coo_matrix, csc_matrix
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


def _create_array(group, name, data):
    """Create a zarr array and populate it, without relying on create_array's
    data= kwarg: not present in zarr 3.0.x (our pinned version), only added in
    later 3.x releases.
    """
    arr = group.create_array(name, shape=data.shape, dtype=data.dtype)
    arr[:] = data
    return arr


def crop_cells_zarr(
    src_zip, dst_zip, px_x0, px_x1, px_y0, px_y1, origin_x_um, origin_y_um
):
    """Crop cells.zarr.zip and return the set of kept cell IDs.

    masks/1 (the cell label raster) is cropped first; the surviving nonzero
    label values are the canonical "kept cells" (label L -> row L-1 of
    cell_id/cell_summary), which keeps the output internally consistent, since
    spatialdata_io cross-checks the cropped labels against the cropped table.
    """
    with tempfile.TemporaryDirectory() as extract_dir:
        with zipfile.ZipFile(src_zip) as zf:
            zf.extractall(extract_dir)
        src = zarr.open(extract_dir, mode="r")

        masks_0 = src["masks"]["0"][px_y0:px_y1, px_x0:px_x1]
        masks_1 = src["masks"]["1"][px_y0:px_y1, px_x0:px_x1]
        transform = src["masks"]["homogeneous_transform"][...]

        kept_label_ids = np.unique(masks_1)
        kept_label_ids = kept_label_ids[kept_label_ids > 0]
        kept_row_idx = kept_label_ids - 1

        # Renumber masks_1 to a fresh, contiguous 1..N range matching the row
        # order of the cropped cell_id/cell_summary below (label L -> row L-1),
        # rather than leaving the original (now sparse) global label values in
        # place. Without this, cropping an already-cropped bundle a second time
        # would look up cell_id/cell_summary rows using stale label values that
        # no longer fit the (smaller) cropped arrays.
        remap = np.zeros(int(masks_1.max()) + 1, dtype=masks_1.dtype)
        remap[kept_label_ids] = np.arange(
            1, len(kept_label_ids) + 1, dtype=masks_1.dtype
        )
        masks_1 = remap[masks_1]

        cell_id_arr = src["cell_id"][...][kept_row_idx]
        cell_summary = src["cell_summary"][...][kept_row_idx].copy()
        column_names = list(src["cell_summary"].attrs["column_names"])
        column_descriptions = list(src["cell_summary"].attrs["column_descriptions"])
        # rebase centroid columns to the crop's own origin, same as every other
        # micron-valued coordinate (see the note in main() for why).
        for col, offset in [
            ("cell_centroid_x", origin_x_um),
            ("nucleus_centroid_x", origin_x_um),
            ("cell_centroid_y", origin_y_um),
            ("nucleus_centroid_y", origin_y_um),
        ]:
            cell_summary[:, column_names.index(col)] -= offset

        cell_id_str = cell_id_str_from_prefix_suffix(
            cell_id_arr[:, 0], cell_id_arr[:, 1]
        )

        with tempfile.TemporaryDirectory() as store_dir:
            dst = zarr.open(store_dir, mode="w")
            dst.attrs.update(dict(src.attrs))
            dst.attrs["number_cells"] = int(len(kept_row_idx))

            masks_grp = dst.create_group("masks")
            _create_array(masks_grp, "0", masks_0)
            _create_array(masks_grp, "1", masks_1)
            _create_array(masks_grp, "homogeneous_transform", transform)

            _create_array(dst, "cell_id", cell_id_arr)
            summary_arr = _create_array(dst, "cell_summary", cell_summary)
            summary_arr.attrs["column_names"] = column_names
            summary_arr.attrs["column_descriptions"] = column_descriptions

            # polygon_sets is unused by the reader beyond a self-consistency check
            # on its own (unmodified) length, so it's copied through untouched.
            src_polygon_sets = src["polygon_sets"]
            dst_polygon_sets = dst.create_group("polygon_sets")
            for key in src_polygon_sets.keys():
                sub_src = src_polygon_sets[key]
                sub_dst = dst_polygon_sets.create_group(key)
                for arr_name in sub_src.keys():
                    arr = _create_array(sub_dst, arr_name, sub_src[arr_name][...])
                    arr.attrs.update(dict(sub_src[arr_name].attrs))
                sub_dst.attrs.update(dict(sub_src.attrs))

            store_root = Path(store_dir)
            dst_zip = Path(dst_zip)
            if dst_zip.exists():
                dst_zip.unlink()
            with zipfile.ZipFile(dst_zip, "w", zipfile.ZIP_STORED) as zf:
                for file_path in sorted(store_root.rglob("*")):
                    if file_path.is_file():
                        zf.write(file_path, arcname=file_path.relative_to(store_root))

        return set(cell_id_str.tolist())


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

    px_x0, px_x1 = int(x_min / pixel_size), int(np.ceil(x_max / pixel_size))
    px_y0, px_y1 = int(y_min / pixel_size), int(np.ceil(y_max / pixel_size))
    logger.info(f"bbox (px): x=[{px_x0},{px_x1}] y=[{px_y0},{px_y1}]")

    # The cropped images/masks are new arrays starting at pixel (0, 0); every
    # micron-valued coordinate (cell/nucleus centroids, boundary vertices,
    # transcript locations) must be rebased by the same origin so they stay
    # aligned with the cropped pixel data, since the raw bundle format has no
    # separate offset/translation field to record a crop's origin.
    origin_x_um = px_x0 * pixel_size
    origin_y_um = px_y0 * pixel_size

    logger.info("Cropping cells.zarr.zip")
    keep_cell_ids = crop_cells_zarr(
        input_dir / "cells.zarr.zip",
        output_dir / "cells.zarr.zip",
        px_x0,
        px_x1,
        px_y0,
        px_y1,
        origin_x_um,
        origin_y_um,
    )
    logger.info(f"{len(keep_cell_ids)} cells kept after cropping cells.zarr.zip")

    cells_kept = cells[cells["cell_id"].isin(keep_cell_ids)].copy()
    cells_kept["x_centroid"] -= origin_x_um
    cells_kept["y_centroid"] -= origin_y_um
    cells_kept.to_parquet(output_dir / "cells.parquet", index=False)
    logger.info(f"cells.parquet: {len(cells_kept)} rows kept (of {len(cells)})")

    for fname in ["cell_boundaries.parquet", "nucleus_boundaries.parquet"]:
        df = pd.read_parquet(input_dir / fname)
        df_kept = df[df["cell_id"].isin(keep_cell_ids)].copy()
        df_kept["vertex_x"] -= origin_x_um
        df_kept["vertex_y"] -= origin_y_um
        df_kept.to_parquet(output_dir / fname, index=False)
        logger.info(f"{fname}: {len(df_kept)} rows kept (of {len(df)})")

    logger.info("Cropping transcripts.parquet")
    transcripts = pd.read_parquet(input_dir / "transcripts.parquet")
    tmask = (
        (transcripts["x_location"] >= x_min)
        & (transcripts["x_location"] <= x_max)
        & (transcripts["y_location"] >= y_min)
        & (transcripts["y_location"] <= y_max)
    )
    transcripts_kept = transcripts[tmask].copy()
    transcripts_kept["x_location"] -= origin_x_um
    transcripts_kept["y_location"] -= origin_y_um
    transcripts_kept.to_parquet(output_dir / "transcripts.parquet", index=False)
    logger.info(f"{int(tmask.sum())} transcripts kept (of {len(transcripts)})")

    logger.info("Cropping cell_feature_matrix.h5")
    crop_cell_feature_matrix(
        input_dir / "cell_feature_matrix.h5",
        output_dir / "cell_feature_matrix.h5",
        keep_cell_ids,
    )

    logger.info("Cropping morphology_focus images")
    morphology_dir = input_dir / "morphology_focus"
    (output_dir / "morphology_focus").mkdir()
    crop_h = px_y1 - px_y0
    crop_w = px_x1 - px_x0
    for tif_path in sorted(morphology_dir.glob("*.ome.tif")):
        # key=0: the OME metadata references sibling channel files, which would
        # otherwise make tifffile stitch them into one (C, Y, X) array; we want
        # just this file's own single physical (Y, X) page.
        image = tifffile.imread(tif_path, key=0)
        cropped = image[px_y0:px_y1, px_x0:px_x1]

        # The reader (for multi-channel datasets) requires a valid OME-XML
        # ImageDescription with named channels, and reconstructs the multi-file
        # series from it; a plain tifffile.imwrite(..., ome=True) auto-generates
        # unnamed channels and breaks that. So each file's own original OME-XML
        # (identical across the channel files bar its own document UUID) is
        # reused verbatim, with only SizeX/SizeY patched to the cropped extent.
        ome_xml = tifffile.tiffcomment(tif_path)
        ome_xml = re.sub(r'SizeX="\d+"', f'SizeX="{crop_w}"', ome_xml)
        ome_xml = re.sub(r'SizeY="\d+"', f'SizeY="{crop_h}"', ome_xml)
        tifffile.imwrite(
            output_dir / "morphology_focus" / tif_path.name,
            cropped,
            description=ome_xml.encode("utf-8"),
            metadata=None,
            photometric="minisblack",
        )
        logger.info(f"{tif_path.name}: {image.shape} -> {cropped.shape}")

    shutil.copy2(input_dir / "experiment.xenium", output_dir / "experiment.xenium")
    shutil.copy2(input_dir / "metrics_summary.csv", output_dir / "metrics_summary.csv")

    logger.info(f"Wrote cropped bundle to {output_dir}")


if __name__ == "__main__":
    main(par)
