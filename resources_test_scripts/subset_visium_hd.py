"""Crop a Visium HD Space Ranger ``outs`` directory down to a tiny test fixture.

A Visium HD slide bins the whole capture area into millions of 2 um squares, so
even a downsampled dataset is far too large for CI. This script keeps only a
small square region of the tissue:

1. On the 8 um grid, find a contiguous ``WINDOW`` x ``WINDOW`` block of bins that
   straddles the tissue edge -- i.e. it contains both in-tissue and background
   bins. Keeping a mix means the *raw* matrix (all bins) ends up larger than the
   *filtered* matrix (in-tissue bins only), which the downstream tests rely on.
2. Convert that block to a full-resolution pixel bounding box and apply the same
   box to every bin size (2 / 8 / 16 um), so the regions line up across bins.
3. For each bin size, keep only the bins inside the box: filter
   ``tissue_positions.parquet`` and column-subset the count matrices to those
   barcodes.

``feature_slice.h5`` is replaced by an attributes-only stub, because
``spatialdata_io.visium_hd`` only reads its ``metadata_json`` attribute. The
tissue images and scale factors are copied unchanged -- the bin coordinates
still map onto them via the (unchanged) scale factor.

Usage:
    subset_visium_hd.py --input OUTS_DIR --output FIXTURE_DIR --dataset-id ID
"""

import argparse
import os
import shutil

import h5py
import numpy as np
import pandas as pd

# Bin sizes produced by Space Ranger HD, smallest to largest.
BIN_DIRS = ["square_002um", "square_008um", "square_016um"]
# The two count matrices kept for each bin size.
MATRIX_FILES = ["filtered_feature_bc_matrix.h5", "raw_feature_bc_matrix.h5"]
# Columns of tissue_positions.parquet, used to sanity-check the input.
TISSUE_POSITION_COLUMNS = [
    "barcode",
    "in_tissue",
    "array_row",
    "array_col",
    "pxl_row_in_fullres",
    "pxl_col_in_fullres",
]
# Side length (in 8 um bins) of the square region we crop to.
WINDOW = 22


def find_crop_box(outs_dir):
    """Find a tissue-straddling window on the 8 um grid.

    Returns the full-resolution pixel bounding box (min/max row, min/max col) of
    a ``WINDOW`` x ``WINDOW`` block of 8 um bins that mixes in-tissue and
    background bins.
    """
    positions = pd.read_parquet(
        f"{outs_dir}/binned_outputs/square_008um/spatial/tissue_positions.parquet"
    )
    n_rows = int(positions.array_row.max()) + 1
    n_cols = int(positions.array_col.max()) + 1

    # Lay the bins out on their (row, col) grid: one grid marking in-tissue bins,
    # one marking which grid cells have a bin at all.
    in_tissue_grid = np.zeros((n_rows, n_cols), np.int64)
    bin_present_grid = np.zeros((n_rows, n_cols), np.int64)
    rows = positions.array_row.values
    cols = positions.array_col.values
    in_tissue_grid[rows, cols] = positions.in_tissue.values == 1
    bin_present_grid[rows, cols] = 1

    # Summed-area tables: cumulative sums let us total any rectangular window in
    # O(1) instead of re-summing its cells.
    in_tissue_cumsum = in_tissue_grid.cumsum(0).cumsum(1)
    bin_present_cumsum = bin_present_grid.cumsum(0).cumsum(1)

    def window_total(cumsum, top, left):
        """Sum of a WINDOW x WINDOW window with top-left corner at (top, left)."""
        bottom, right = top + WINDOW - 1, left + WINDOW - 1
        total = cumsum[bottom, right]
        total -= cumsum[top - 1, right] if top > 0 else 0
        total -= cumsum[bottom, left - 1] if left > 0 else 0
        total += cumsum[top - 1, left - 1] if top > 0 and left > 0 else 0
        return int(total)

    # Scan windows (stride 2 for speed) and pick the one whose in-tissue fraction
    # is closest to 0.6 -- a balanced mix of tissue and background.
    best = None  # (distance_from_0.6, top, left)
    for top in range(0, n_rows - WINDOW, 2):
        for left in range(0, n_cols - WINDOW, 2):
            n_bins = window_total(bin_present_cumsum, top, left)
            if n_bins < WINDOW * WINDOW:
                continue  # window must be fully populated with bins
            n_in_tissue = window_total(in_tissue_cumsum, top, left)
            tissue_fraction = n_in_tissue / n_bins
            if 0.45 <= tissue_fraction <= 0.80 and n_in_tissue >= 150:
                distance = abs(tissue_fraction - 0.6)
                if best is None or distance < best[0]:
                    best = (distance, top, left)
    assert best is not None, "no tissue-straddling window found"

    _, top, left = best
    window_bins = positions[
        (positions.array_row >= top)
        & (positions.array_row < top + WINDOW)
        & (positions.array_col >= left)
        & (positions.array_col < left + WINDOW)
    ]
    print(
        f"  8um window at row[{top},{top + WINDOW}) col[{left},{left + WINDOW}): "
        f"{len(window_bins)} bins, {int(window_bins.in_tissue.sum())} in-tissue"
    )
    return (
        window_bins.pxl_row_in_fullres.min(),
        window_bins.pxl_row_in_fullres.max(),
        window_bins.pxl_col_in_fullres.min(),
        window_bins.pxl_col_in_fullres.max(),
    )


def subset_matrix(src_path, dst_path, keep_barcodes):
    """Write a copy of a 10x feature-barcode matrix keeping only some barcodes.

    The matrix is stored as CSC (one column per barcode), so subsetting barcodes
    means selecting columns: gather each kept column's value range from ``indptr``
    and rebuild the ``data`` / ``indices`` / ``indptr`` arrays. Returns the number
    of barcodes kept.
    """
    with h5py.File(src_path, "r") as src:
        barcodes = src["matrix/barcodes"][:]
        kept_columns = np.flatnonzero(np.isin(barcodes, list(keep_barcodes)))
        indptr = src["matrix/indptr"][:]

        # Each kept column occupies indptr[col]:indptr[col + 1] in data/indices.
        column_ranges = [(int(indptr[c]), int(indptr[c + 1])) for c in kept_columns]
        if column_ranges:
            kept_value_indices = np.concatenate(
                [np.arange(start, end) for start, end in column_ranges]
            )
        else:
            kept_value_indices = np.array([], dtype=indptr.dtype)

        new_indptr = np.zeros(len(kept_columns) + 1, indptr.dtype)
        if column_ranges:
            new_indptr[1:] = np.cumsum([end - start for start, end in column_ranges])

        n_features = int(src["matrix/shape"][0])
        with h5py.File(dst_path, "w") as dst:
            for key, value in src.attrs.items():
                dst.attrs[key] = value
            matrix = dst.create_group("matrix")
            for key, value in src["matrix"].attrs.items():
                matrix.attrs[key] = value
            matrix.create_dataset(
                "barcodes", data=barcodes[kept_columns], compression="gzip"
            )
            matrix.create_dataset(
                "data", data=src["matrix/data"][kept_value_indices], compression="gzip"
            )
            matrix.create_dataset(
                "indices",
                data=src["matrix/indices"][kept_value_indices],
                compression="gzip",
            )
            matrix.create_dataset("indptr", data=new_indptr, compression="gzip")
            matrix.create_dataset(
                "shape", data=np.array([n_features, len(kept_columns)], np.int32)
            )
            src.copy(src["matrix/features"], matrix, name="features")
    return len(kept_columns)


def stub_feature_slice(src_path, dst_path):
    """Copy only the root attributes of feature_slice.h5 (drop the heavy data)."""
    with h5py.File(src_path, "r") as src, h5py.File(dst_path, "w") as dst:
        for key, value in src.attrs.items():
            dst.attrs[key] = value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", required=True, help="extracted Space Ranger outs dir"
    )
    parser.add_argument("--output", required=True, help="output tiny_spaceranger dir")
    parser.add_argument(
        "--dataset-id", required=True, help="feature_slice prefix / dataset id"
    )
    args = parser.parse_args()

    if os.path.exists(args.output):
        shutil.rmtree(args.output)
    os.makedirs(args.output)

    stub_feature_slice(
        f"{args.input}/feature_slice.h5",
        f"{args.output}/{args.dataset_id}_feature_slice.h5",
    )

    print(">>> choosing crop region from 8um grid")
    row_min, row_max, col_min, col_max = find_crop_box(args.input)

    for bin_dir in BIN_DIRS:
        src_dir = f"{args.input}/binned_outputs/{bin_dir}"
        out_dir = f"{args.output}/binned_outputs/{bin_dir}"
        os.makedirs(f"{out_dir}/spatial")

        positions = pd.read_parquet(f"{src_dir}/spatial/tissue_positions.parquet")
        assert list(positions.columns) == TISSUE_POSITION_COLUMNS, list(
            positions.columns
        )
        region = positions[
            (positions.pxl_row_in_fullres >= row_min)
            & (positions.pxl_row_in_fullres <= row_max)
            & (positions.pxl_col_in_fullres >= col_min)
            & (positions.pxl_col_in_fullres <= col_max)
        ].reset_index(drop=True)
        region.to_parquet(f"{out_dir}/spatial/tissue_positions.parquet", index=False)

        # Barcodes are stored as bytes in the matrix h5, so encode the strings.
        keep_barcodes = {
            bc.encode() if isinstance(bc, str) else bc for bc in region.barcode
        }
        n_kept = {
            matrix_file: subset_matrix(
                f"{src_dir}/{matrix_file}", f"{out_dir}/{matrix_file}", keep_barcodes
            )
            for matrix_file in MATRIX_FILES
        }

        # Images and scale factors are copied unchanged; bin coordinates still
        # map onto the full image via the unchanged scale factor.
        for spatial_file in (
            "scalefactors_json.json",
            "tissue_hires_image.png",
            "tissue_lowres_image.png",
        ):
            shutil.copy(f"{src_dir}/spatial/{spatial_file}", f"{out_dir}/spatial/")

        print(
            f"  {bin_dir}: region {len(region)} bins | "
            f"filtered={n_kept['filtered_feature_bc_matrix.h5']} "
            f"raw={n_kept['raw_feature_bc_matrix.h5']}"
        )


if __name__ == "__main__":
    main()
