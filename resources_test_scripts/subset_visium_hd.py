"""Crop a Visium HD SpaceRanger `outs` to a tiny tissue-straddling region.

Picks a contiguous WxW window of 8um bins that mixes in/out of tissue (so raw
keeps more bins than filtered), then applies that fullres pixel bbox to every
bin size, column-subsetting the count matrices to the region's barcodes.
feature_slice.h5 is reduced to an attrs-only stub (the reader only needs its
metadata_json). Images are kept as-is (coords stay aligned via the unchanged
scalefactor)."""

import os
import shutil
import sys

import h5py
import numpy as np
import pandas as pd

SRC = sys.argv[1]  # extracted outs dir
OUT = sys.argv[2]  # output tiny_spaceranger dir
DATASET_ID = sys.argv[3]  # feature_slice prefix / dataset id
BINS = ["square_002um", "square_008um", "square_016um"]
H5_FILES = ["filtered_feature_bc_matrix.h5", "raw_feature_bc_matrix.h5"]
SPATIAL = [
    "tissue_positions.parquet",
    "scalefactors_json.json",
    "tissue_hires_image.png",
    "tissue_lowres_image.png",
]
TP_COLS = [
    "barcode",
    "in_tissue",
    "array_row",
    "array_col",
    "pxl_row_in_fullres",
    "pxl_col_in_fullres",
]
W = 22


def choose_bbox():
    tp = pd.read_parquet(
        f"{SRC}/binned_outputs/square_008um/spatial/tissue_positions.parquet"
    )
    R, C = int(tp.array_row.max()) + 1, int(tp.array_col.max()) + 1
    in_t = np.zeros((R, C), np.int64)
    exists = np.zeros((R, C), np.int64)
    in_t[tp.array_row.values, tp.array_col.values] = tp.in_tissue.values == 1
    exists[tp.array_row.values, tp.array_col.values] = 1
    It, Iv = in_t.cumsum(0).cumsum(1), exists.cumsum(0).cumsum(1)

    def winsum(integral, r0, c0):
        r1, c1 = r0 + W - 1, c0 + W - 1
        tot = integral[r1, c1]
        tot -= integral[r0 - 1, c1] if r0 > 0 else 0
        tot -= integral[r1, c0 - 1] if c0 > 0 else 0
        tot += integral[r0 - 1, c0 - 1] if r0 > 0 and c0 > 0 else 0
        return int(tot)

    best = None
    for r0 in range(0, R - W, 2):
        for c0 in range(0, C - W, 2):
            v = winsum(Iv, r0, c0)
            if v < W * W:
                continue
            t = winsum(It, r0, c0)
            frac = t / v
            if 0.45 <= frac <= 0.80 and t >= 150:
                score = abs(frac - 0.6)
                if best is None or score < best[0]:
                    best = (score, r0, c0)
    assert best is not None, "no mixed-tissue window found"
    _, r0, c0 = best
    blk = tp[
        (tp.array_row >= r0)
        & (tp.array_row < r0 + W)
        & (tp.array_col >= c0)
        & (tp.array_col < c0 + W)
    ]
    print(
        f"  8um window at row[{r0},{r0 + W}) col[{c0},{c0 + W}): "
        f"{len(blk)} bins, {int(blk.in_tissue.sum())} in-tissue"
    )
    return (
        blk.pxl_row_in_fullres.min(),
        blk.pxl_row_in_fullres.max(),
        blk.pxl_col_in_fullres.min(),
        blk.pxl_col_in_fullres.max(),
    )


def subset_h5(src, dst, keep):
    with h5py.File(src, "r") as f:
        bc = f["matrix/barcodes"][:]
        cols = np.flatnonzero(np.isin(bc, list(keep)))
        indptr = f["matrix/indptr"][:]
        spans = [(int(indptr[j]), int(indptr[j + 1])) for j in cols]
        points = (
            np.concatenate([np.arange(s, e) for s, e in spans])
            if spans
            else np.array([], dtype=indptr.dtype)
        )
        new_indptr = np.zeros(len(cols) + 1, indptr.dtype)
        if spans:
            new_indptr[1:] = np.cumsum([e - s for s, e in spans])
        n_feat = int(f["matrix/shape"][0])
        with h5py.File(dst, "w") as g:
            for k, v in f.attrs.items():
                g.attrs[k] = v
            m = g.create_group("matrix")
            for k, v in f["matrix"].attrs.items():
                m.attrs[k] = v
            m.create_dataset("barcodes", data=bc[cols], compression="gzip")
            m.create_dataset("data", data=f["matrix/data"][points], compression="gzip")
            m.create_dataset(
                "indices", data=f["matrix/indices"][points], compression="gzip"
            )
            m.create_dataset("indptr", data=new_indptr, compression="gzip")
            m.create_dataset("shape", data=np.array([n_feat, len(cols)], np.int32))
            f.copy(f["matrix/features"], m, name="features")
    return len(cols)


def main():
    if os.path.exists(OUT):
        shutil.rmtree(OUT)
    os.makedirs(OUT)
    # feature_slice stub at top level
    with (
        h5py.File(f"{SRC}/feature_slice.h5", "r") as s,
        h5py.File(f"{OUT}/{DATASET_ID}_feature_slice.h5", "w") as d,
    ):
        for k, v in s.attrs.items():
            d.attrs[k] = v

    print(">>> choosing crop region from 8um grid")
    rmin, rmax, cmin, cmax = choose_bbox()

    for b in BINS:
        sdir, odir = f"{SRC}/binned_outputs/{b}", f"{OUT}/binned_outputs/{b}"
        os.makedirs(f"{odir}/spatial")
        tp = pd.read_parquet(f"{sdir}/spatial/tissue_positions.parquet")
        assert list(tp.columns) == TP_COLS, list(tp.columns)
        region = tp[
            (tp.pxl_row_in_fullres >= rmin)
            & (tp.pxl_row_in_fullres <= rmax)
            & (tp.pxl_col_in_fullres >= cmin)
            & (tp.pxl_col_in_fullres <= cmax)
        ].reset_index(drop=True)
        region.to_parquet(f"{odir}/spatial/tissue_positions.parquet", index=False)
        keep = {x.encode() if isinstance(x, str) else x for x in region.barcode}
        nobs = {h: subset_h5(f"{sdir}/{h}", f"{odir}/{h}", keep) for h in H5_FILES}
        # images + scalefactors kept as-is (coords align via unchanged scalefactor)
        shutil.copy(f"{sdir}/spatial/scalefactors_json.json", f"{odir}/spatial/")
        shutil.copy(f"{sdir}/spatial/tissue_hires_image.png", f"{odir}/spatial/")
        shutil.copy(f"{sdir}/spatial/tissue_lowres_image.png", f"{odir}/spatial/")
        print(
            f"  {b}: region {len(region)} bins | "
            f"filtered={nobs['filtered_feature_bc_matrix.h5']} "
            f"raw={nobs['raw_feature_bc_matrix.h5']}"
        )
    os.system(f"du -sh {OUT}")


if __name__ == "__main__":
    main()
