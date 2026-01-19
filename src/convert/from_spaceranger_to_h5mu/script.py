from pathlib import Path
import mudata
import scanpy as sc
import sys
import pandas as pd

## VIASH START
par = {
    "input": "spaceranger_test",
    "uns_metrics": "metrics_cellranger",
    "uns_probe_set": "probe_set",
    "obsm_coordinates": "spatial",
    "output": "foo.h5mu",
    "min_genes": None,
    "min_counts": None,
    "output_compression": "gzip",
    "output_type": "filtered"
}
meta = {
    "resources_dir": "src/utils"
}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()


def retrieve_input_data(spaceranger_output_bundle, input_type="filtered"):
    # Expected folder structure (showing only relevant files):
    # ├── Spatial/
    # │   └── tissue_positions.csv
    # ├── filtered_feature_bc_matrix.h5 OR raw_feature_bc_matrix.h5
    # ├── metrics_summary.csv
    # └── probe_set.csv

    matrix_pattern = "**/filtered_feature_bc_matrix.h5" if input_type == "filtered" else "**/raw_feature_bc_matrix.h5"
    spaceranger_file_patterns = {
      "count_matrix": matrix_pattern,
      "metrics_summary": "**/metrics_summary.csv",
      "probe_set": "**/probe_set.csv",
      "spatial_coords": "**/Spatial/tissue_positions.csv"
    }

    spaceranger_output_bundle = Path(spaceranger_output_bundle)

    spaceranger_files = {}

    for key, pattern in spaceranger_file_patterns.items():
        file = list(spaceranger_output_bundle.glob(pattern))
        assert len(file) == 1, (
          f"Expected exactly one file for pattern '{pattern}', found {len(file)}."
        )
        spaceranger_files[key] = file[0]

    return spaceranger_files


def main():
    spaceranger_files = retrieve_input_data(par["input"], input_type=par["output_type"])

    logger.info("Reading count matrix...")
    adata = sc.read_10x_h5(spaceranger_files["count_matrix"], gex_only=False)

    # set the gene ids as var_names
    logger.info("Renaming var columns")
    adata.var = adata.var.rename_axis("gene_symbol").reset_index().set_index("gene_ids")

    if par["uns_metrics"]:
        logger.info("Reading metrics summary file...")
        metrics_summary = pd.read_csv(
            spaceranger_files["metrics_summary"], decimal=".", quotechar='"', thousands=","
        )

        logger.info("Storing metrics summary in .uns slot...")
        adata.uns[par["uns_metrics"]] = metrics_summary

    if par["uns_probe_set"]:
        logger.info("Reading probe set file...")
        def read_hash_metadata(path):
            meta = {}
            with open(path, "r", encoding="utf-8") as f:
                for i, line in enumerate(f):
                    if not line.startswith("#"):
                        break
                    line = line[1:].strip()
                    if "=" in line:
                        k, v = line.split("=", 1)
                        meta[k.strip()] = v.strip()
            return meta

        meta = read_hash_metadata(spaceranger_files["probe_set"])
        probe_set = pd.read_csv(
            spaceranger_files["probe_set"], comment="#"
            )

        logger.info("Storing probe set in .uns slot...")
        adata.uns[par["uns_probe_set"]] = probe_set
        adata.uns[par["uns_probe_set"] + "_meta"] = meta

    logger.info("Reading spatial coordinates...")
    spatial_coords = pd.read_csv(
        spaceranger_files["spatial_coords"], decimal=".", thousands=","
    )

    spatial_coords_aligned = spatial_coords.set_index("barcode").reindex(adata.obs_names)
    logger.info("Storing spatial coordinates in .obsm slot...")
    adata.obsm[par["obsm_coordinates"]] = spatial_coords_aligned[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy()

    # might perform basic filtering to get rid of some data
    # applicable when starting from the raw counts
    if par["min_genes"]:
        logger.info("Filtering with min_genes=%d", par["min_genes"])
        sc.pp.filter_cells(adata, min_genes=par["min_genes"])

    if par["min_counts"]:
        logger.info("Filtering with min_counts=%d", par["min_counts"])
        sc.pp.filter_cells(adata, min_counts=par["min_counts"])

    # generate output
    logger.info("Convert to mudata")
    mdata = mudata.MuData(adata)

    # override root .obs and .uns
    mdata.obs = adata.obs
    mdata.uns = adata.uns

    # write output
    logger.info("Writing %s", par["output"])
    mdata.write_h5mu(par["output"], compression=par["output_compression"])


if __name__ == "__main__":
    main()
