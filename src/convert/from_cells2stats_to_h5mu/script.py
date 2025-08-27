import sys
from pathlib import Path
import scipy.sparse as sp
import pandas as pd
import mudata as mu
import anndata as ad
import re
import json
import zipfile
import os

## VIASH START
par = {
    "input": "./resources_test/aviti/aviti_teton_tiny_2",
    "modality": "rna",
    "output": "aviti_tiny_test.h5mu",
    "output_compression": "gzip",
    "layer_nuclear_counts": "nuclear_counts",
    "obsm_coordinates": "spatial",
    "obsm_cell_paint": "cell_paint",
    "obsm_cell_paint_nuclear": "cell_paint_nuclear",
    "obsm_cell_profiler": "cell_profiler",
    "obsm_unassigned_targets": "unassigned_targets",
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from unzip_archived_folder import extract_selected_files_from_zip


logger = setup_logger()


def assert_matching_order(var_names, count_columns, split_pattern=None):
    for var, col in zip(var_names, count_columns):
        count_var = col if not split_pattern else col.split("_Nuclear")[0]
        assert var == count_var, "Orders do not match"


def categorize_columns(column_list, target_panel):
    # Extract imaging and barcoding information from Panel.json
    imaging_batches = [tube["BatchName"] for tube in target_panel["ImagingPrimerTubes"]]
    barcoding_batches = [
        tube["BatchName"] for tube in target_panel["BarcodingPrimerTubes"]
    ]

    # Extract target information
    cellpaint_targets = [target["Target"] for target in target_panel["ImagingTargets"]]
    barcoding_targets = [
        target["Target"] for target in target_panel["BarcodingTargets"]
    ]

    # METADATA (for .obs and .obsm)
    # Fixed columns
    columns_fixed = [
        "Area",
        "AreaUm",
        "Cell",
        "NuclearArea",
        "NuclearAreaUm",
        "Tile",
        "Well",
        "WellLabel",
    ]
    obs_columns_fixed = list(set(columns_fixed) & set(column_list))

    # Coordinate columns
    coordinate_columns = ["X", "Y", "Xum", "Yum"]
    obsm_coordinate_columns = list(set(coordinate_columns) & set(column_list))

    # Cell Paint target intensity columns (format: {cell_paint_target.batch})
    cell_paint_columns = [
        col
        for col in column_list
        if any(
            col.startswith(f"{target}.") and col.endswith(f".{batch}")
            for target in cellpaint_targets
            for batch in imaging_batches
        )
    ]

    # Cell Paint nuclear target intensity columns (format: {cell_paint_target_Nuclear.batch})
    cell_paint_nuclear_columns = [
        col
        for col in column_list
        if any(
            col.startswith(f"{target}_Nuclear") and col.endswith(f".{batch}")
            for target in cellpaint_targets
            for batch in imaging_batches
        )
    ]

    # CellProfiler morphology metrics
    morphology_patterns = [
        r"^AreaShape_",
        r"^Granularity_",
        r"^Texture_",
        r"^Intensity_",
        r"^Location_",
        r"^RadialDistribution_",
    ]
    cell_profiler_columns = [
        col
        for col in column_list
        for pattern in morphology_patterns
        if re.match(pattern, col)
    ]

    # COUNT MATRICES (for .X and layers)
    # Feature Count Matrix - barcoding targets (format: {target.batch})
    # Includes cellular and nuclear counts
    count_columns = [
        col
        for col in column_list
        if any(
            col.startswith(f"{target}.") and col.endswith(f".{batch}")
            for target in barcoding_targets
            for batch in barcoding_batches
        )
    ]

    # Nuclear Feature Count Matrix - barcoding targets (format: {target_Nuclear.batch})
    # Includes only nuclear counts
    nuclear_count_columns = [
        col
        for col in column_list
        if any(
            col.startswith(f"{target}_Nuclear") and col.endswith(f".{batch}")
            for target in barcoding_targets
            for batch in barcoding_batches
        )
    ]

    # Unassigned columns (format: {Unassigned_*.*})
    unassigned_columns = [col for col in column_list if col.startswith("Unassigned")]

    # Make sure all columns have been categorized and have expected sizes
    assert len(count_columns) == len(nuclear_count_columns), (
        "Cellular and nuclear count columns do not match."
    )
    all_categorized_columns = (
        obs_columns_fixed
        + obsm_coordinate_columns
        + cell_paint_columns
        + cell_paint_nuclear_columns
        + cell_profiler_columns
        + count_columns
        + nuclear_count_columns
        + unassigned_columns
    )
    assert len(column_list) == len(all_categorized_columns), (
        "Column categorization incomplete."
    )

    return (
        obs_columns_fixed,
        obsm_coordinate_columns,
        cell_paint_columns,
        cell_paint_nuclear_columns,
        cell_profiler_columns,
        count_columns,
        nuclear_count_columns,
        unassigned_columns,
    )


def retrieve_input_data(cells2stats_output_bundle):
    # Expected folder structure (showing only relevant files):
    # ├── Cytoprofiling/
    # │   └── Instrument/
    # │       └── RawCellStats.parquet
    # └── Panel.json

    required_files = [
        "Panel.json",
        "Cytoprofiling/Instrument/RawCellStats.parquet"
    ]

    if zipfile.is_zipfile(cells2stats_output_bundle):
        cells2stats_output_bundle = extract_selected_files_from_zip(
            cells2stats_output_bundle,
            members=required_files
        )

    assert os.path.isdir(cells2stats_output_bundle), "Input is expected to be a (compressed) directory."
    input_dir = Path(cells2stats_output_bundle)
    input_data = dict(zip(
        ["target_panel", "count_matrix"],
        [input_dir / file for file in required_files]
    ))

    assert all([file.exists() for file in input_data.values()]), (
        f"Not all required input files are found. Make sure that {par['input']} contains {input_data.values()}."
    )

    return input_data


def main():

    logger.info("Reading input data...")
    input_data = retrieve_input_data(par["input"])
    with open(input_data["target_panel"], "r") as f:
        target_panel = json.load(f)
    df = pd.read_parquet(input_data["count_matrix"], engine="pyarrow")
    df_columns = df.columns.tolist()

    logger.info("Categorizing input data...")
    (
        obs_columns_fixed,
        coordinate_columns,
        cell_paint_columns,
        cell_paint_nuclear_columns,
        cell_profiler_columns,
        count_columns,
        nuclear_count_columns,
        unassigned_columns,
    ) = categorize_columns(df_columns, target_panel)

    df = df.set_index(df["Cell"].astype(str), drop=False)
    df.index_name = None

    # var and obs names
    var_names = [var.split(".")[0] for var in count_columns]
    obs_names = df["Cell"].astype(str).tolist()

    # Count matrix
    logger.info("Creating count matrix...")
    count_df = df[count_columns].copy()
    count_matrix_sparse = sp.csr_matrix(count_df.values)

    # Obs field
    logger.info(f"Creating obs field with columns {obs_columns_fixed}")
    obs_df = df[obs_columns_fixed].copy()

    # Create AnnData object
    logger.info("Creating AnnData object...")
    adata = ad.AnnData(
        X=count_matrix_sparse,
        obs=obs_df,
        var=pd.DataFrame(index=var_names),
    )

    adata.obs_names = obs_names
    adata.var_names = var_names

    # Spatial coordinates
    coordinate_sets = {
        par["obsm_coordinates"]: ["X", "Y"],
        f"{par['obsm_coordinates']}_um": ["Xum", "Yum"],
    }

    for obsm_key, coord_cols in coordinate_sets.items():
        if all(col in coordinate_columns for col in coord_cols):
            coordinates = df[coord_cols].copy()
            adata.obsm[obsm_key] = coordinates.values
            adata.uns[obsm_key] = coord_cols
            logger.info(f"Added {obsm_key} coordinates ({coord_cols}) to obsm")
        else:
            missing_cols = [col for col in coord_cols if col not in coordinate_columns]
            logger.warning(
                f"Skipping {obsm_key}: missing coordinate columns {missing_cols}"
            )

    # Add (optional) .obsm fields
    if par["obsm_cell_paint"]:
        logger.info(f"Adding {par['obsm_cell_paint']} to obsm")
        adata.obsm[par["obsm_cell_paint"]] = df[cell_paint_columns].copy()
        adata.uns[par["obsm_cell_paint"]] = cell_paint_columns
    if par["obsm_cell_paint_nuclear"]:
        logger.info(f"Adding {par['obsm_cell_paint_nuclear']} to obsm")
        adata.obsm[par["obsm_cell_paint_nuclear"]] = df[
            cell_paint_nuclear_columns
        ].copy()
        adata.uns[par["obsm_cell_paint_nuclear"]] = cell_paint_nuclear_columns
    if par["obsm_cell_profiler"]:
        logger.info(f"Adding {par['obsm_cell_profiler']} to obsm")
        adata.obsm[par["obsm_cell_profiler"]] = df[cell_profiler_columns].copy()
        adata.uns[par["obsm_cell_profiler"]] = cell_profiler_columns
    if par["obsm_unassigned_targets"]:
        logger.info(f"Adding {par['obsm_unassigned_targets']} to obsm")
        adata.obsm["unassigned_targets"] = df[unassigned_columns].copy()
        adata.uns["unassigned_targets"] = unassigned_columns

    # Add (optional) nuclear count layer
    if par["layer_nuclear_counts"]:
        assert_matching_order(
            var_names, nuclear_count_columns, split_pattern="_Nuclear"
        )
        logger.info(f"Adding {par['layer_nuclear_counts']} to layers")
        nuclear_count_df = df[nuclear_count_columns].copy()
        nuclear_count_matrix_sparse = sp.csr_matrix(nuclear_count_df.values)
        adata.layers[par["layer_nuclear_counts"]] = nuclear_count_matrix_sparse

    # Write output MuData
    logger.info("Writing MuData object...")
    mdata = mu.MuData({par["modality"]: adata})
    mdata.write_h5mu(par["output"], compression=par["output_compression"])


if __name__ == "__main__":
    main()
