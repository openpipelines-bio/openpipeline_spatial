import json
import os
import shutil
import sys
from pathlib import Path

import pandas as pd

from nichecompass.utils import (
    extract_gp_dict_from_mebocost_ms_interactions,
    extract_gp_dict_from_nichenet_lrt_interactions,
    extract_gp_dict_from_omnipath_lr_interactions,
    filter_and_combine_gp_dict_gps_v2,
    extract_gp_dict_from_collectri_tf_network,
)


## VIASH START
par = {
    "species": "mouse",
    "create_omnipath_gene_program_mask": False,
    "create_nichenet_gene_program_mask": False,
    "create_mebocost_gene_program_mask": False,
    "create_collectri_tf_gene_program_mask": True,
    # omnipath params
    "input_gene_orthologs_mapping_file": "resources_test/niche/human_mouse_gene_orthologs.csv",
    "omnipath_min_curation_effort": 2,
    "input_omnipath_lr_network": "resources_test/niche/omnipath_lr_network.csv",
    # "input_omnipath_lr_network": None,
    # nichenet params
    "input_nichenet_lrt_network": "resources_test/niche/nichenet_lr_network.csv",
    "input_nichenet_ligand_target_matrix": "resources_test/niche/nichenet_ligand_target_matrix_v2_mouse.csv",
    "nichenet_version": "v2",
    "nichenet_keep_target_genes_ratio": 1.0,
    "nichenet_max_n_target_genes_per_gp": 250,
    # mebocost_gene_program_mask
    "input_metabolite_enzymes": "resources_test/niche/mouse_metabolite_enzymes.tsv",
    "input_metabolite_sensors": "resources_test/niche/mouse_metabolite_sensors.tsv",
    # collectri params
    "input_collectri_tf_network": "resources_test/niche/collectri_tf_network.csv",
    # filter and combine programs
    "overlap_thresh_target_genes": 1.0,
    # output paths
    "output": "collectri_gp_mask.json",
    "output_omnipath_lr_network": "resources_test/niche/omnipath_lr_network.csv",
    "output_nichenet_lrt_network": "resources_test/niche/nichenet_lr_network.csv",
    "output_nichenet_ligand_target_matrix": "nichenet_ligand_target_matrix_v2_mouse.csv",
    "output_collectri_tf_network": "collectri_tf_network.csv",
    "output_omnipath_gp_gene_count_distributions": "omnipath_gp_gene_count_distributions.svg",
    "output_nichenet_gp_gene_count_distributions": "nichenet_gp_gene_count_distributions.svg",
    "output_mebocost_gp_gene_count_distributions": "mebocost_gp_gene_count_distributions.svg",
    "output_collectri_tf_gp_gene_count_distributions": "collectri_tf_gp_gene_count_distributions.svg",
}

meta = {"temp_dir": "tmp/", "resources_dir": "src/utils/"}
## VIASH END
sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()


def _sanitize_omnipath_csv(file_path: str) -> None:
    """
    Sanitize the OmniPath CSV file format for compatibility with NicheCompass
    `load_from_disk`.

    Bug 1 - CSV index inconsistency:
    The library saves with `to_csv(path, index=False)` but loads with
    `read_csv(path, index_col=0)`. This causes the first data column to be
    incorrectly treated as the index when loading.

    Bug 2 - Missing NaN handling:
    OmniPath contains proteins (TrEMBL/unreviewed entries like A0A2R8YE73) that
    don't have gene symbol mappings. These appear as NaN in genesymbol columns.
    The `resolve_protein_complexes` function doesn't handle NaN, causing
    `TypeError`. When fetching from the API, groupby operations implicitly
    filter some NaN rows, but `load_from_disk` doesn't have this filtering.

    This function sanitizes the loaded CSV by:
     1. Dropping rows with NaN in genesymbol columns
         (unusable for gene programs)
    2. Rewriting with an index column (to fix the load inconsistency)

    See:
    https://github.com/Lotfollahi-lab/nichecompass/blob/main/src/nichecompass/utils/gene_programs.py

    Args:
        file_path: Path to the OmniPath CSV file to fix.
    """
    if not file_path or not os.path.exists(file_path):
        return

    df = pd.read_csv(file_path)

    # Find all genesymbol columns dynamically and drop rows where any are NaN.
    # These proteins don't have gene mappings and cannot be used in gene programs.
    genesymbol_cols = df.columns.str.lower().str.contains("genesymbol")
    if genesymbol_cols.any():
        df = df.dropna(subset=df.columns[genesymbol_cols])

    # Rewrite with an index column so downstream `read_csv(..., index_col=0)` is stable.
    os.makedirs(meta["temp_dir"], exist_ok=True)

    lr_network_file_path = os.path.join(
        meta["temp_dir"],
        "input_omnipath_lr_network_sanitized.csv",
    )
    df.to_csv(lr_network_file_path, index=True)
    return lr_network_file_path


def create_omnipath_gene_program_mask(
    output_lr_network: str | None,
    output_count_distr: str | None,
    input_lr_network: str | None,
    input_orthologs: str | None,
) -> dict:
    # Generate omnipath gene program mask
    # Determine output distribution
    plot_gp_gene_count_distributions = bool(output_count_distr)

    # Determine load_from_disk and save_to_disk from I/O params.
    load_from_disk = bool(input_lr_network)
    save_to_disk = bool(output_lr_network) and (not load_from_disk)

    # Warn if both input and output are provided
    if input_lr_network and output_lr_network:
        logger.warning(
            "Both Omnipath input and output paths are provided. "
            "Using input file; output will not be saved."
        )

    # Use input file path if provided, otherwise use output file path.
    if load_from_disk:
        lr_network_file_path = _sanitize_omnipath_csv(input_lr_network)
    else:
        lr_network_file_path = output_lr_network

    omnipath_gp_dict = extract_gp_dict_from_omnipath_lr_interactions(
        species=par["species"],
        min_curation_effort=par["omnipath_min_curation_effort"],
        load_from_disk=load_from_disk,
        save_to_disk=save_to_disk,
        lr_network_file_path=lr_network_file_path,
        gene_orthologs_mapping_file_path=input_orthologs,
        plot_gp_gene_count_distributions=plot_gp_gene_count_distributions,
        gp_gene_count_distributions_save_path=output_count_distr,
    )

    return omnipath_gp_dict


def create_nichenet_gene_program_mask(
    output_lrt_network: str | None,
    output_lt_matrix: str | None,
    output_count_distr: str | None,
    input_lrt_network: str | None,
    input_lt_matrix: str | None,
    input_orthologs: str | None,
) -> dict:
    plot_gp_gene_count_distributions = bool(output_count_distr)

    # Validate NicheNet I/O.
    load_from_disk = bool(input_lrt_network) and bool(input_lt_matrix)

    save_to_disk = (
        bool(output_lrt_network) or bool(output_lt_matrix) and not load_from_disk
    )

    # Warn if both input and output are provided.
    if load_from_disk and (bool(output_lrt_network) or bool(output_lt_matrix)):
        logger.warning(
            "Both NicheNet input and output paths are provided. "
            "Using input files; outputs will not be saved."
        )

    # Use input file path if provided, otherwise use output file path.
    if load_from_disk:
        lr_network_file_path = input_lrt_network
        ligand_target_matrix_file_path = input_lt_matrix
    else:
        lr_network_file_path = output_lrt_network
        ligand_target_matrix_file_path = output_lt_matrix

    nichenet_gp_dict = extract_gp_dict_from_nichenet_lrt_interactions(
        species=par["species"],
        version=par["nichenet_version"],
        keep_target_genes_ratio=par["nichenet_keep_target_genes_ratio"],
        max_n_target_genes_per_gp=par["nichenet_max_n_target_genes_per_gp"],
        load_from_disk=load_from_disk,
        save_to_disk=save_to_disk,
        lr_network_file_path=lr_network_file_path,
        ligand_target_matrix_file_path=ligand_target_matrix_file_path,
        gene_orthologs_mapping_file_path=input_orthologs,
        plot_gp_gene_count_distributions=plot_gp_gene_count_distributions,
        gp_gene_count_distributions_save_path=output_count_distr,
    )
    return nichenet_gp_dict


def create_mebocost_gene_program_mask(
    output_count_distr: str | None,
    input_metabolite_enzymes: str,
    input_metabolite_sensors: str,
) -> dict:
    os.makedirs(meta["temp_dir"], exist_ok=True)

    metabolite_enzymes_path = os.path.join(
        meta["temp_dir"],
        f"{par['species']}_metabolite_enzymes.tsv",
    )
    metabolite_sensors_path = os.path.join(
        meta["temp_dir"],
        f"{par['species']}_metabolite_sensors.tsv",
    )

    shutil.copy2(
        input_metabolite_enzymes,
        metabolite_enzymes_path,
    )
    shutil.copy2(
        input_metabolite_sensors,
        metabolite_sensors_path,
    )
    plot_gp_gene_count_distributions = bool(output_count_distr)

    mebocost_gp_dict = extract_gp_dict_from_mebocost_ms_interactions(
        dir_path=meta["temp_dir"],
        species=par["species"],
        plot_gp_gene_count_distributions=plot_gp_gene_count_distributions,
        gp_gene_count_distributions_save_path=output_count_distr,
    )
    return mebocost_gp_dict


def create_collectri_tf_gene_program_mask(
    output_count_distr: str | None,
    output_tf_network: str | None,
    input_tf_network: str | None,
) -> dict:
    plot_gp_gene_count_distributions = bool(output_count_distr)

    # Determine load_from_disk and save_to_disk from I/O params.
    load_from_disk = bool(input_tf_network)
    save_to_disk = bool(output_tf_network) and not load_from_disk

    # Warn if both input and output are provided
    if input_tf_network and output_tf_network:
        logger.warning(
            "Both CollecTRI input and output paths are provided. "
            "Using input file; output will not be saved."
        )

    # Use input file path if provided, otherwise use output file path
    tf_network_file_path = input_tf_network if load_from_disk else output_tf_network

    collectri_gp_dict = extract_gp_dict_from_collectri_tf_network(
        species=par["species"],
        load_from_disk=load_from_disk,
        save_to_disk=save_to_disk,
        tf_network_file_path=tf_network_file_path,
        plot_gp_gene_count_distributions=plot_gp_gene_count_distributions,
        gp_gene_count_distributions_save_path=output_count_distr,
    )
    return collectri_gp_dict


def main():
    # Validate that inputs are provided correctly
    if not any(
        [
            par["create_omnipath_gene_program_mask"],
            par["create_nichenet_gene_program_mask"],
            par["create_mebocost_gene_program_mask"],
            par["create_collectri_tf_gene_program_mask"],
        ]
    ):
        raise ValueError("At least one gene program mask must be set to True")
    if (
        par["create_omnipath_gene_program_mask"]
        and par["species"] == "mouse"
        and not par["input_gene_orthologs_mapping_file"]
    ):
        raise ValueError(
            "Mouse species requires --input_gene_orthologs_mapping_file "
            "to generate the Omnipath mask."
        )
    if (
        par["create_nichenet_gene_program_mask"]
        and par["species"] == "mouse"
        and not par["input_gene_orthologs_mapping_file"]
    ):
        raise ValueError(
            "Mouse species requires --input_gene_orthologs_mapping_file "
            "to generate the NicheNet mask."
        )
    if par["create_mebocost_gene_program_mask"] and (
        (not par["input_metabolite_enzymes"]) or (not par["input_metabolite_sensors"])
    ):
        raise ValueError(
            "MeBocost mask requires --input_metabolite_enzymes "
            "and --input_metabolite_sensors."
        )

    # Assemble gene program dictionaries
    gp_dicts = []

    masks = {
        "create_omnipath_gene_program_mask": (
            "Omnipath",
            create_omnipath_gene_program_mask,
        ),
        "create_nichenet_gene_program_mask": (
            "NicheNet",
            create_nichenet_gene_program_mask,
        ),
        "create_mebocost_gene_program_mask": (
            "MeBocost",
            create_mebocost_gene_program_mask,
        ),
        "create_collectri_tf_gene_program_mask": (
            "CollecTRI TF",
            create_collectri_tf_gene_program_mask,
        ),
    }

    mask_args = {
        "create_omnipath_gene_program_mask": (
            par["output_omnipath_lr_network"],
            par["output_omnipath_gp_gene_count_distributions"],
            par["input_omnipath_lr_network"],
            par["input_gene_orthologs_mapping_file"],
        ),
        "create_nichenet_gene_program_mask": (
            par["output_nichenet_lrt_network"],
            par["output_nichenet_ligand_target_matrix"],
            par["output_nichenet_gp_gene_count_distributions"],
            par["input_nichenet_lrt_network"],
            par["input_nichenet_ligand_target_matrix"],
            par["input_gene_orthologs_mapping_file"],
        ),
        "create_mebocost_gene_program_mask": (
            par["output_mebocost_gp_gene_count_distributions"],
            par["input_metabolite_enzymes"],
            par["input_metabolite_sensors"],
        ),
        "create_collectri_tf_gene_program_mask": (
            par["output_collectri_tf_gp_gene_count_distributions"],
            par["output_collectri_tf_network"],
            par["input_collectri_tf_network"],
        ),
    }

    for mask, (mask_name, mask_function) in masks.items():
        if par[mask]:
            logger.info(f"Generating {mask_name} gene program mask...")
            gp_dict = mask_function(*mask_args[mask])
            gp_dicts.append(gp_dict)

    # Filter and combine GPs
    assert len(gp_dicts) > 0, "No gene program dictionaries were created."

    combined_gp_dict = filter_and_combine_gp_dict_gps_v2(
        gp_dicts,
        overlap_thresh_target_genes=par["overlap_thresh_target_genes"],
        verbose=True,
    )

    logger.info("Gene program mask generation completed.")
    logger.info(
        "Number of gene programs after filtering and combining: %s.",
        len(combined_gp_dict),
    )

    output_path = Path(par["output"])
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info("Saving combined gene program mask to: %s", str(output_path))
    with output_path.open("w", encoding="utf-8") as f:
        json.dump(combined_gp_dict, f)


if __name__ == "__main__":
    main()
