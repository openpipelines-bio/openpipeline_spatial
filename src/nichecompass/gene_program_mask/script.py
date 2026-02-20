import os
import sys
import shutil
import json

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
    "output": "prior_knowledge_gene_program_mask.json",
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


def sanitize_omnipath_csv(file_path: str) -> None:
    """
    Fix the OmniPath CSV file format for compatibility with NicheCompass load_from_disk.

    WORKAROUND for NicheCompass library bugs:

    Bug 1 - CSV index inconsistency:
    The library saves with `to_csv(path, index=False)` but loads with `read_csv(path, index_col=0)`.
    This causes the first data column to be incorrectly treated as the index when loading.

    Bug 2 - Missing NaN handling:
    OmniPath contains proteins (TrEMBL/unreviewed entries like A0A2R8YE73) that don't have
    gene symbol mappings. These appear as NaN in genesymbol columns.
    The resolve_protein_complexes function doesn't handle NaN, causing TypeError.
    When fetching from API, groupby operations implicitly filter some NaN rows, but
    load_from_disk doesn't have this filtering.

    This function fixes the saved CSV by:
    1. Dropping rows with NaN in genesymbol columns (these are unusable for gene programs)
    2. Rewriting with an index column (to fix the load inconsistency)

    See: https://github.com/Lotfollahi-lab/nichecompass/blob/main/src/nichecompass/utils/gene_programs.py

    Args:
        file_path: Path to the OmniPath CSV file to fix.
    """
    if file_path and os.path.exists(file_path):
        df = pd.read_csv(file_path)
        # Find all genesymbol columns dynamically and drop rows where any are NaN
        # These proteins don't have gene mappings and cannot be used in gene programs
        genesymbol_cols = [c for c in df.columns if "genesymbol" in c.lower()]
        if genesymbol_cols:
            df = df.dropna(subset=genesymbol_cols)
        df.to_csv(file_path, index=True)  # Rewrite with index column


def create_omnipath_gene_program_mask():
    # Generate omnipath gene program mask
    # Determine output distribution
    plot_gp_gene_count_distributions = (
        True if par["output_omnipath_gp_gene_count_distributions"] else False
    )

    # Determine load_from_disk and save_to_disk based on input/output parameters
    load_from_disk = bool(par["input_omnipath_lr_network"])
    save_to_disk = bool(par["output_omnipath_lr_network"]) and not load_from_disk

    # Warn if both input and output are provided
    if par["input_omnipath_lr_network"] and par["output_omnipath_lr_network"]:
        logger.warning(
            "Both --input_omnipath_lr_network and --output_omnipath_lr_network are provided. "
            "Using input file (load_from_disk=True), output will not be saved."
        )

    # Use input file path if provided, otherwise use output file path
    lr_network_file_path = (
        par["input_omnipath_lr_network"] or par["output_omnipath_lr_network"]
    )

    omnipath_gp_dict = extract_gp_dict_from_omnipath_lr_interactions(
        species=par["species"],
        min_curation_effort=par["omnipath_min_curation_effort"],
        load_from_disk=load_from_disk,
        save_to_disk=save_to_disk,
        lr_network_file_path=lr_network_file_path,
        gene_orthologs_mapping_file_path=par["input_gene_orthologs_mapping_file"],
        plot_gp_gene_count_distributions=plot_gp_gene_count_distributions,
        gp_gene_count_distributions_save_path=par[
            "output_omnipath_gp_gene_count_distributions"
        ],
    )

    return omnipath_gp_dict


def create_nichenet_gene_program_mask():
    plot_gp_gene_count_distributions = (
        True if par["output_nichenet_gp_gene_count_distributions"] else False
    )
    load_from_disk = bool(par["input_nichenet_lrt_network"])
    save_to_disk = bool(par["output_nichenet_lrt_network"]) and not load_from_disk

    # Warn if both input and output are provided
    if par["input_nichenet_lrt_network"] and par["output_nichenet_lrt_network"]:
        logger.warning(
            "Both --input_nichenet_lrt_network and --output_nichenet_lrt_network are provided. "
            "Using input file (load_from_disk=True), output will not be saved."
        )

    # Use input file path if provided, otherwise use output file path
    lr_network_file_path = (
        par["input_nichenet_lrt_network"] or par["output_nichenet_lrt_network"]
    )
    lf_target_matrix_file_path = (
        par["input_nichenet_ligand_target_matrix"]
        or par["output_nichenet_ligand_target_matrix"]
    )

    nichenet_gp_dict = extract_gp_dict_from_nichenet_lrt_interactions(
        species=par["species"],
        version=par["nichenet_version"],
        keep_target_genes_ratio=par["nichenet_keep_target_genes_ratio"],
        max_n_target_genes_per_gp=par["nichenet_max_n_target_genes_per_gp"],
        load_from_disk=load_from_disk,
        save_to_disk=save_to_disk,
        lr_network_file_path=lr_network_file_path,
        ligand_target_matrix_file_path=lf_target_matrix_file_path,
        gene_orthologs_mapping_file_path=par["input_gene_orthologs_mapping_file"],
        plot_gp_gene_count_distributions=plot_gp_gene_count_distributions,
        gp_gene_count_distributions_save_path=par[
            "output_nichenet_gp_gene_count_distributions"
        ],
    )
    return nichenet_gp_dict


def create_mobocost_gene_program_mask():
    shutil.copy2(
        par["input_metabolite_enzymes"],
        os.path.join(meta["temp_dir"], f"{par['species']}_metabolite_enzymes.tsv"),
    )
    shutil.copy2(
        par["input_metabolite_sensors"],
        os.path.join(meta["temp_dir"], f"{par['species']}_metabolite_sensors.tsv"),
    )
    plot_gp_gene_count_distributions = (
        True if par["output_mebocost_gp_gene_count_distributions"] else False
    )

    mebocost_gp_dict = extract_gp_dict_from_mebocost_ms_interactions(
        dir_path=meta["temp_dir"],
        species=par["species"],
        plot_gp_gene_count_distributions=plot_gp_gene_count_distributions,
        gp_gene_count_distributions_save_path=par[
            "output_mebocost_gp_gene_count_distributions"
        ],
    )
    return mebocost_gp_dict


def create_collectri_tf_gene_program_mask():
    plot_gp_gene_count_distributions = (
        True if par["output_collectri_tf_gp_gene_count_distributions"] else False
    )

    # Determine load_from_disk and save_to_disk based on input/output parameters
    load_from_disk = bool(par["input_collectri_tf_network"])
    save_to_disk = bool(par["output_collectri_tf_network"]) and not load_from_disk

    # Warn if both input and output are provided
    if par["input_collectri_tf_network"] and par["output_collectri_tf_network"]:
        logger.warning(
            "Both --input_collectri_tf_network and --output_collectri_tf_network are provided. "
            "Using input file (load_from_disk=True), output will not be saved."
        )

    # Use input file path if provided, otherwise use output file path
    tf_network_file_path = (
        par["input_collectri_tf_network"] or par["output_collectri_tf_network"]
    )

    collectri_gp_dict = extract_gp_dict_from_collectri_tf_network(
        species=par["species"],
        load_from_disk=load_from_disk,
        save_to_disk=save_to_disk,
        tf_network_file_path=tf_network_file_path,
        plot_gp_gene_count_distributions=plot_gp_gene_count_distributions,
        gp_gene_count_distributions_save_path=par[
            "output_collectri_tf_gp_gene_count_distributions"
        ],
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
            "For mouse species, a --input_gene_orthologs_mapping_file file must be provided for generating the omnipath mask."
        )
    if (
        par["create_nichenet_gene_program_mask"]
        and par["species"] == "mouse"
        and not par["input_gene_orthologs_mapping_file"]
    ):
        raise ValueError(
            "For mouse species, a --input_gene_orthologs_mapping_file file must be provided for generating the nichenet mask."
        )
    if par["create_mebocost_gene_program_mask"] and (
        not par["input_metabolite_enzymes"] or not par["input_metabolite_sensors"]
    ):
        raise ValueError(
            "For mebocost gene program mask, both --input_metabolite_enzymes and --input_metabolite_sensors files must be provided."
        )

    # Assemble gene program dictionaries
    gp_dicts = []

    if par["create_omnipath_gene_program_mask"]:
        logger.info("Generating Omnipath gene program mask...")
        omnipath_gp_dict = create_omnipath_gene_program_mask()
        gp_dicts.append(omnipath_gp_dict)

    if par["create_nichenet_gene_program_mask"]:
        logger.info("Generating NicheNet gene program mask...")
        nichenet_gp_dict = create_nichenet_gene_program_mask()
        gp_dicts.append(nichenet_gp_dict)

    if par["create_mebocost_gene_program_mask"]:
        logger.info("Generating MeBocost gene program mask...")
        mebocost_gp_dict = create_mobocost_gene_program_mask()
        gp_dicts.append(mebocost_gp_dict)

    if par["create_collectri_tf_gene_program_mask"]:
        logger.info("Generating CollecTRI TF gene program mask...")
        collectri_gp_dict = create_collectri_tf_gene_program_mask()
        gp_dicts.append(collectri_gp_dict)

    # Filter and combine GPs
    assert len(gp_dicts) > 0, "No gene program dictionaries were created."

    combined_gp_dict = filter_and_combine_gp_dict_gps_v2(
        gp_dicts,
        overlap_thresh_target_genes=par["overlap_thresh_target_genes"],
        verbose=True,
    )

    logger.info("Gene program mask generation completed.")
    logger.info(
        f"Number of gene programs after filtering and combining: {len(combined_gp_dict)}."
    )

    logger.info(f"Saving combined gene program mask to: {par['output']}")
    with open(par["output"], "w") as f:
        json.dump(combined_gp_dict, f)


if __name__ == "__main__":
    main()
