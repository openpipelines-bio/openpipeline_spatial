import os
import sys
import shutil
import json

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
    "create_omnipath_gene_program_mask": True,
    "create_nichenet_gene_program_mask": True,
    "create_mebocost_gene_program_mask": True,
    "create_collectri_tf_gene_program_mask": False,
    # omnipath params
    "input_gene_orthologs_mapping_file": "resources_test/niche/human_mouse_gene_orthologs.csv",
    "omnipath_min_curation_effort": 2,
    # nichenet params
    "nichenet_version": "v2",
    "nichenet_keep_target_genes_ratio": 1.0,
    "nichenet_max_n_target_genes_per_gp": 250,
    # mebocost_gene_program_mask
    "input_metabolite_enzymes": "resources_test/niche/mouse_metabolite_enzymes.tsv",
    "input_metabolite_sensors": "resources_test/niche/mouse_metabolite_sensors.tsv",
    # filter and combine programs
    "overlap_thresh_target_genes": 1.0,
    # output paths
    "output": "prior_knowledge_gene_program_mask.json",
    "output_omnipath_lr_network": "omnipath_lr_network.csv",
    "output_nichenet_lr_network": "nichenet_lr_network.csv",
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

    plot_gp_gene_count_distributions = (
        True if par["output_omnipath_gp_gene_count_distributions"] else False
    )
    save_to_disk = True if par["output_omnipath_lr_network"] else False

    omnipath_gp_dict = extract_gp_dict_from_omnipath_lr_interactions(
        species=par["species"],
        min_curation_effort=par["omnipath_min_curation_effort"],
        load_from_disk=False,
        save_to_disk=True,
        lr_network_file_path=par["output_omnipath_lr_network"],
        gene_orthologs_mapping_file_path=par["input_gene_orthologs_mapping_file"],
        plot_gp_gene_count_distributions=plot_gp_gene_count_distributions,
        gp_gene_count_distributions_save_path=par[
            "output_omnipath_gp_gene_count_distributions"
        ],
    )

    gp_dicts.append(omnipath_gp_dict)

if par["create_nichenet_gene_program_mask"]:
    logger.info("Generating NicheNet gene program mask...")

    plot_gp_gene_count_distributions = (
        True if par["output_nichenet_gp_gene_count_distributions"] else False
    )
    save_to_disk = (
        True
        if (
            par["output_nichenet_lr_network"]
            or par["output_nichenet_ligand_target_matrix"]
        )
        else False
    )

    nichenet_gp_dict = extract_gp_dict_from_nichenet_lrt_interactions(
        species=par["species"],
        version=par["nichenet_version"],
        keep_target_genes_ratio=par["nichenet_keep_target_genes_ratio"],
        max_n_target_genes_per_gp=par["nichenet_max_n_target_genes_per_gp"],
        load_from_disk=False,
        save_to_disk=save_to_disk,
        lr_network_file_path=par["output_nichenet_lr_network"],
        ligand_target_matrix_file_path=par["output_nichenet_ligand_target_matrix"],
        gene_orthologs_mapping_file_path=par["input_gene_orthologs_mapping_file"],
        plot_gp_gene_count_distributions=plot_gp_gene_count_distributions,
        gp_gene_count_distributions_save_path=par[
            "output_nichenet_gp_gene_count_distributions"
        ],
    )

    gp_dicts.append(nichenet_gp_dict)

if par["create_mebocost_gene_program_mask"]:
    logger.info("Generating MeBocost gene program mask...")

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

    gp_dicts.append(mebocost_gp_dict)

if par["create_collectri_tf_gene_program_mask"]:
    logger.info("Generating CollecTRI TF gene program mask...")

    plot_gp_gene_count_distributions = (
        True if par["output_collectri_tf_gp_gene_count_distributions"] else False
    )
    save_to_disk = True if par["output_collectri_tf_network"] else False

    collectri_gp_dict = extract_gp_dict_from_collectri_tf_network(
        species=par["species"],
        save_to_disk=save_to_disk,
        tf_network_file_path=par["output_collectri_tf_network"],
        plot_gp_gene_count_distributions=plot_gp_gene_count_distributions,
        gp_gene_count_distributions_save_path=par[
            "output_collectri_tf_gp_gene_count_distributions"
        ],
    )

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
