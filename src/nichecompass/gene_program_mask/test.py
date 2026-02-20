import pytest
from pathlib import Path
import json

## VIASH START
meta = {
    "executable": "./target/executable/nichecompass/gene_program_mask/gene_program_mask",
    "resources_dir": "./resources_test/niche/",
}
## VIASH END

import sys

sys.path.append(meta["resources_dir"])

ortholog_file = f"{meta['resources_dir']}/niche/human_mouse_gene_orthologs.csv"
enzymes_file = f"{meta['resources_dir']}/niche/mouse_metabolite_enzymes.tsv"
sensors_file = f"{meta['resources_dir']}/niche/mouse_metabolite_sensors.tsv"
omnipath_lr_network_file = f"{meta['resources_dir']}/niche/omnipath_lr_network.csv"
nichenet_lr_network_file = f"{meta['resources_dir']}/niche/nichenet_lr_network.csv"
nichenet_matrix_file = (
    f"{meta['resources_dir']}/niche/nichenet_ligand_target_matrix_v2_mouse.csv"
)
collectri_tf_network_file = f"{meta['resources_dir']}/niche/collectri_tf_network.csv"


def test_api_execution(run_component, tmp_path):
    output = tmp_path / "output.json"

    args = [
        "--input_gene_orthologs_mapping_file",
        ortholog_file,
        "--input_metabolite_enzymes",
        enzymes_file,
        "--input_metabolite_sensors",
        sensors_file,
        "--species",
        "mouse",
        "--output",
        output,
    ]

    run_component(args)

    # check files
    assert output.is_file(), "Output file does not exist"

    # Read gene program mask
    with open(output, "r") as f:
        gp_mask = json.load(f)

    expected_gp_keys = [
        "ligand_receptor_GP",
        "ligand_receptor_target_gene_GP",
        "metabolite_enzyme_sensor_GP",
        "TF_target_genes_GP",
        "combined_GP",
    ]
    matching_gp = []
    for key in expected_gp_keys:
        assert any(key in gp for gp in gp_mask.keys()), (
            f"No gene programs containing '{key}' found"
        )

        gp = next(gp for gp in gp_mask.keys() if key in gp)
        matching_gp.append(gp)

    for gp in matching_gp:
        expected_keys = [
            "sources",
            "targets",
            "sources_categories",
            "targets_categories",
        ]
        assert all([key in gp_mask[gp] for key in expected_keys]), (
            f"Gene program {gp} is missing expected keys"
        )


def test_omnipath_gene_program_mask(run_component, tmp_path):
    output = tmp_path / "output.json"
    omnipath_distr = tmp_path / "omnipath_distr.svg"

    args = [
        "--create_omnipath_gene_program_mask",
        "True",
        "--create_nichenet_gene_program_mask",
        "False",
        "--create_mebocost_gene_program_mask",
        "False",
        "--create_collectri_tf_gene_program_mask",
        "False",
        "--input_gene_orthologs_mapping_file",
        ortholog_file,
        "--input_omnipath_lr_network",
        omnipath_lr_network_file,
        "--species",
        "mouse",
        "--output",
        output,
        "--output_omnipath_gp_gene_count_distributions",
        omnipath_distr,
    ]

    run_component(args)

    expected_outputs = [output, omnipath_distr]

    for output_file in expected_outputs:
        assert output_file.is_file(), (
            f"Expected output file {output_file} does not exist"
        )


def test_nichenet_gene_program_mask(run_component, tmp_path):
    output = tmp_path / "output.json"
    nichenet_distr = tmp_path / "nichenet_distr.svg"

    args = [
        "--create_omnipath_gene_program_mask",
        "False",
        "--create_nichenet_gene_program_mask",
        "True",
        "--create_mebocost_gene_program_mask",
        "False",
        "--create_collectri_tf_gene_program_mask",
        "False",
        "--input_gene_orthologs_mapping_file",
        ortholog_file,
        "--input_nichenet_lrt_network",
        nichenet_lr_network_file,
        "--input_nichenet_ligand_target_matrix",
        nichenet_matrix_file,
        "--species",
        "mouse",
        "--output",
        output,
        "--output_nichenet_gp_gene_count_distributions",
        nichenet_distr,
    ]
    assert Path(nichenet_lr_network_file).is_file(), (
        f"NicheNet LR network file {nichenet_lr_network_file} does not exist"
    )

    run_component(args)

    expected_outputs = [
        output,
        nichenet_distr,
        # nichenet_lt_matrix
    ]

    for output_file in expected_outputs:
        assert output_file.is_file(), (
            f"Expected output file {output_file} does not exist"
        )


def test_mebocost_gene_program_mask(run_component, tmp_path):
    output = tmp_path / "output.json"
    mebocost_distr = tmp_path / "omnipath_distr.svg"

    args = [
        "--create_omnipath_gene_program_mask",
        "False",
        "--create_nichenet_gene_program_mask",
        "False",
        "--create_mebocost_gene_program_mask",
        "True",
        "--create_collectri_tf_gene_program_mask",
        "False",
        "--input_metabolite_enzymes",
        enzymes_file,
        "--input_metabolite_sensors",
        sensors_file,
        "--species",
        "mouse",
        "--output",
        output,
        "--output_mebocost_gp_gene_count_distributions",
        mebocost_distr,
    ]

    run_component(args)

    expected_outputs = [output, mebocost_distr]

    for output_file in expected_outputs:
        assert output_file.is_file(), (
            f"Expected output file {output_file} does not exist"
        )


def test_collectri_tf_gene_program_mask(run_component, tmp_path):
    output = tmp_path / "output.json"
    collectri_distr = tmp_path / "collectri_distr.svg"

    args = [
        "--create_omnipath_gene_program_mask",
        "False",
        "--create_nichenet_gene_program_mask",
        "False",
        "--create_mebocost_gene_program_mask",
        "False",
        "--create_collectri_tf_gene_program_mask",
        "True",
        "--input_gene_orthologs_mapping_file",
        ortholog_file,
        "--input_collectri_tf_network",
        collectri_tf_network_file,
        "--input_collectri_ligand_target_matrix",
        collectri_tf_network_file,
        "--species",
        "mouse",
        "--output",
        output,
        "--output_collectri_tf_gp_gene_count_distributions",
        collectri_distr,
    ]

    run_component(args)

    expected_outputs = [output, collectri_distr]

    for output_file in expected_outputs:
        assert output_file.is_file(), (
            f"Expected output file {output_file} does not exist"
        )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
