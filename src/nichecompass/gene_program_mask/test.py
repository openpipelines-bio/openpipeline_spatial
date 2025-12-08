import pytest

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


def test_simple_execution(run_component, tmp_path):
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
        "igand_receptor_GP",
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


def test_outputs(run_component, tmp_path):
    output = tmp_path / "output.json"
    omnipath_lr = tmp_path / "omnipath_lr_network.tsv"
    nichenet_lr = tmp_path / "nichenet_lr_network.tsv"
    nichenet_lt = tmp_path / "nichenet_ligand_target_matrix.csv"
    collectri_tf = tmp_path / "output_collectri_tf_network.csv"
    omnipath_distr = tmp_path / "omnipath_distr.svg"
    nichenet_distr = tmp_path / "nichenet_distr.svg"
    mebocost_distr = tmp_path / "mebocost_distr.svg"
    collectri_distr = tmp_path / "collectri_distr.svg"

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
        "--output_omnipath_lr_network",
        omnipath_lr,
        "--output_nichenet_lr_network",
        nichenet_lr,
        "--output_nichenet_ligand_target_matrix",
        nichenet_lt,
        "--output_collectri_tf_network",
        collectri_tf,
        "--output_omnipath_gp_gene_count_distributions",
        omnipath_distr,
        "--output_nichenet_gp_gene_count_distributions",
        nichenet_distr,
        "--output_mebocost_gp_gene_count_distributions",
        mebocost_distr,
        "--output_collectri_tf_gp_gene_count_distributions",
        collectri_distr,
    ]

    run_component(args)

    expected_outputs = [
        omnipath_lr,
        nichenet_lr,
        nichenet_lt,
        collectri_tf,
        omnipath_distr,
        nichenet_distr,
        mebocost_distr,
        collectri_distr,
    ]

    for output in expected_outputs:
        assert output.is_file(), f"Expected output file {output} does not exist"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
