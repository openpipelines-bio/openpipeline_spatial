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
omnipath_lr_network_file = f"{meta['resources_dir']}/niche/omnipath_lr_network.csv"
collectri_tf_network_file = f"{meta['resources_dir']}/niche/collectri_tf_network.csv"


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


def test_inputs_and_outputs(run_component, tmp_path):
    """Test loading from input files instead of querying APIs.
    
    This test uses pre-downloaded omnipath and collectri network files
    to avoid API rate limits when running tests in parallel.
    """
    output = tmp_path / "output.json"
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
        "--input_omnipath_lr_network",
        omnipath_lr_network_file,
        "--input_collectri_tf_network",
        collectri_tf_network_file,
        "--species",
        "mouse",
        "--output",
        output,
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
        output,
        omnipath_distr,
        nichenet_distr,
        mebocost_distr,
        collectri_distr,
    ]

    for output_file in expected_outputs:
        assert output_file.is_file(), f"Expected output file {output_file} does not exist"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
