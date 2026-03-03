import pytest
import json
import subprocess
import re

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


@pytest.fixture
def nichenet_io(tmp_path):
    nichenet_distr = tmp_path / "nichenet_distr.svg"
    return {
        "--create_nichenet_gene_program_mask": "True",
        "--create_omnipath_gene_program_mask": "False",
        "--create_mebocost_gene_program_mask": "False",
        "--create_collectri_tf_gene_program_mask": "False",
        "--input_gene_orthologs_mapping_file": ortholog_file,
        "--input_nichenet_lrt_network": nichenet_lr_network_file,
        "--input_nichenet_ligand_target_matrix": nichenet_matrix_file,
        "--output_nichenet_gp_gene_count_distributions": nichenet_distr,
    }


@pytest.fixture
def omnipath_io(tmp_path):
    omnipath_distr = tmp_path / "omnipath_distr.svg"
    return {
        "--create_nichenet_gene_program_mask": "False",
        "--create_omnipath_gene_program_mask": "True",
        "--create_mebocost_gene_program_mask": "False",
        "--create_collectri_tf_gene_program_mask": "False",
        "--input_gene_orthologs_mapping_file": ortholog_file,
        "--input_omnipath_lr_network": omnipath_lr_network_file,
        "--output_omnipath_gp_gene_count_distributions": omnipath_distr,
    }


@pytest.fixture
def mebocost_io(tmp_path):
    mebocost_distr = tmp_path / "mebocost_distr.svg"
    return {
        "--create_nichenet_gene_program_mask": "False",
        "--create_omnipath_gene_program_mask": "False",
        "--create_mebocost_gene_program_mask": "True",
        "--create_collectri_tf_gene_program_mask": "False",
        "--input_metabolite_enzymes": enzymes_file,
        "--input_metabolite_sensors": sensors_file,
        "--output_mebocost_gp_gene_count_distributions": mebocost_distr,
    }


@pytest.fixture
def collectri_io(tmp_path):
    collectri_distr = tmp_path / "collectri_distr.svg"
    return {
        "--create_nichenet_gene_program_mask": "False",
        "--create_omnipath_gene_program_mask": "False",
        "--create_mebocost_gene_program_mask": "False",
        "--create_collectri_tf_gene_program_mask": "True",
        "--input_gene_orthologs_mapping_file": ortholog_file,
        "--input_collectri_tf_network": collectri_tf_network_file,
        "--input_collectri_ligand_target_matrix": collectri_tf_network_file,
        "--output_collectri_tf_gp_gene_count_distributions": collectri_distr,
    }


@pytest.fixture
def gene_program_io(nichenet_io, omnipath_io, mebocost_io, collectri_io):
    return {
        "nichenet": nichenet_io,
        "omnipath": omnipath_io,
        "mebocost": mebocost_io,
        "collectri": collectri_io,
    }


@pytest.mark.parametrize(
    "mask_type,expected_gp_keys",
    [
        ("nichenet", ["ligand_receptor_target_gene_GP"]),
        ("omnipath", ["ligand_receptor_GP", "combined_GP"]),
        ("mebocost", ["metabolite_enzyme_sensor_GP", "combined_GP"]),
        ("collectri", ["TF_target_genes_GP"]),
    ],
)
def test_api_execution(run_component, tmp_path):
    output = tmp_path / "output.json"

    args = [
        "--input_gene_orthologs_mapping_file",
        ortholog_file,
        "--input_metabolite_enzymes",
        enzymes_file,
        "--input_metabolite_sensors",
        sensors_file,
        "--create_omnipath_gene_program_mask",
        "True",
        "--create_nichenet_gene_program_mask",
        "True",
        "--create_mebocost_gene_program_mask",
        "True",
        "--create_collectri_tf_gene_program_mask",
        "True",
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
        "ligand_receptor_GP",  # omnipath
        "ligand_receptor_target_gene_GP",  # nichenet
        "metabolite_enzyme_sensor_GP",  # mebocost
        "TF_target_genes_GP",  # collectri
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


def test_io_gene_program_mask(
    run_component, tmp_path, mask_type, expected_gp_keys, gene_program_io
):
    output = tmp_path / "output.json"
    io_args = gene_program_io[mask_type]

    args = [
        "--species",
        "mouse",
        "--output",
        output,
    ]

    for flag, value in io_args.items():
        args.extend([flag, value])

    run_component(args)

    expected_outputs = [output]
    for flag, value in io_args.items():
        if flag.startswith("--output"):
            expected_outputs.append(value)

    for output_file in expected_outputs:
        assert output_file.is_file(), (
            f"Expected output file {output_file} does not exist"
        )

    # Read gene program mask
    with open(output, "r") as f:
        gp_mask = json.load(f)

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


def test_fail_no_masks(run_component, tmp_path):
    output = tmp_path / "output.json"
    # fails because input data are not correctly lognormalized
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input_gene_orthologs_mapping_file",
                ortholog_file,
                "--input_metabolite_enzymes",
                enzymes_file,
                "--input_metabolite_sensors",
                sensors_file,
                "--create_omnipath_gene_program_mask",
                "False",
                "--create_nichenet_gene_program_mask",
                "False",
                "--create_mebocost_gene_program_mask",
                "False",
                "--create_collectri_tf_gene_program_mask",
                "False",
                "--species",
                "mouse",
                "--output",
                output,
            ]
        )
    assert re.search(
        r"At least one gene program mask must be set to True",
        err.value.stdout.decode("utf-8"),
    )


def test_fail_missing_omnipath_orthologs(run_component, tmp_path):
    output = tmp_path / "output.json"

    # fails because omnipath mask creation requires gene ortholog mapping
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--create_omnipath_gene_program_mask",
                "True",
                "--create_nichenet_gene_program_mask",
                "False",
                "--create_mebocost_gene_program_mask",
                "False",
                "--create_collectri_tf_gene_program_mask",
                "False",
                "--species",
                "mouse",
                "--output",
                output,
            ]
        )
    assert re.search(
        r"Mouse species requires --input_gene_orthologs_mapping_file to generate the Omnipath mask.",
        err.value.stdout.decode("utf-8"),
    )


def test_fail_missing_nichenet_orthologs(run_component, tmp_path):
    output = tmp_path / "output.json"

    # fails because nichenet mask creation requires gene ortholog mapping
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--create_omnipath_gene_program_mask",
                "False",
                "--create_nichenet_gene_program_mask",
                "True",
                "--create_mebocost_gene_program_mask",
                "False",
                "--create_collectri_tf_gene_program_mask",
                "False",
                "--species",
                "mouse",
                "--output",
                output,
            ]
        )
    assert re.search(
        r"Mouse species requires --input_gene_orthologs_mapping_file to generate the NicheNet mask.",
        err.value.stdout.decode("utf-8"),
    )


def test_fail_missing_mebocost_metabolites(run_component, tmp_path):
    output = tmp_path / "output.json"

    # fails because mebocost mask creation requires metabolite files
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--create_omnipath_gene_program_mask",
                "False",
                "--create_nichenet_gene_program_mask",
                "False",
                "--create_mebocost_gene_program_mask",
                "True",
                "--create_collectri_tf_gene_program_mask",
                "False",
                "--species",
                "mouse",
                "--output",
                output,
            ]
        )
    assert re.search(
        r"MeBocost mask requires --input_metabolite_enzymes and --input_metabolite_sensors.",
        err.value.stdout.decode("utf-8"),
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
