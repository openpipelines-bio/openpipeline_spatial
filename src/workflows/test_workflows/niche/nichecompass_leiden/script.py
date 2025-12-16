from mudata import read_h5mu
import sys
import pytest

##VIASH START
par = {"input": "nichecompass_leiden/output.h5mu"}

meta = {"resources_dir": "resources_test"}
##VIASH END


def test_run():
    input_mudata = read_h5mu(par["input"])

    expected_mod = ["rna"]
    expected_obsm = ["X_leiden_nichecompass_umap", "nichecompass_latent"]
    expected_obs = ["sample_id", "nichecompass_leiden_1.0"]
    expected_obsp = [
        "spatial_distances",
        "spatial_connectivities",
        "nichecompass_connectivities",
        "nichecompass_distances",
    ]
    expected_varm = [
        "nichecompass_gp_sources",
        "nichecompass_gp_targets",
        "nichecompass_gp_sources_categories",
        "nichecompass_gp_targets_categories",
    ]
    expected_uns = [
        "nichecompass_sources_categories_label_encoder",
        "nichecompass_targets_categories_label_encoder",
        "nichecompass_source_genes_idx",
        "nichecompass_target_genes_idx",
        "nichecompass_genes_idx",
        "nichecompass_gp_names",
        "nichecompass_active_gp_names",
        "nichecompass_neighbors",
        "spatial",
        "xenium_spatial_neighbors",
    ]

    assert all(key in list(input_mudata.mod) for key in expected_mod), (
        f"Input modalities should be: {expected_mod}, found: {input_mudata.mod.keys()}."
    )
    assert all(key in list(input_mudata.mod["rna"].obsm) for key in expected_obsm), (
        f"Input mod['rna'] obsm columns should be: {expected_obsm}, found: {input_mudata.mod['rna'].obsm.keys()}."
    )
    assert all(key in list(input_mudata.mod["rna"].obs) for key in expected_obs), (
        f"Input mod['rna'] obs columns should be: {expected_obs}, found: {input_mudata.mod['rna'].obs.keys()}."
    )
    assert all(key in list(input_mudata.mod["rna"].obsp) for key in expected_obsp), (
        f"Input mod['rna'] obsp columns should be: {expected_obsp}, found: {input_mudata.mod['rna'].obsp.keys()}."
    )
    assert all(key in list(input_mudata.mod["rna"].uns) for key in expected_uns), (
        f"Input mod['rna'] uns columns should be: {expected_uns}, found: {input_mudata.mod['rna'].uns.keys()}."
    )
    assert all(key in list(input_mudata.mod["rna"].varm) for key in expected_varm), (
        f"Input mod['rna'] varm columns should be: {expected_varm}, found: {input_mudata.mod['rna'].varm.keys()}."
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "--import-mode=importlib"]))
