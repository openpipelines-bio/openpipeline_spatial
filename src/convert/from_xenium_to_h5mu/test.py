import pytest
import sys
import subprocess
import mudata as mu

## VIASH START
meta = {
    "executable": "./target/executable/convert/from_xenium_to_h5mu/from_xenium_to_h5mu",
    "resources_dir": "resources_test/",
    "config": "src/convert/from_xenium_to_h5mu/config.vsh.yaml",
}
## VIASH END

input = f"{meta['resources_dir']}/xenium_tiny"


def test_simple_execution(run_component, tmp_path):
    output = tmp_path / "xenium.h5mu"

    # run component
    run_component(
        ["--input", input, "--output", str(output), "--output_compression", "gzip"]
    )

    assert output.is_file(), "output file was not created"

    mdata = mu.read_h5mu(output)
    assert list(mdata.mod.keys()) == ["rna"], "Expected modality rna"
    adata = mdata.mod["rna"]

    assert list(adata.obs.keys()) == [
        "transcript_counts",
        "control_probe_counts",
        "genomic_control_counts",
        "control_codeword_counts",
        "unassigned_codeword_counts",
        "deprecated_codeword_counts",
        "total_counts",
        "cell_area",
        "nucleus_area",
        "nucleus_count",
        "segmentation_method",
    ]

    assert list(adata.uns.keys()) == ["xenium_experiment", "xenium_metrics"]
    assert list(adata.obsm.keys()) == ["spatial"]
    assert list(adata.var.keys()) == ["gene_ids", "feature_types", "genome"]

    assert adata.X.dtype.kind == "f"
    assert all(adata.var["feature_types"] == "Gene Expression")
    assert adata.obsm["spatial"].dtype == "float"
    obs_counts = [
        "transcript_counts",
        "control_probe_counts",
        "genomic_control_counts",
        "unassigned_codeword_counts",
        "deprecated_codeword_counts",
        "total_counts",
        "nucleus_count",
    ]
    assert all([adata.obs[obs].dtype == "int" for obs in obs_counts])
    obs_areas = ["cell_area", "nucleus_area"]
    assert all([adata.obs[obs].dtype == "float" for obs in obs_areas])


def test_compressed_input(run_component, tmp_path):
    output = tmp_path / "xenium.h5mu"
    zipped_input = tmp_path / "xenium_tiny.zip"

    subprocess.run(
        ["zip", "-r", str(zipped_input), "."],
        cwd=input,
        check=True
    )

    # run component
    run_component(
        ["--input", zipped_input, "--output", str(output), "--output_compression", "gzip"]
    )

    assert output.is_file(), "output file was not created"

    mdata = mu.read_h5mu(output)
    assert list(mdata.mod.keys()) == ["rna"], "Expected modality rna"
    adata = mdata.mod["rna"]

    assert list(adata.obs.keys()) == [
        "transcript_counts",
        "control_probe_counts",
        "genomic_control_counts",
        "control_codeword_counts",
        "unassigned_codeword_counts",
        "deprecated_codeword_counts",
        "total_counts",
        "cell_area",
        "nucleus_area",
        "nucleus_count",
        "segmentation_method",
    ]

    assert list(adata.uns.keys()) == ["xenium_experiment", "xenium_metrics"]
    assert list(adata.obsm.keys()) == ["spatial"]
    assert list(adata.var.keys()) == ["gene_ids", "feature_types", "genome"]

    assert adata.X.dtype.kind == "f"
    assert all(adata.var["feature_types"] == "Gene Expression")
    assert adata.obsm["spatial"].dtype == "float"
    obs_counts = [
        "transcript_counts",
        "control_probe_counts",
        "genomic_control_counts",
        "unassigned_codeword_counts",
        "deprecated_codeword_counts",
        "total_counts",
        "nucleus_count",
    ]
    assert all([adata.obs[obs].dtype == "int" for obs in obs_counts])
    obs_areas = ["cell_area", "nucleus_area"]
    assert all([adata.obs[obs].dtype == "float" for obs in obs_areas])


def test_rename_fields(run_component, tmp_path):
    output = tmp_path / "xenium.h5mu"

    # run component
    run_component(
        [
            "--input",
            input,
            "--output",
            str(output),
            "--obsm_coordinates",
            "test_coord",
            "--uns_experiment",
            "test_experiment",
            "--uns_metrics",
            "test_metrics",
            "--output_compression",
            "gzip",
        ]
    )

    assert output.is_file(), "output file was not created"

    mdata = mu.read_h5mu(output)
    assert list(mdata.mod.keys()) == ["rna"]
    adata = mdata.mod["rna"]

    assert list(adata.uns.keys()) == ["test_experiment", "test_metrics"]
    assert list(adata.obsm.keys()) == ["test_coord"]


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
