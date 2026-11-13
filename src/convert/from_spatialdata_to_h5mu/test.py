import pytest
import sys
import mudata as mu
import spatialdata as sd


def test_simple_execution(run_component, tmp_path):
    output = tmp_path / "output.h5mu"

    run_component(
        [
            "--input",
            meta["resources_dir"] + "/xenium_tiny.zarr",
            "--output",
            output,
        ]
    )
    assert output.is_file(), "output file was not created"

    mdata = mu.read_h5mu(output)
    assert list(mdata.mod.keys()) == ["rna"], "Expected modality rna"

    adata = mdata.mod["rna"]

    # TODO: update what is checked here when spatialdata from other experimental set-ups are tested (e.g. cosmx, visium)
    assert list(adata.obs.keys()) == [
        "cell_id",
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
        "region",
        "z_level",
        "cell_labels",
    ]

    assert list(adata.uns.keys()) == ["spatialdata_attrs"]
    assert list(adata.obsm.keys()) == ["spatial"]
    assert list(adata.var.keys()) == ["gene_ids", "feature_types", "genome"]

    assert all(adata.var["feature_types"] == "Gene Expression")
    assert adata.obsm["spatial"].dtype == "float"


def test_input_table(run_component, tmp_path):
    sdata = sd.read_zarr(meta["resources_dir"] + "/xenium_tiny.zarr")
    sdata["cells"] = sdata["table"]
    del sdata["table"]
    input = tmp_path / "input.zarr"
    sdata.write(input)
    output = tmp_path / "output.h5mu"

    run_component(["--input", input, "--input_table", "cells", "--output", output])

    mdata = mu.read_h5mu(output)
    assert mdata.mod["rna"].n_obs == sdata["cells"].n_obs


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
