import pytest
import sys
import mudata as mu

## VIASH START
meta = {
    "executable": "./target/executable/convert/from_spaceranger_to_h5mu/from_spaceranger_to_h5mu",
    "resources_dir": "resources_test/",
    "config": "src/convert/from_spaceranger_to_h5mu/config.vsh.yaml",
}
## VIASH END

input = f"{meta['resources_dir']}/Visium_FFPE_Human_Ovarian_Cancer_tiny_spaceranger"


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

    assert list(adata.uns.keys()) == [
        "metrics_spaceranger",
        "probe_set",
        "probe_set_meta",
    ]
    assert list(adata.obsm.keys()) == ["spatial"]
    assert list(adata.var.keys()) == ["gene_symbol", "feature_types", "genome"]

    assert adata.X.dtype.kind == "f"
    assert all(adata.var["feature_types"] == "Gene Expression")
    assert adata.obsm["spatial"].dtype == "float"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
