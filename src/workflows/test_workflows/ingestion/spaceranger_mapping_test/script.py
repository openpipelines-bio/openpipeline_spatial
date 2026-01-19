from mudata import read_h5mu
import sys
import pytest

##VIASH START
par = {"input": "input.h5mu"}

meta = {"resources_dir": "resources_test"}
##VIASH END


def test_run():
    input_mudata = read_h5mu(par["input"])
    expected_var_columns = ["gene_symbol", "feature_types", "genome"]

    assert list(input_mudata.mod.keys()) == ["rna"], (
        "Input should contain rna modality."
    )
    assert list(input_mudata.var.columns) == expected_var_columns, (
        f"Input var columns should be: {expected_var_columns}."
    )
    assert list(input_mudata.mod["rna"].var.columns) == expected_var_columns, (
        f"Input mod['rna'] var columns should be: {expected_var_columns}."
    )

    assert list(input_mudata.mod["rna"].obs.keys()) == ["spatial"]


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "--import-mode=importlib"]))
