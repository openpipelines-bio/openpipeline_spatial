import os
import sys
import pytest
import pandas as pd


def test_simple_execution(run_component, tmp_path):
    output_path = tmp_path / "output"
    dataset_id = "Lung5_Rep2"
    run_component(
        [
            "--input",
            meta["resources_dir"] + "/Lung5_Rep2_tiny",
            "--dataset_id",
            dataset_id,
            "--num_fovs",
            "2",
            "--output",
            output_path,
        ]
    )

    assert os.path.exists(output_path), "Output folder was not created"

    counts_file = output_path / f"{dataset_id}_exprMat_file.csv"
    fov_file = output_path / f"{dataset_id}_fov_positions_file.csv"
    meta_file = output_path / f"{dataset_id}_metadata_file.csv"
    tx_file = output_path / f"{dataset_id}_tx_file.csv"

    matrices = [counts_file, fov_file, meta_file, tx_file]
    images = ["CellComposite", "CellLabels", "CellOverlay", "CompartmentLabels"]

    for image in images:
        assert os.path.exists(output_path / image), f"{image} folder was not created"
        assert len(os.listdir(output_path / image)) == 2, (
            f"{image} folder should contain 2 files"
        )

    for matrix in matrices:
        assert os.path.exists(matrix), f"{matrix} file was not created"
        data = pd.read_csv(matrix)
        data["fov"].value_counts().shape[0] == 2, f"{matrix} should contain 2 fovs"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
