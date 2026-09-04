import re
import subprocess
import sys
import pytest
import numpy as np
import spatialdata as sd
from cellpose.models import CellposeModel

## VIASH START
meta = {
    "executable": "./target/executable/segmentation/cellpose3/cellpose3",
    "resources_dir": "resources_test/xenium/",
}
## VIASH END

input_file = f"{meta['resources_dir']}/xenium_tiny.zarr"


def _get_image_array(image_element):
    # Multiscale images (DataTree) expose full-resolution pixel data under
    # the "scale0" node; single-scale images expose it directly via .data.
    if hasattr(image_element, "data"):
        return np.asarray(image_element.data)
    return np.asarray(image_element["scale0"]["image"].data)


@pytest.fixture(scope="module")
def pretrained_model_path():
    # Reuse a cached built-in model's checkpoint file as a stand-in for a
    # "custom" pretrained model, rather than shipping a separate model file
    # as a test resource: it goes through the exact same file-based loading
    # path (`--pretrained_model`) that a real user-trained model would.
    return CellposeModel(gpu=False, model_type="nuclei").pretrained_model


def test_default_execution(run_component, tmp_path):
    output = tmp_path / "segmented.zarr"

    run_component(
        [
            "--input",
            input_file,
            "--output",
            str(output),
        ]
    )

    assert output.is_dir(), "Output Zarr store was not created."
    sdata = sd.read_zarr(output)

    assert "cellpose_labels" in sdata.labels, (
        "Expected default output labels key to be present."
    )

    labels_arr = np.asarray(sdata.labels["cellpose_labels"].data)
    image_arr = _get_image_array(sdata.images["morphology_focus"])

    assert labels_arr.shape == image_arr.shape[-2:], (
        "Labels shape should match the (y, x) shape of the input image."
    )
    assert labels_arr.dtype.kind in ("u", "i"), (
        "Labels should be stored as an integer array."
    )
    n_objects = len(np.unique(labels_arr)) - 1
    assert n_objects > 0, (
        "Expected at least one segmented object with the default settings "
        "(the default normalization percentiles matter here: Cellpose's own "
        "default of [1, 99] detects 0 objects on this mostly-background "
        "fluorescence image)."
    )


def test_custom_output_labels(run_component, tmp_path):
    output = tmp_path / "segmented_custom.zarr"

    run_component(
        [
            "--input",
            input_file,
            "--output",
            str(output),
            "--output_labels",
            "nuclei_masks",
        ]
    )

    sdata = sd.read_zarr(output)
    assert "nuclei_masks" in sdata.labels, (
        "Expected custom output labels key to be present."
    )
    assert "cellpose_labels" not in sdata.labels, (
        "Default labels key should not be present when a custom key is used."
    )


def test_fail_missing_image_key(run_component, tmp_path):
    output = tmp_path / "should_not_exist.zarr"

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_file,
                "--output",
                str(output),
                "--input_image",
                "nonexistent_image",
            ]
        )
    assert re.search(
        r"Image key 'nonexistent_image' not found",
        err.value.stdout.decode("utf-8"),
    )


def test_fail_invalid_normalize_percentiles(run_component, tmp_path):
    output = tmp_path / "should_not_exist.zarr"

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_file,
                "--output",
                str(output),
                "--normalize_percentile_low",
                "99.9",
                "--normalize_percentile_high",
                "1.0",
            ]
        )
    assert re.search(
        r"'--normalize_percentile_low' \(99.9\) must be lower than "
        r"'--normalize_percentile_high' \(1.0\)",
        err.value.stdout.decode("utf-8"),
    )


def test_pretrained_model_file(run_component, tmp_path, pretrained_model_path):
    output = tmp_path / "segmented_pretrained.zarr"

    stdout = run_component(
        [
            "--input",
            input_file,
            "--output",
            str(output),
            "--pretrained_model",
            pretrained_model_path,
        ]
    ).decode("utf-8")

    assert "Loading custom pretrained model" in stdout, (
        "Expected the component to report that it is using the custom "
        "pretrained model."
    )

    assert output.is_dir(), "Output Zarr store was not created."
    sdata = sd.read_zarr(output)
    assert "cellpose_labels" in sdata.labels

    labels_arr = np.asarray(sdata.labels["cellpose_labels"].data)
    n_objects = len(np.unique(labels_arr)) - 1
    assert n_objects > 0, (
        "Expected at least one segmented object using the custom pretrained model."
    )


def test_pretrained_model_takes_precedence_over_model_type(
    run_component, tmp_path, pretrained_model_path
):
    output = tmp_path / "segmented_precedence.zarr"

    stdout = run_component(
        [
            "--input",
            input_file,
            "--output",
            str(output),
            "--pretrained_model",
            pretrained_model_path,
            # A model_type other than the component default, so that if this
            # branch were (incorrectly) taken instead, it would be visible in
            # the logs.
            "--model_type",
            "cyto3",
        ]
    ).decode("utf-8")

    assert "Loading custom pretrained model" in stdout
    assert "Loading built-in model" not in stdout, (
        "'--pretrained_model' should take precedence over '--model_type', per "
        "the documented behavior of these two arguments."
    )


def test_fail_invalid_pretrained_model_file(run_component, tmp_path):
    bad_model_file = tmp_path / "not_a_model.pt"
    bad_model_file.write_bytes(b"this is not a valid cellpose model file")
    output = tmp_path / "should_not_exist.zarr"

    with pytest.raises(subprocess.CalledProcessError):
        run_component(
            [
                "--input",
                input_file,
                "--output",
                str(output),
                "--pretrained_model",
                str(bad_model_file),
            ]
        )

    assert not output.exists(), (
        "No (partial) output should be written when the pretrained model "
        "file fails to load."
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
