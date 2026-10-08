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
    "resources_dir": "resources_test/xenium_multichannel/",
}
## VIASH END

input_file = f"{meta['resources_dir']}/xenium_multicellseg_tiny.zarr"


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


# Each case is one inference run over the 4-channel image, so the other
# argument checks piggyback on them instead of paying for separate runs:
#  * channel 1 also sets a custom '--output_labels' key.
#  * channel 2 also sets '--nuclear_channel' (ignored by the default `nuclei`
#    model, but exercises the two-channel selection in the ndim == 3 branch).
@pytest.mark.parametrize(
    "cytoplasm_channel,nuclear_channel,output_labels",
    [
        (1, 0, "nuclei_masks"),
        (2, 1, "cellpose_labels"),
        (3, 0, "cellpose_labels"),
        (4, 0, "cellpose_labels"),
    ],
)
def test_each_channel_of_multichannel_image(
    run_component, tmp_path, cytoplasm_channel, nuclear_channel, output_labels
):
    output = tmp_path / f"segmented_channel_{cytoplasm_channel}.zarr"

    run_component(
        [
            "--input",
            input_file,
            "--output",
            str(output),
            "--cytoplasm_channel",
            str(cytoplasm_channel),
            "--nuclear_channel",
            str(nuclear_channel),
            "--output_labels",
            output_labels,
        ]
    )

    assert output.is_dir(), "Output Zarr store was not created."
    sdata = sd.read_zarr(output)

    image_arr = _get_image_array(sdata.images["morphology_focus"])
    assert image_arr.ndim == 3 and image_arr.shape[0] == 4, (
        "Expected the multichannel test image to have 4 channels."
    )

    assert output_labels in sdata.labels, (
        f"Expected output labels key '{output_labels}' to be present."
    )
    if output_labels != "cellpose_labels":
        assert "cellpose_labels" not in sdata.labels, (
            "Default labels key should not be present when a custom key is used."
        )

    labels_arr = np.asarray(sdata.labels[output_labels].data)
    assert labels_arr.shape == image_arr.shape[-2:], (
        "Labels shape should match the (y, x) shape of the input image."
    )
    n_objects = len(np.unique(labels_arr)) - 1
    assert n_objects > 0, (
        f"Expected at least one segmented object when selecting channel "
        f"{cytoplasm_channel} of the multichannel image."
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


def test_pretrained_model_file_takes_precedence_over_model_type(
    run_component, tmp_path, pretrained_model_path
):
    # Also covers plain `--pretrained_model` usage (custom pretrained model
    # loading + successful segmentation): combined with the precedence check
    # into a single run to avoid paying for a second inference pass.
    output = tmp_path / "segmented_pretrained.zarr"

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

    assert "Loading custom pretrained model" in stdout, (
        "Expected the component to report that it is using the custom pretrained model."
    )
    assert "Loading built-in model" not in stdout, (
        "'--pretrained_model' should take precedence over '--model_type', per "
        "the documented behavior of these two arguments."
    )

    assert output.is_dir(), "Output Zarr store was not created."
    sdata = sd.read_zarr(output)
    assert "cellpose_labels" in sdata.labels

    labels_arr = np.asarray(sdata.labels["cellpose_labels"].data)
    n_objects = len(np.unique(labels_arr)) - 1
    assert n_objects > 0, (
        "Expected at least one segmented object using the custom pretrained model."
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
