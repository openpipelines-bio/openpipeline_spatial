vvvvvvvvvvvvimport re
import subprocess
import sys
import pytest
import numpy as np
import spatialdata as sd
from spatialdata.models import Image2DModel
from cellpose.models import CellposeModel

## VIASH START
meta = {
    "executable": "./target/executable/segmentation/cellpose_sam/cellpose_sam",
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
    # path (`--pretrained_model_file`) that a real user-trained model would.
    return CellposeModel(gpu=False, pretrained_model="cpsam_v2").pretrained_model


@pytest.fixture(scope="module")
def small_input_file(tmp_path_factory):
    # Cellpose-SAM tiles its input into fixed 256x256 patches and runs a full
    # transformer forward pass per tile, so segmenting the full ~3500x5800px
    # xenium_tiny image test file (300+ tiles) takes on the order of hours on a
    # CPU-only CI runner. Crop to a small, still cell-containing region so
    # the tests that actually run inference stay within the CI time budget
    # (this crop takes ~30s on CPU).
    sdata = sd.read_zarr(input_file)
    image_arr = _get_image_array(sdata.images["morphology_focus"])
    crop = image_arr[:, 1536:2048, 4096:4608]

    cropped_sdata = sd.SpatialData(
        images={"morphology_focus": Image2DModel.parse(crop, dims=("c", "y", "x"))}
    )
    path = tmp_path_factory.mktemp("small_input") / "small_input.zarr"
    cropped_sdata.write(path)
    return str(path)


def test_default_execution(run_component, tmp_path, small_input_file):
    output = tmp_path / "segmented.zarr"

    run_component(
        [
            "--input",
            small_input_file,
            "--output",
            str(output),
        ]
    )

    assert output.is_dir(), "Output Zarr store was not created."
    sdata = sd.read_zarr(output)

    assert "cellpose_sam_labels" in sdata.labels, (
        "Expected default output labels key to be present."
    )

    labels_arr = np.asarray(sdata.labels["cellpose_sam_labels"].data)
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


def test_custom_output_labels_and_channels(run_component, tmp_path, small_input_file):
    # Combined into a single run (rather than separate tests per option) to
    # avoid paying for Cellpose-SAM's (comparatively heavy, transformer-based)
    # inference more than once per behavior under test.
    output = tmp_path / "segmented_custom.zarr"

    run_component(
        [
            "--input",
            small_input_file,
            "--output",
            str(output),
            "--output_labels",
            "nuclei_masks",
            "--channels",
            "0",
        ]
    )

    sdata = sd.read_zarr(output)
    assert "nuclei_masks" in sdata.labels, (
        "Expected custom output labels key to be present."
    )
    assert "cellpose_sam_labels" not in sdata.labels, (
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


def test_pretrained_model_file_takes_precedence_over_pretrained_model_name(
    run_component, tmp_path, pretrained_model_path, small_input_file
):
    # Also covers plain `--pretrained_model_file` usage (custom pretrained
    # model loading + successful segmentation): combined with the precedence
    # check into a single run to avoid paying for a second inference pass.
    output = tmp_path / "segmented_precedence.zarr"

    stdout = run_component(
        [
            "--input",
            small_input_file,
            "--output",
            str(output),
            "--pretrained_model_file",
            pretrained_model_path,
            # A pretrained_model_name other than the component default, so
            # that if this branch were (incorrectly) taken instead, it would
            # be visible in the logs.
            "--pretrained_model_name",
            "cpsam",
        ]
    ).decode("utf-8")

    assert "Loading custom pretrained model" in stdout, (
        "Expected the component to report that it is using the custom pretrained model."
    )
    assert "Loading built-in model" not in stdout, (
        "'--pretrained_model_file' should take precedence over "
        "'--pretrained_model_name', per the documented behavior of these "
        "two arguments."
    )

    assert output.is_dir(), "Output Zarr store was not created."
    sdata = sd.read_zarr(output)
    assert "cellpose_sam_labels" in sdata.labels

    labels_arr = np.asarray(sdata.labels["cellpose_sam_labels"].data)
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
                "--pretrained_model_file",
                str(bad_model_file),
            ]
        )

    assert not output.exists(), (
        "No (partial) output should be written when the pretrained model "
        "file fails to load."
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
