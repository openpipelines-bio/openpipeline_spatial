import re
import subprocess
import sys
import pytest
import numpy as np
import spatialdata as sd

## VIASH START
meta = {
    "executable": "./target/executable/segmentation/cellposeSAM/cellposeSAM",
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

    assert "cellposeSAM_labels" in sdata.labels, (
        "Expected default output labels key to be present."
    )

    labels_arr = np.asarray(sdata.labels["cellposeSAM_labels"].data)
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


def test_custom_output_labels_and_channels(run_component, tmp_path):
    # Combined into a single run (rather than separate tests per option) to
    # avoid paying for Cellpose-SAM's (comparatively heavy, transformer-based)
    # inference more than once per behavior under test.
    output = tmp_path / "segmented_custom.zarr"

    run_component(
        [
            "--input",
            input_file,
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
    assert "cellposeSAM_labels" not in sdata.labels, (
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


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
