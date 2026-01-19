import sys
import pytest
from pathlib import Path

## VIASH START
meta = {"name": "spaceranger_count", "resources_dir": "resources_test"}
## VIASH END

input = meta["resources_dir"] + "/visium/Visium_FFPE_Human_Ovarian_Cancer_tiny/"
probe_set = meta["resources_dir"] + "/visium/Visium_FFPE_Human_Ovarian_Cancer_probe_set.csv"
image = meta["resources_dir"] + "/visium/Visium_FFPE_Human_Ovarian_Cancer_image_tiny.jpg"
reference = meta["resources_dir"] + "/GRCh38"


def test_simple_execution(run_component, random_path):
    output = random_path()
    run_component(
        [
            "--input",
            input,
            "--gex_reference",
            reference,
            "--probe_set",
            probe_set,
            "--image",
            image,
            "--area",
            "D1",
            "--slide",
            "V10L13-020",
            "--create_bam",
            "false",
            "--output",
            output,
        ]
    )

    assert (output / "filtered_feature_bc_matrix.h5").is_file(), (
        "No filtered .h5 count matrix was created."
    )

    assert (output / "raw_feature_bc_matrix.h5").is_file(), (
        "No raw .h5 count matrix was created."
    )

    assert (output / "metrics_summary.csv").is_file(), (
        "No metrics summary was created."
    )

    assert (output / "web_summary.html").is_file(), (
        "No web summary was created."
    )


def test_with_fastqs(run_component, random_path):
    output = random_path()
    run_component(
        [
            "--input",
            Path(meta["resources_dir"] + "/visium/Visium_FFPE_Human_Ovarian_Cancer_tiny/Visium_FFPE_Human_Ovarian_Cancer_S1_L001_R1_001.fastq.gz"),
            "--input",
            Path(meta["resources_dir"] + "/visium/Visium_FFPE_Human_Ovarian_Cancer_tiny/Visium_FFPE_Human_Ovarian_Cancer_S1_L001_R2_001.fastq.gz"),
            "--gex_reference",
            reference,
            "--probe_set",
            probe_set,
            "--image",
            image,
            "--area",
            "D1",
            "--slide",
            "V10L13-020",
            "--create_bam",
            "false",
            "--output",
            output,
        ]
    )

    assert (output / "filtered_feature_bc_matrix.h5").is_file(), (
        "No filtered .h5 count matrix was created."
    )

    assert (output / "raw_feature_bc_matrix.h5").is_file(), (
        "No raw .h5 count matrix was created."
    )

    assert (output / "metrics_summary.csv").is_file(), (
        "No metrics summary was created."
    )

    assert (output / "web_summary.html").is_file(), (
        "No web summary was created."
    )


def test_with_optional_params(run_component, random_path):
    output = random_path()
    run_component(
        [
            "--input",
            input,
            "--gex_reference",
            reference,
            "--probe_set",
            probe_set,
            "--image",
            image,
            "--area",
            "D1",
            "--slide",
            "V10L13-020",
            "--nosecondary",
            "true",
            "--r1_length",
            "100",
            "--r2_length",
            "100",
            "filter_probes",
            "false",
            "--create_bam",
            "true",
            "--output",
            output,
        ]
    )

    assert (output / "filtered_feature_bc_matrix.h5").is_file(), (
        "No filtered .h5 count matrix was created."
    )

    assert (output / "raw_feature_bc_matrix.h5").is_file(), (
        "No raw .h5 count matrix was created."
    )

    assert (output / "metrics_summary.csv").is_file(), (
        "No metrics summary was created."
    )

    assert (output / "web_summary.html").is_file(), (
        "No web summary was created."
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
