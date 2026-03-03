import mudata as mu
import sys
import pytest

## VIASH START
meta = {
    "resources_dir": "resources_test/xenium",
    "executable": "./xenium_spatial_statistics",
}
## VIASH END


def test_calculate_spatial_statistics(run_component, tmp_path):
    input_path = meta["resources_dir"] + "/xenium_tiny.qc.neighbors.h5mu"
    output_path = tmp_path / "output.h5mu"

    # Run the component directly on the input file
    run_component(
        [
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--modality",
            "rna",
            "--obsm_spatial_coordinates",
            "spatial",
            "--density_bandwidth",
            "50.0",
            "--calculate_ripley_l",
            "false",
        ]
    )

    assert output_path.exists(), "Output file not created"

    mdata_out = mu.read_h5mu(output_path)
    adata_out = mdata_out.mod["rna"]

    # Check for expected columns in .obs
    expected_cols = [
        "spatial_nucleus_cell_ratio",
        "spatial_cell_area_percentile",
        "spatial_distance_to_centroid",
        "spatial_norm_x",
        "spatial_norm_y",
        "spatial_kernel_density",
        "spatial_graph_degree",
        "spatial_voronoi_area",
        "spatial_voronoi_neighbors",
    ]

    for col in expected_cols:
        assert col in adata_out.obs.columns, f"Column {col} missing from output .obs"

    # Check for expected metrics in .uns['spatial_stats']
    assert "spatial_stats" in adata_out.uns, "spatial_stats missing from .uns"
    stats = adata_out.uns["spatial_stats"]

    expected_stats = [
        "cell_density",
        "total_area",
        "n_cells",
        "spatial_extent_x",
        "spatial_extent_y",
        "centroid_x",
        "centroid_y",
    ]

    for stat in expected_stats:
        assert stat in stats, f"Stat {stat} missing from spatial_stats"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
