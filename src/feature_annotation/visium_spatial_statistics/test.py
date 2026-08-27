import pytest
import sys
import mudata as md

## VIASH START
meta = {
    "resources_dir": "resources_test/visium",
    "executable": "./visium_spatial_statistics",
}
## VIASH END


@pytest.fixture
def input_path(meta):
    return meta["resources_dir"] + "/visium_tiny_neighbours.h5mu"


@pytest.fixture
def output_path(tmp_path):
    return str(tmp_path / "output.h5mu")


def test_default(run_component, input_path, output_path):
    run_component(
        [
            "--input",
            input_path,
            "--output",
            output_path,
        ]
    )

    result = md.read_h5mu(output_path)
    adata = result["rna"]

    # --- .obs columns ---
    for col in [
        "spatial_n_neighbors",
        "spatial_tissue_edge",
        "spatial_distance_to_centroid",
        "spatial_norm_x",
        "spatial_norm_y",
        "spatial_distance_to_boundary",
        "spatial_local_expression_density",
    ]:
        assert col in adata.obs.columns, f"Missing obs column: {col}"

    # n_neighbors: between 0 and 6 (hexagonal grid)
    n_neighbors = adata.obs["spatial_n_neighbors"]
    assert (n_neighbors >= 0).all()
    assert (n_neighbors <= 6).all()

    # tissue_edge: bool, at least some peripheral spots
    edge = adata.obs["spatial_tissue_edge"]
    assert edge.dtype == bool
    assert edge.any(), "Expected at least some peripheral spots"

    # distance_to_centroid: non-negative
    assert (adata.obs["spatial_distance_to_centroid"] >= 0).all()

    # normalized coordinates: within [0, 1]
    assert (adata.obs["spatial_norm_x"] >= 0).all()
    assert (adata.obs["spatial_norm_x"] <= 1).all()
    assert (adata.obs["spatial_norm_y"] >= 0).all()
    assert (adata.obs["spatial_norm_y"] <= 1).all()

    # distance_to_boundary: non-negative
    assert (adata.obs["spatial_distance_to_boundary"] >= 0).all()

    # local_expression_density: non-negative (sum of neighbour counts)
    assert (adata.obs["spatial_local_expression_density"] >= 0).all()

    # --- .uns["spatial_stats"] ---
    assert "spatial_stats" in adata.uns
    stats = adata.uns["spatial_stats"]

    assert stats["area_calculation_method"] == "convex_hull"
    assert stats["spot_density"] > 0
    assert stats["total_area"] > 0
    assert stats["n_spots"] == len(adata)
    assert stats["spatial_extent_x"] > 0
    assert stats["spatial_extent_y"] > 0
    assert "centroid_x" in stats
    assert "centroid_y" in stats


def test_tissue_edge_max_neighbors(run_component, input_path, output_path):
    """A higher neighbour threshold classifies more spots as tissue edge.

    The default threshold (6) suits the hexagonal Visium grid; Visium HD uses 8
    for its square-bin (Moore) neighbourhood. On this hexagonal resource (max 6
    neighbours) a threshold of 8 marks every spot as tissue edge.
    """
    run_component(
        [
            "--input",
            input_path,
            "--output",
            output_path,
            "--tissue_edge_max_neighbors",
            "8",
        ]
    )

    adata = md.read_h5mu(output_path)["rna"]
    edge = adata.obs["spatial_tissue_edge"]
    n_neighbors = adata.obs["spatial_n_neighbors"]

    assert edge.dtype == bool
    # threshold 8 > hexagonal max of 6, so all spots are below threshold
    assert (edge == (n_neighbors < 8)).all()
    assert edge.all(), "With threshold 8 every hexagonal-grid spot is tissue edge"


def test_custom_prefix(run_component, input_path, output_path):
    run_component(
        [
            "--input",
            input_path,
            "--output",
            output_path,
            "--output_prefix",
            "vis_",
        ]
    )

    adata = md.read_h5mu(output_path)["rna"]
    for col in [
        "vis_n_neighbors",
        "vis_tissue_edge",
        "vis_distance_to_centroid",
        "vis_norm_x",
        "vis_norm_y",
        "vis_distance_to_boundary",
        "vis_local_expression_density",
    ]:
        assert col in adata.obs.columns, f"Missing obs column: {col}"


def test_custom_obs_total_counts(run_component, input_path, output_path, tmp_path):
    # Create a copy of the input with total_counts renamed to n_counts
    mdata = md.read_h5mu(input_path)
    mdata["rna"].obs["n_counts"] = mdata["rna"].obs["total_counts"]
    del mdata["rna"].obs["total_counts"]
    custom_input = str(tmp_path / "custom_input.h5mu")
    mdata.write_h5mu(custom_input)

    run_component(
        [
            "--input",
            custom_input,
            "--output",
            output_path,
            "--obs_total_counts",
            "n_counts",
        ]
    )

    adata = md.read_h5mu(output_path)["rna"]
    assert "spatial_local_expression_density" in adata.obs.columns, (
        "local_expression_density should be present when custom column is found"
    )
    assert (adata.obs["spatial_local_expression_density"] >= 0).all()


def test_missing_obs_total_counts(run_component, input_path, output_path, tmp_path):
    # Create a copy of the input without total_counts
    mdata = md.read_h5mu(input_path)
    if "total_counts" in mdata["rna"].obs.columns:
        del mdata["rna"].obs["total_counts"]
    no_counts_input = str(tmp_path / "no_counts_input.h5mu")
    mdata.write_h5mu(no_counts_input)

    # Should complete without error, skipping local_expression_density
    run_component(
        [
            "--input",
            no_counts_input,
            "--output",
            output_path,
        ]
    )

    adata = md.read_h5mu(output_path)["rna"]
    assert "spatial_local_expression_density" not in adata.obs.columns, (
        "local_expression_density should be absent when total_counts is missing"
    )


def test_custom_obsm_spatial_coordinates(
    run_component, input_path, output_path, tmp_path
):
    # Rename .obsm["spatial"] to a custom key
    mdata = md.read_h5mu(input_path)
    mdata["rna"].obsm["my_coords"] = mdata["rna"].obsm["spatial"]
    del mdata["rna"].obsm["spatial"]
    custom_input = str(tmp_path / "custom_coords_input.h5mu")
    mdata.write_h5mu(custom_input)

    run_component(
        [
            "--input",
            custom_input,
            "--output",
            output_path,
            "--obsm_spatial_coordinates",
            "my_coords",
        ]
    )

    adata = md.read_h5mu(output_path)["rna"]
    for col in [
        "spatial_n_neighbors",
        "spatial_tissue_edge",
        "spatial_distance_to_centroid",
        "spatial_norm_x",
        "spatial_norm_y",
        "spatial_distance_to_boundary",
    ]:
        assert col in adata.obs.columns, f"Missing obs column: {col}"


def test_custom_obsp_spatial_graph(run_component, input_path, output_path, tmp_path):
    # Rename .obsp["spatial_connectivities"] to a custom key
    mdata = md.read_h5mu(input_path)
    mdata["rna"].obsp["my_graph"] = mdata["rna"].obsp["spatial_connectivities"]
    del mdata["rna"].obsp["spatial_connectivities"]
    custom_input = str(tmp_path / "custom_graph_input.h5mu")
    mdata.write_h5mu(custom_input)

    run_component(
        [
            "--input",
            custom_input,
            "--output",
            output_path,
            "--obsp_spatial_graph",
            "my_graph",
        ]
    )

    adata = md.read_h5mu(output_path)["rna"]
    for col in [
        "spatial_n_neighbors",
        "spatial_tissue_edge",
        "spatial_local_expression_density",
    ]:
        assert col in adata.obs.columns, f"Missing obs column: {col}"

    # n_neighbors should still be in valid range for hexagonal grid
    n_neighbors = adata.obs["spatial_n_neighbors"]
    assert (n_neighbors >= 0).all()
    assert (n_neighbors <= 6).all()


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
