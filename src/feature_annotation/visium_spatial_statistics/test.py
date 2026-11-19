import sys

import mudata as md
import pytest

## VIASH START
meta = {
    "executable": "./target/executable/feature_annotation/visium_spatial_statistics/visium_spatial_statistics",
    "resources_dir": "resources_test/visium",
    "config": "src/feature_annotation/visium_spatial_statistics/config.vsh.yaml",
}
## VIASH END

input = f"{meta['resources_dir']}/visium_tiny_neighbours.h5mu"


def test_default(run_component, tmp_path):
    output = tmp_path / "output.h5mu"
    run_component(["--input", input, "--output", str(output)])

    assert output.is_file(), "output file was not created"

    adata = md.read_h5mu(output)["rna"]

    # .obs columns
    for col in [
        "spatial_n_neighbors",
        "spatial_tissue_edge",
        "spatial_distance_to_centroid",
        "spatial_norm_x",
        "spatial_norm_y",
        "spatial_distance_to_boundary",
        "spatial_local_expression_density",
    ]:
        assert col in adata.obs.columns, f"missing .obs column: {col}"

    # n_neighbors is between 0 and 6 (hexagonal grid)
    n_neighbors = adata.obs["spatial_n_neighbors"]
    assert (n_neighbors >= 0).all()
    assert (n_neighbors <= 6).all()

    # tissue_edge is boolean, with at least some peripheral spots
    edge = adata.obs["spatial_tissue_edge"]
    assert edge.dtype == bool
    assert edge.any(), "expected at least some peripheral spots"

    assert (adata.obs["spatial_distance_to_centroid"] >= 0).all()

    # Normalised coordinates are in [0, 1]
    for col in ["spatial_norm_x", "spatial_norm_y"]:
        assert (adata.obs[col] >= 0).all()
        assert (adata.obs[col] <= 1).all()

    assert (adata.obs["spatial_distance_to_boundary"] >= 0).all()

    # local_expression_density is a sum of neighbour counts
    assert (adata.obs["spatial_local_expression_density"] >= 0).all()

    # .uns["spatial_stats"]
    assert "spatial_stats" in adata.uns, "missing .uns['spatial_stats']"
    stats = adata.uns["spatial_stats"]

    assert stats["area_calculation_method"] == "convex_hull"
    assert stats["spot_density"] > 0
    assert stats["total_area"] > 0
    assert stats["n_spots"] == adata.n_obs
    assert stats["spatial_extent_x"] > 0
    assert stats["spatial_extent_y"] > 0
    assert "centroid_x" in stats
    assert "centroid_y" in stats


def test_tissue_edge_max_neighbors(run_component, tmp_path):
    """A higher neighbour threshold classifies more spots as tissue edge.

    The default threshold (6) suits the hexagonal Visium grid. Visium HD uses 8
    for its square-bin (Moore) neighbourhood. On this hexagonal resource (at
    most 6 neighbours) a threshold of 8 marks every spot as tissue edge.
    """
    output = tmp_path / "output.h5mu"
    run_component(
        [
            "--input",
            input,
            "--output",
            str(output),
            "--tissue_edge_max_neighbors",
            "8",
        ]
    )

    adata = md.read_h5mu(output)["rna"]
    edge = adata.obs["spatial_tissue_edge"]
    n_neighbors = adata.obs["spatial_n_neighbors"]

    assert edge.dtype == bool
    assert (edge == (n_neighbors < 8)).all()
    assert edge.all(), "with threshold 8 every hexagonal-grid spot is tissue edge"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
