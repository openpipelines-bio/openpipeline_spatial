import sys
from pathlib import Path

import anndata as ad
import spatialdata as sd

## VIASH START
par = {
    "inputs": ["input1.zarr", "input2.zarr"],
    "output": "output.zarr",
    "attrs_merge": None,
    "anndata_join": "inner",
    "anndata_merge": None,
    "anndata_uns_merge": None,
    "anndata_label": None,
    "anndata_pairwise": False,
}
## VIASH END

sys.path.insert(0, meta["resources_dir"])
from setup_logger import setup_logger  # noqa: E402


def get_unique_names_from_inputs(inputs):
    """Derive unique dataset names from input file stems.

    Duplicate stems are disambiguated with numeric suffixes (_2, _3, …).
    """
    names = []
    names_seen = set()
    for input_path in inputs:
        name = Path(input_path).stem
        if name in names_seen:
            count = 2
            candidate = f"{name}_{count}"
            while candidate in names_seen:
                count += 1
                candidate = f"{name}_{count}"
            name = candidate
        names.append(name)
        names_seen.add(name)
    return names


def _find_renamed_region(original_region, sample_key, concatenated):
    all_elements = (
        set(concatenated.shapes) | set(concatenated.labels) | set(concatenated.points)
    )
    candidate = f"{original_region}-{sample_key}"
    if candidate in all_elements:
        return candidate
    if original_region in all_elements:
        return original_region
    raise ValueError(
        f"Cannot find region '{original_region}' for sample '{sample_key}' "
        f"in concatenated SpatialData. Available elements: {sorted(all_elements)}"
    )


def main(par):
    logger = setup_logger()
    logger.info("Concatenate SpatialData (spatialdata v%s)", sd.__version__)

    if len(par["inputs"]) == 1:
        logger.warning(
            "Only one input provided – writing to output without concatenation."
        )
        logger.info("Reading single input '%s'", par["inputs"][0])
        name = get_unique_names_from_inputs(par["inputs"])[0]
        sdata = sd.read_zarr(par["inputs"][0])

        logger.info("Applying table merge settings via anndata.concat…")
        original_attrs = sdata["table"].uns.get("spatialdata_attrs")
        sdata["table"] = ad.concat(
            {name: sdata["table"]},
            join=par["anndata_join"],
            merge=par["anndata_merge"],
            uns_merge=par["anndata_uns_merge"],
            label=par["anndata_label"],
            index_unique="-",
            pairwise=par["anndata_pairwise"],
        )
        if original_attrs is not None:
            sdata["table"].uns["spatialdata_attrs"] = original_attrs

        logger.info("Writing SpatialData to '%s'", par["output"])
        sdata.write(par["output"])
        logger.info("Done.")
        return 0

    logger.info("Reading %d SpatialData objects…", len(par["inputs"]))
    names = get_unique_names_from_inputs(par["inputs"])
    sdatas = {}
    for name, input_path in zip(names, par["inputs"]):
        logger.info("  %s  ←  %s", name, input_path)
        sdatas[name] = sd.read_zarr(input_path)

    logger.info("Concatenating SpatialData objects…")
    concatenated = sd.concatenate(
        sdatas,
        concatenate_tables=False,
        attrs_merge=par["attrs_merge"],
    )

    logger.info("Concatenating main tables…")
    region_key = None
    instance_key = None
    all_new_regions = []
    updated_tables = {}

    for name, sdata_item in sdatas.items():
        table = sdata_item["table"]
        attrs = table.uns.get("spatialdata_attrs")

        if attrs is None:
            updated_tables[name] = table
            continue

        if region_key is None:
            region_key = attrs.get("region_key")
        if instance_key is None:
            instance_key = attrs.get("instance_key")

        original_regions = attrs.get("region", [])
        if hasattr(original_regions, "tolist"):
            original_regions = original_regions.tolist()
        elif isinstance(original_regions, str):
            original_regions = [original_regions]
        else:
            original_regions = list(original_regions)

        region_mapping = {
            r: _find_renamed_region(r, name, concatenated) for r in original_regions
        }
        all_new_regions.extend(region_mapping.values())

        if region_key and region_key in table.obs.columns:
            table.obs[region_key] = table.obs[region_key].map(
                lambda x, m=region_mapping: m.get(x, x)
            )
        updated_tables[name] = table

    concatenated["table"] = ad.concat(
        updated_tables,
        join=par["anndata_join"],
        merge=par["anndata_merge"],
        uns_merge=par["anndata_uns_merge"],
        label=par["anndata_label"],
        index_unique="-",
        pairwise=par["anndata_pairwise"],
    )

    if region_key and instance_key and all_new_regions:
        concatenated["table"].uns["spatialdata_attrs"] = {
            "region": all_new_regions,
            "region_key": region_key,
            "instance_key": instance_key,
        }

    logger.info("Writing concatenated SpatialData to '%s'…", par["output"])
    concatenated.write(par["output"])
    logger.info("Done.")


if __name__ == "__main__":
    sys.exit(main(par))
