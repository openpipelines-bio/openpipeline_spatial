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


def main(par):
    logger = setup_logger()
    logger.info("Concatenate SpatialData (spatialdata v%s)", sd.__version__)

    if len(par["inputs"]) == 1:
        logger.warning(
            "Only one input provided – writing to output without concatenation."
        )
        logger.info("Reading single input '%s'", par["inputs"][0])
        sdata = sd.read_zarr(par["inputs"][0])
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
    tables = {name: sdata["table"] for name, sdata in sdatas.items()}
    concatenated["table"] = ad.concat(
        tables,
        join=par["anndata_join"],
        merge=par["anndata_merge"],
        uns_merge=par["anndata_uns_merge"],
        label=par["anndata_label"],
        index_unique="-",
        pairwise=par["anndata_pairwise"],
    )

    logger.info("Writing concatenated SpatialData to '%s'…", par["output"])
    concatenated.write(par["output"])
    logger.info("Done.")


if __name__ == "__main__":
    sys.exit(main(par))
