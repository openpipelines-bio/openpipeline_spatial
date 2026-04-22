import sys

################################################################################
# VIASH
################################################################################

## VIASH START
par = {
    "inputs": ["input1.zarr", "input2.zarr"],
    "output": "output.zarr",
    "attrs_merge": "same",
    "anndata_join": "inner",
    "anndata_merge": "same",
    "anndata_uns_merge": "same",
    "anndata_label": None,
    "anndata_pairwise": False,
}
## VIASH END


def main(par):
    # Placeholder scaffold. Logic will be ported in the next phase.
    print(
        "concatenate_spatialdata scaffold created; implementation pending.",
        flush=True,
    )
    print(f"inputs={len(par['inputs'])}, output={par['output']}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(par))
