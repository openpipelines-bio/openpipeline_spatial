import mudata as mu
import squidpy as sq

## VIASH START
par = {
    "input": "resources_test/xenium/xenium_tiny_qc_graph.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "graph_key": "spatial_connectivities",
    "mode": "moran",
    "genes": None,
    "attr": "X",
    "n_perms": 100,
    "layer": None,
    "use_raw": False,
}
## VIASH END


def main():
    print("Reading input MuData...", flush=True)
    mdata = mu.read_h5mu(par["input"])
    adata = mdata.mod[par["modality"]]

    # Check for connectivity key
    if par["graph_key"] not in adata.obsp:
        raise ValueError(
            f"Connectivity key '{par['graph_key']}' not found in .obsp of modality '{par['modality']}'."
        )

    genes = par["genes"]
    if genes and len(genes) == 0:
        genes = None

    # Handle layer
    layer = par["layer"]
    if layer == "None":
        layer = None

    print(f"Calculating spatial autocorrelation ({par['mode']})...", flush=True)

    # Run Squidpy spatial_autocorr
    # Note: sq.gr.spatial_autocorr modifies adata in-place, adding results to .uns
    sq.gr.spatial_autocorr(
        adata,
        connectivity_key=par["graph_key"],
        genes=genes,
        mode=par["mode"],
        attr=par["attr"],
        n_perms=par["n_perms"],
        layer=layer,
        use_raw=par["use_raw"],
    )

    result_key = f"{par['mode']}I" if par["mode"] == "moran" else f"{par['mode']}C"

    if result_key in adata.uns:
        # Log top spatially variable genes
        df = adata.uns[result_key]
        if not df.empty:
            sort_col = "I" if par["mode"] == "moran" else "C"
            print("Top spatially variable genes:", flush=True)
            print(df.sort_values(by=sort_col, ascending=False).head(), flush=True)
    else:
        print(
            f"Warning: Expected key '{result_key}' not found in .uns after calculation.",
            flush=True,
        )

    print("Writing output...", flush=True)
    mdata.write_h5mu(par["output"])
    print("Done!", flush=True)


if __name__ == "__main__":
    main()
