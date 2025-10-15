library(testthat, warn.conflicts = FALSE)
library(hdf5r)


## VIASH START
meta <- list(
  executable = "target/executable/convert/from_h5ad_to_seurat",
  resources_dir = "resources_test_sc",
  name = "from_h5ad_to_seurat"
)
## VIASH END

## Simple conversion
cat("> Checking conversion of single-modality of h5mu file\n")

in_h5mu <- paste0(
  meta[["resources_dir"]],
  "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"
)
out_rds <- "output.rds"

cat("> Running ", meta[["name"]], "\n", sep = "")
out <- processx::run(
  meta[["executable"]],
  c(
    "--input", in_h5mu,
    "--output", out_rds
  )
)

cat("> Checking whether output file exists\n")
expect_equal(out$status, 0)
expect_true(file.exists(out_rds))

cat("> Reading output file\n")
obj <- readRDS(file = out_rds)

cat("> Checking whether Seurat object is in the right format\n")
expect_is(obj, "Seurat")
expect_equal(names(slot(obj, "assays")), "RNA")

dim_rds <- dim(obj)
mu_in <- H5File$new(in_h5mu, mode = "r")
dim_ad <- mu_in[["/mod/rna/X"]]$attr_open("shape")$read()

expect_equal(dim_rds[1], dim_ad[2])
expect_equal(dim_rds[2], dim_ad[1])

## Multi-modal conversion
cat("> Checking conversion of multi-modal h5mu file\n")

in_h5mu <- paste0(
  meta[["resources_dir"]],
  "/10x_5k_anticmv/5k_human_antiCMV_T_TBNK_connect.h5mu"
)
out_rds <- "output.rds"

cat("> Running ", meta[["name"]], "\n", sep = "")
out <- processx::run(
  meta[["executable"]],
  c(
    "--input", in_h5mu,
    "--modality", "rna",
    "--modality", "prot",
    "--modality", "vdj_t",
    "--assay", "RNA",
    "--assay", "ADT",
    "--assay", "TCR",
    "--output", out_rds
  )
)

cat("> Checking whether output file exists\n")
expect_equal(out$status, 0)
expect_true(file.exists(out_rds))

cat("> Reading output file\n")
obj <- readRDS(file = out_rds)

cat("> Checking whether Seurat object is in the right format\n")
expect_is(obj, "Seurat")
expect_true(all(names(slot(obj, "assays")) %in% c("RNA", "ADT")))

dim_rds <- dim(obj)
mu_in <- H5File$new(in_h5mu, mode = "r")
dim_ad <- mu_in[["/mod/rna/X"]]$attr_open("shape")$read()

expect_equal(dim_rds[1], dim_ad[2])
expect_equal(dim_rds[2], dim_ad[1])

vdj_t_cols <- c(
  "TCR_is_cell", "TCR_high_confidence", "TCR_multi_chain", "TCR_extra_chains",
  "TCR_IR_VJ_1_c_call", "TCR_IR_VJ_2_c_call", "TCR_IR_VDJ_1_c_call",
  "TCR_IR_VDJ_2_c_call", "TCR_IR_VJ_1_consensus_count",
  "TCR_IR_VJ_2_consensus_count", "TCR_IR_VDJ_1_consensus_count",
  "TCR_IR_VDJ_2_consensus_count", "TCR_IR_VJ_1_d_call",
  "TCR_IR_VJ_2_d_call", "TCR_IR_VDJ_1_d_call", "TCR_IR_VDJ_2_d_call",
  "TCR_IR_VJ_1_duplicate_count", "TCR_IR_VJ_2_duplicate_count",
  "TCR_IR_VDJ_1_duplicate_count", "TCR_IR_VDJ_2_duplicate_count",
  "TCR_IR_VJ_1_j_call", "TCR_IR_VJ_2_j_call", "TCR_IR_VDJ_1_j_call",
  "TCR_IR_VDJ_2_j_call", "TCR_IR_VJ_1_junction", "TCR_IR_VJ_2_junction",
  "TCR_IR_VDJ_1_junction", "TCR_IR_VDJ_2_junction", "TCR_IR_VJ_1_junction_aa",
  "TCR_IR_VJ_2_junction_aa", "TCR_IR_VDJ_1_junction_aa", "TCR_IR_VDJ_2_junction_aa",
  "TCR_IR_VJ_1_locus", "TCR_IR_VJ_2_locus", "TCR_IR_VDJ_1_locus",
  "TCR_IR_VDJ_2_locus", "TCR_IR_VJ_1_productive", "TCR_IR_VJ_2_productive",
  "TCR_IR_VDJ_1_productive", "TCR_IR_VDJ_2_productive", "TCR_IR_VJ_1_v_call",
  "TCR_IR_VJ_2_v_call", "TCR_IR_VDJ_1_v_call", "TCR_IR_VDJ_2_v_call", "TCR_has_ir"
)
obs_cols <- c("orig.ident", "nCount_RNA", "nFeature_RNA")

expect_true(all(vdj_t_cols %in% colnames(obj@meta.data)))
expect_true(all(obs_cols %in% colnames(obj@meta.data)))