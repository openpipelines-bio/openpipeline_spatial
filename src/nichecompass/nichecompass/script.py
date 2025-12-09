import sys
import json
import mudata as mu

from nichecompass.models import NicheCompass
from nichecompass.utils import add_gps_from_gp_dict_to_adata
from torch.cuda import is_available as cuda_is_available

## VIASH START
par = {
    # Inputs
    "input": "resources_test/cosmx/Lung5_Rep2_tiny.h5mu",
    "modality": "rna",
    "layer": None,
    "input_gp_mask": "resources_test/niche/prior_knowledge_gp_mask.json",
    "input_obs_covariates": None,
    "input_obsm_spatial_connectivities": "spatial_connectivities",
    ## GP Mask
    "min_genes_per_gp": 2,
    "min_source_genes_per_gp": 1,
    "min_target_genes_per_gp": 1,
    "max_genes_per_gp": None,
    "max_source_genes_per_gp": None,
    "max_target_genes_per_gp": None,
    "filter_genes_not_in_masks": False,
    # outputs
    "output": "nichecompass_output.h5mu",
    "output_model": "nichecompass_model/",
    "output_uns_gp_names": "nichecompass_gp_names",
    "output_varm_gp_targets_mask": "nichecompass_gp_targets",
    "output_varm_gp_sources_mask": "nichecompass_gp_sources",
    "output_obsm_embedding": "nichecompass_latent",
    "output_uns_covariate_embeddings": None,
    "output_obsp_reconstructed_adj_edge_proba": "nichecompass_recon_connectivities",
    "output_uns_active_gp_names": "nichecompass_active_gp_names",
    "output_uns_gene_index": "nichecompass_genes_idx",
    "output_uns_target_genes_index": "nichecompass_target_genes_idx",
    "output_uns_source_genes_index": "nichecompass_source_genes_idx",
    "output_obsp_agg_weights": "nichecompass_agg_weights",
    # model architecture
    "include_edge_recon_loss": True,
    "include_gene_expr_recon_loss": True,
    "include_cat_covariates_contrastive_loss": False,
    "covariates_edges": None,
    "covariate_embedding_injection_layers": ["gene_expr_decoder"],
    "gene_expr_recon_dist": "nb",
    "log_variational": True,
    "node_label_method": "one-hop-norm",
    "active_gp_thresh_ratio": 0.01,
    "active_gp_type": "separate",
    "n_fc_layers_encoder": 1,
    "n_layers_encoder": 1,
    "n_hidden_encoder": None,
    "conv_layer_encoder": "gatv2conv",
    "encoder_n_attention_heads": 4,
    "encoder_use_bn": False,
    "dropout_rate_encoder": 0.0,
    "dropout_rate_graph_decoder": 0.0,
    "cat_covariates_cats": None,
    "n_addon_gp": 100,
    "cat_covariates_embeds_nums": None,
    "random_state": 0,
    # model training
    "n_epochs": 1,
    "n_epochs_all_gps": 0,
    "n_epochs_no_edge_recon": 0,
    "n_epochs_no_cat_covariates_contrastive": 0,
    "lr": 0.001,
    "weight_decay": 0.0,
    "lambda_edge_recon": 500000.0,
    "lambda_gene_expr_recon": 300.0,
    "lambda_cat_covariates_contrastive": 0.0,
    "contrastive_logits_pos_ratio": 0.0,
    "contrastive_logits_neg_ratio": 0.0,
    "lambda_group_lasso": 0.0,
    "lambda_l1_masked": 0.0,
    "l1_targets_categories": ["target_gene"],
    "l1_sources_categories": ["source_gene"],
    "lambda_l1_addon": 30.0,
    "edge_val_ratio": 0.1,
    "node_val_ratio": 0.1,
    "edge_batch_size": 256,
    "node_batch_size": None,
    "n_sampled_neighbors": -1,
}

meta = {"resources_dir": "src/utils/"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

use_gpu = cuda_is_available()
logger.info("GPU enabled? %s", use_gpu)

## Read in data
adata = mu.read_h5ad(par["input"], mod=par["modality"])

## Add GP mask to data
logger.info("Adding prior knowledge gene program mask to data...")
with open(par["input_gp_mask"], "r") as f:
    prior_knowledge_gp_mask = json.load(f)

add_gps_from_gp_dict_to_adata(
    gp_dict=prior_knowledge_gp_mask,
    adata=adata,
    gp_targets_mask_key=par["output_varm_gp_targets_mask"],
    gp_sources_mask_key=par["output_varm_gp_sources_mask"],
    gp_names_key=par["output_uns_gp_names"],
    genes_idx_key=par["output_uns_genes_index"],
    target_genes_idx_key=par["output_uns_target_genes_index"],
    source_genes_idx_key=par["output_uns_source_genes_index"],
    min_genes_per_gp=par["min_genes_per_gp"],
    min_source_genes_per_gp=par["min_source_genes_per_gp"],
    min_target_genes_per_gp=par["min_target_genes_per_gp"],
    max_genes_per_gp=par["max_genes_per_gp"],
    max_source_genes_per_gp=par["max_source_genes_per_gp"],
    max_target_genes_per_gp=par["max_target_genes_per_gp"],
    filter_genes_not_in_masks=par["filter_genes_not_in_masks"],
)

logger.info("Initializing NicheCompass model...")
model = NicheCompass(
    adata,
    counts_key=par["layer"],
    adj_key=par["input_obsm_spatial_connectivities"],
    gp_names_key=par["output_uns_gp_names"],
    active_gp_names_key=par["output_uns_active_gp_names"],
    gp_targets_mask_key=par["output_varm_gp_targets_mask"],
    gp_sources_mask_key=par["output_varm_gp_sources_mask"],
    latent_key=par["output_obsm_embedding"],
    cat_covariates_keys=par["input_obs_covariates"],
    cat_covariates_no_edges=par["covariate_edges"],
    cat_covariates_embeds_keys=par["output_uns_covariate_embeddings"],
    cat_covariates_embeds_injection=par["covariate_embedding_injection_layers"],
    genes_idx_key=par["output_uns_genes_index"],
    target_genes_idx_key=par["output_uns_target_genes_index"],
    source_genes_idx_key=par["output_uns_source_genes_index"],
    recon_adj_key=par["output_obsp_reconstructed_adj_edge_proba"],
    agg_weights_key=par["output_obsp_agg_weights"],
    include_edge_recon_loss=par["include_edge_recon_loss"],
    include_gene_expr_recon_loss=par["include_gene_expr_recon_loss"],
    include_cat_covariates_contrastive_loss=par[
        "include_cat_covariates_contrastive_loss"
    ],
    gene_expr_recon_dist=par["gene_expr_recon_dist"],
    log_variational=par["log_variational"],
    node_label_method=par["node_label_method"],
    active_gp_thresh_ratio=par["active_gp_thresh_ratio"],
    active_gp_type=par["active_gp_type"],
    n_fc_layers_encoder=par["n_fc_layers_encoder"],
    n_layers_encoder=par["n_layers_encoder"],
    n_hidden_encoder=par["n_hidden_encoder"],
    conv_layer_encoder=par["conv_layer_encoder"],
    encoder_n_attention_heads=par["encoder_n_attention_heads"],
    encoder_use_bn=par["encoder_use_bn"],
    dropout_rate_encoder=par["dropout_rate_encoder"],
    dropout_rate_graph_decoder=par["dropout_rate_graph_decoder"],
    n_addon_gp=par["n_addon_gp"],
    cat_covariates_embeds_nums=par["cat_covariates_embeds_nums"],
    seed=par["random_state"],
    use_cuda_if_available=use_gpu,
)

logger.info("Training NicheCompass model...")
model.train(
    n_epochs=par["n_epochs"],
    n_epochs_all_gps=par["n_epochs_all_gps"],
    n_epochs_no_edge_recon=par["n_epochs_no_edge_recon"],
    n_epochs_no_cat_covariates_contrastive=par[
        "n_epochs_no_cat_covariates_contrastive"
    ],
    lr=par["lr"],
    weight_decay=par["weight_decay"],
    lambda_edge_recon=par["lambda_edge_recon"],
    lambda_gene_expr_recon=par["lambda_gene_expr_recon"],
    lambda_cat_covariates_contrastive=par["lambda_cat_covariates_contrastive"],
    contrastive_logits_pos_ratio=par["contrastive_logits_pos_ratio"],
    contrastive_logits_neg_ratio=par["contrastive_logits_neg_ratio"],
    lambda_group_lasso=par["lambda_group_lasso"],
    lambda_l1_masked=par["lambda_l1_masked"],
    l1_targets_categories=par["l1_targets_categories"],
    l1_sources_categories=par["l1_sources_categories"],
    lambda_l1_addon=par["lambda_l1_addon"],
    edge_val_ratio=par["edge_val_ratio"],
    node_val_ratio=par["node_val_ratio"],
    edge_batch_size=par["edge_batch_size"],
    node_batch_size=par["node_batch_size"],
    n_sampled_neighbors=par["n_sampled_neighbors"],
    use_cuda_if_available=use_gpu,
)

## Save model and data
logger.info("Saving NicheCompass model and data...")
mdata = mu.MuData({par["modality"]: adata})
mdata.write_h5mu(par["output"], compression=par["output_compression"])

model.save(par["output_model"], save_adata=False)
