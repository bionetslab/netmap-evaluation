"""Run Netmap on the shared ~100-gene PBMC subset (see prepare_common_pbmc_100genes.py),
then sparsify the resulting GRN the same way as the blood manuscript pipeline
(src/manuscript_blood/netmap_downstream/batch_analyse_output.py).

train_model() mirrors train_models_blood.py's train_model() (same
create_model_zoo/inferrence calls, same NBAutoencoder config).

sparsify_grn() replicates batch_analyse_output.py's core sparsification
sequence (the parts that produce a sparsified GRN object, not its
CollecTRI-annotation/regulon-sweep/UCell reporting extras, which are
manuscript-reporting steps rather than inputs the ATAC mask or benchmark
need):
  1. retrieve_top_edges() to load real attribution values for the top
     `top_edges_percentage` globally-important edges (needed to cluster the
     GRN at all -- the raw inferrence() output is parquet-backed with an
     empty .X).
  2. Cluster the GRN itself (PCA(zero_center=False) -> neighbors -> UMAP ->
     Leiden) and align its clusters to `celltype_col` via
     unify_group_labelling(), gated on `min_score` -- same diagnostic blood
     uses to confirm the GRN's own structure recovers real cell types before
     trusting cluster-wise sparsification.
  3. add_neighbourhood_expression_mask() + add_cluster_based_candidate_edges()
     to flag, per Leiden cluster, edges whose source+target are co-expressed
     in >= `mask_threshold` of that cluster's cells -- the actual
     sparsification step.
  4. retrieve_edges_by_index() to load real attribution values for just the
     surviving ("candidate") edges.

Both steps are timed and written to runtime_all.yaml via write_config, the
same convention run_linger.py uses, so total runtime is directly comparable.
"""
import os
import os.path as op
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy.sparse as scs
import torch

from netmap.model.train_model import create_model_zoo
from netmap.grn.inferrence import inferrence
from netmap.utils.data_utils import retrieve_top_edges, retrieve_edges_by_index
from netmap.downstream.clustering import unify_group_labelling
from netmap.masking.internal import add_neighbourhood_expression_mask, add_cluster_based_candidate_edges

from src.utils import write_config


def train_model(adata, outdir, model_name):
    gene_names = np.array(adata.var.index)
    data_tensor = adata.X  # Log normalized, but not standardized data.

    if scs.issparse(data_tensor):
        data_tensor = torch.tensor(data_tensor.todense(), dtype=torch.float32)
    else:
        data_tensor = torch.tensor(data_tensor, dtype=torch.float32)

    model_zoo = create_model_zoo(data_tensor, n_models=10, n_epochs=10000, model_type='NBAutoencoder',
                                  latent_dim=8, dropout_rate=0.1, hidden_dim=[64])

    grn_adata = inferrence(model_zoo, data_tensor.cuda(), gene_names, xai_method='GradientShap',
                            background_type='zeros', backing_file=op.join(outdir, f'{model_name}.parquet'),
                            return_in_memory=False)

    grn_adata.obs = adata.obs
    grn_adata.write_h5ad(op.join(outdir, f'{model_name}_grn.h5ad'))
    grn_adata.var.to_csv(op.join(outdir, f'{model_name}_var.tsv'), sep='\t')
    adata.obs.to_csv(op.join(outdir, f'{model_name}_obs.tsv'), sep='\t')


def sparsify_grn(adata_raw, outdir, model_name, celltype_col='label', min_score=0.85,
                  top_edges_percentage=0.1, leiden_resolution=0.29, n_neighbors=50, mask_threshold=0.5):
    output_dir_grn = Path(op.join(outdir, 'grn'))
    var = pd.read_csv(op.join(outdir, f'{model_name}_var.tsv'), sep='\t')
    obs = pd.read_csv(op.join(outdir, f'{model_name}_obs.tsv'), sep='\t')

    # --- cluster the GRN itself on its top edges, align clusters to real cell types ---
    grn_adata = ad.AnnData(shape=(obs.shape[0], var.shape[0]), var=var.copy(), obs=obs.copy())
    grn_adata.obs[celltype_col] = adata_raw.obs[celltype_col].values
    grn_adata = retrieve_top_edges(grn_adata, output_dir_grn, percentage=top_edges_percentage)
    grn_adata.var = grn_adata.var.set_index('index')

    sc.tl.pca(grn_adata, svd_solver='randomized', zero_center=False)
    sc.pp.neighbors(grn_adata, n_neighbors=n_neighbors)
    sc.tl.umap(grn_adata)
    sc.tl.leiden(grn_adata, resolution=leiden_resolution)

    score, _mapping = unify_group_labelling(adata_raw, grn_adata, celltype_col, 'leiden', True)
    if score < min_score:
        raise ValueError(
            f"GRN-cluster/{celltype_col} alignment score {score:.3f} is below min_score={min_score}; "
            "the GRN clustering doesn't recover real cell types well enough to trust "
            "cluster-wise sparsification. Lower --min_score to proceed anyway."
        )

    # --- reload all edges (not just the top %), sparsify by per-cluster co-expression support ---
    grn_adata2 = ad.AnnData(shape=(obs.shape[0], var.shape[0]), var=var.copy(), obs=obs.copy())
    grn_adata2.var = grn_adata2.var.set_index('index')
    grn_adata2.obs = grn_adata.obs

    if 'X_pca' not in adata_raw.obsm:
        sc.pp.pca(adata_raw)
    add_neighbourhood_expression_mask(adata_raw, grn_adata2, strict=False, layer='counts')
    grn_adata2 = add_cluster_based_candidate_edges(grn_adata2, cluster_column='leiden_remap',
                                                    threshold=mask_threshold)

    index_list = np.where(grn_adata2.var['candidate_edge'])[0]
    grn_adata3 = retrieve_edges_by_index(grn_adata2, output_dir_grn, index_list)
    grn_adata3.obs = grn_adata2.obs
    grn_adata3.obsm['X_pca'] = grn_adata.obsm['X_pca']
    grn_adata3.obsm['X_umap'] = grn_adata.obsm['X_umap']

    grn_adata3.write_h5ad(op.join(outdir, f'{model_name}_sparsified.h5ad'))
    return grn_adata3


def run(rna_h5ad, outdir, model_name='pbmc-100genes', celltype_col='label', min_score=0.85,
        top_edges_percentage=0.1, leiden_resolution=0.29, n_neighbors=50, mask_threshold=0.5,
        skip_sparsify=False):
    os.makedirs(outdir, exist_ok=True)

    adata = sc.read_h5ad(rna_h5ad)

    start = time.monotonic()
    train_model(adata, outdir, model_name)
    train_time = time.monotonic() - start

    sparsify_time = None
    if not skip_sparsify:
        start_sparsify = time.monotonic()
        sparsify_grn(adata, outdir, model_name, celltype_col=celltype_col, min_score=min_score,
                     top_edges_percentage=top_edges_percentage, leiden_resolution=leiden_resolution,
                     n_neighbors=n_neighbors, mask_threshold=mask_threshold)
        sparsify_time = time.monotonic() - start_sparsify

    write_config(
        {
            'total_time': train_time + (sparsify_time or 0),
            'train_time': train_time,
            'sparsify_time': sparsify_time,
            'n_cells': adata.n_obs,
            'n_genes': adata.n_vars,
        },
        file=op.join(outdir, 'runtime_all.yaml'),
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run Netmap (+ sparsification) on the shared PBMC ~100-gene subset.")
    parser.add_argument("-r", "--rna_h5ad", type=str, required=True,
                         help="adata_netmap_100genes.h5ad from prepare_common_pbmc_100genes.py")
    parser.add_argument("-o", "--outdir", type=str, required=True)
    parser.add_argument("-m", "--model_name", type=str, default="pbmc-100genes")
    parser.add_argument("--celltype_col", type=str, default="label",
                         help="obs column with reference cell-type labels, used to validate GRN clustering")
    parser.add_argument("--min_score", type=float, default=0.85,
                         help="min GRN-cluster/celltype alignment score required before sparsifying")
    parser.add_argument("--top_edges_percentage", type=float, default=0.1)
    parser.add_argument("--leiden_resolution", type=float, default=0.29)
    parser.add_argument("--n_neighbors", type=int, default=50)
    parser.add_argument("--mask_threshold", type=float, default=0.5,
                         help="min fraction of a cluster's cells an edge must be co-expressed in to survive")
    parser.add_argument("--skip_sparsify", action="store_true", help="only run raw inference, skip sparsification")
    args = parser.parse_args()

    run(args.rna_h5ad, args.outdir, model_name=args.model_name, celltype_col=args.celltype_col,
        min_score=args.min_score, top_edges_percentage=args.top_edges_percentage,
        leiden_resolution=args.leiden_resolution, n_neighbors=args.n_neighbors,
        mask_threshold=args.mask_threshold, skip_sparsify=args.skip_sparsify)
