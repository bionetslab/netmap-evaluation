#!/usr/bin/env python
"""
netmap GRN analysis pipeline (CLI).

Combines the exploratory notebook (blood case-study, cell-type-mapping,
regulon/marker diagnostics) with the batch/argparse version of the same
analysis, so the whole thing can be run from the command line for any
experiment / cell-type set.

Example:
    python run_netmap_analysis.py \
        --adata_raw /data/.../bd-rhap-rep2_with_markers.h5ad \
        --collectri /data/.../collectri.tsv \
        --results_dir /data/.../case_studies/blood \
        --analysis_dir /data/.../results/case_studies/blood \
        --experiment_name bd-rhap-rep2-markers \
        --apply_ct_mapping \
        --extra_plots
"""

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
import os
import os.path as op
import sys
import argparse
from itertools import combinations
from pathlib import Path

import yaml
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scipy.sparse as scs
from matplotlib import pyplot as plt

from netmap.downstream.edge_selection import *
from netmap.downstream.clustering import *
from netmap.downstream.regulon import *
from netmap.downstream.celltype_random_walk import *

from netmap.plotting.pienet import *
from netmap.plotting.grn_scatter import *
from netmap.plotting.entropy import *
from netmap.plotting.regulon_dotplot import *
from netmap.plotting.tf_subgraph import *

from netmap.utils.data_utils import *
from netmap.utils.tf_utils import *
from netmap.utils.netmap_config import NetmapConfig

from netmap.model.train_model import create_model_zoo
from netmap.grn.inferrence import inferrence, inferrence_model_wise
from netmap.masking.internal import *
from netmap.masking.external import *


# Raw -> display cell-type label mapping used in the original blood
# case-study notebook. Only applied when --apply_ct_mapping is passed.
DEFAULT_CT_MAPPING = {
    "b_cells": "B cell",
    "cd14+_monocytes": "CD14+ Mono",
    "cd14-_monocytes": "CD14- Mono",
    "cd4+_tcells": "CD4+ T cell",
    "cd8+_tcells": "CD8+ T cell",
    "DC": "DC",
    "pDC": "pDC",
    "nk_cells": "NK cell",
    'megakaryocytes': 'Megakaryocytes'
}

# Marker genes used for the optional random-walk cell-type validation.
DEFAULT_MARKER_DICT = {
    "B cell": ["CD19", "CD79A", "CD21", "MS4A1", "PAX5"],
    "CD4 T cell": ["CD4", 'TCF7', 'LEF1','ITGB1' ],
    "CD8 T cell": ["CD8A", "CD8B", "GZMA", "GZMK"],
    "NK cell": ["NKG7", "NCAM1", "FCGR3A", "GZMB"],
    "CD14 monocyte": ["VCAN", "S100A8", "S100A9"],
    "CD16 monocyte": ["FCGR3A", "CDKN1C"],
    "DC": ["CD1C", "CST3", "FLT3"],
    "pDC": ["PTCRA", "SMIM5", "LAMP5", "JCHAIN"],
}

pal =   {'CD14+ Mono': (0.12156862745098039, 0.47058823529411764, 0.7058823529411765),  
'CD14- Mono': (0.6509803921568628, 0.807843137254902, 0.8901960784313725), 
'CD4+ T cell' : (0.8901960784313725, 0.10196078431372549, 0.10980392156862745),
 'CD8+ T cell': (0.984313725490196, 0.6039215686274509, 0.6),
  'B cell': (0.9921568627450981, 0.7490196078431373, 0.43529411764705883),
   'NK cell': (1.0, 0.4980392156862745, 0.0),
   'pDC': (0.6980392156862745, 0.8745098039215686, 0.5411764705882353),
   'DC':(0.2, 0.6274509803921569, 0.17254901960784313)}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def get_overlapping_signatures(all_regulons, focus="target"):
    """Pairwise comparison of regulon gene sets between every pair of clusters."""
    genes_by_ct = all_regulons.groupby("cluster")[focus].apply(set).to_dict()
    all_cts = sorted(genes_by_ct.keys())
    pairwise_results = []

    for ct1, ct2 in combinations(all_cts, 2):
        genes1, genes2 = genes_by_ct[ct1], genes_by_ct[ct2]
        for gene in genes1.union(genes2):
            if gene in genes1 and gene in genes2:
                status = "both"
            elif gene in genes1:
                status = f"only {ct1}"
            else:
                status = f"only {ct2}"
            pairwise_results.append(
                {"celltype_1": ct1, "celltype_2": ct2, "status": status, "gene": gene}
            )

    return pd.DataFrame(pairwise_results)


def save_clustering_score_to_yaml(analysis_output_dir, score):
    with open(op.join(analysis_output_dir, "clustering_score.yaml"), "w") as fh:
        yaml.dump({"clustering_score": float(score)}, fh, default_flow_style=False)


def compute_regulon_tables(
    keep_edges, all_regulons, grn_adata, adata_raw, cluster_column,
    analysis_output_dir, tp, suffix=""
):
    """Overlap signatures + wilcoxon regulon-rank tables for one top_per_source value."""
    for focus in ["target", "source"]:
        overlap = get_overlapping_signatures(all_regulons, focus=focus)
        overlap.to_csv(
            op.join(analysis_output_dir, f"overlapping_signatures_{tp}_{focus}{suffix}.tsv"),
            sep="\t",
        )

    for focus in ["source", "target"]:
        regus = aggregate_edges(keep_edges, grn_adata, key="unique")
        all_signatures = ad.AnnData(X=regus, obs=adata_raw.obs, obsm=adata_raw.obsm)
        all_signatures.obs["leiden_remap"] = grn_adata.obs["leiden_remap"].values

        sc.tl.rank_genes_groups(all_signatures, cluster_column, method="wilcoxon", key_added="wilcoxon")

        regulons_collector, bc_rank_collector = [], []
        for group in all_signatures.obs[cluster_column].unique():
            bcrank = sc.get.rank_genes_groups_df(all_signatures, group=group, key="wilcoxon")
            bcrank[["celltype", "gene"]] = bcrank["names"].str.rsplit("_", n=1, expand=True)
            bcrank["group"] = group

            merged = all_regulons.merge(
                bcrank, left_on=[focus, "cluster"], right_on=["gene", "celltype"]
            )
            merged["group"] = group
            regulons_collector.append(merged)
            bc_rank_collector.append(bcrank)

        regulons_collector = pd.concat(regulons_collector)
        bc_rank_collector = pd.concat(bc_rank_collector)

        regulons_collector.to_csv(
            op.join(analysis_output_dir, f"regulons_{tp}_{focus}{suffix}.tsv"), sep="\t"
        )
        bc_rank_collector.to_csv(
            op.join(analysis_output_dir, f"regulon_rank_{tp}_{focus}{suffix}.tsv"), sep="\t"
        )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Netmap GRN analysis pipeline")
    parser.add_argument("--adata_raw", required=True, help="Path to raw h5ad")
    parser.add_argument("--collectri", required=True, help="Path to collectri.tsv")
    parser.add_argument("--results_dir", required=True,
                         help="Folder containing <experiment_name>/ (var/obs tsvs + grn/)")
    parser.add_argument("--analysis_dir", required=True,
                         help="Folder under which analysis/<experiment_name> outputs are written")
    parser.add_argument("--experiment_name", required=True)
    parser.add_argument("--min_score", type=float, default=0.85,
                         help="Minimum acceptable clustering concordance score")
    parser.add_argument("--cluster_column", type=str, default="leiden_remap",
                         help="Cluster column used for regulon selection / DE testing")
    parser.add_argument("--top_per_source", type=int, nargs="+", default=[10, 20, 30, 40, 50],
                         help="List of top_per_source values to sweep")
    parser.add_argument("--min_reg_size", type=int, default=1)
    parser.add_argument("--top_edge_percentage", type=float, default=0.1,
                         help="Fraction of top global edges to load initially")
    parser.add_argument("--n_neighbors", type=int, default=50)
    parser.add_argument("--leiden_resolution", type=float, default=0.29)
    parser.add_argument("--apply_ct_mapping", action="store_true",
                         help="Remap raw celltype_semi_manual labels using DEFAULT_CT_MAPPING")
    parser.add_argument("--extra_plots", action="store_true",
                         help="Run additional diagnostics: marker-based random walk, "
                              "regulon dotplot, celltype network plot")
    parser.add_argument("--save_objects", action="store_true",
                         help="Write processed AnnData object(s) to disk")
    parser.add_argument("--rerun", action="store_true",
                         help="Re-run even if clustering_score.yaml already exists")
    args = parser.parse_args()

    # ---- paths -----------------------------------------------------------
    experiment_name = args.experiment_name
    output_dir = op.join(args.results_dir, experiment_name)
    output_dir_grn = Path(op.join(output_dir, "grn"))
    analysis_output_dir = op.join(args.analysis_dir, "analysis", experiment_name)
    os.makedirs(analysis_output_dir, exist_ok=True)

    if not args.rerun and op.exists(op.join(analysis_output_dir, "clustering_score.yaml")):
        sys.exit("Analysis already exists. Skipping...")

    # ---- load inputs -------------------------------------------------------
    collectri = pd.read_csv(args.collectri, sep="\t")
    adata_raw = sc.read_h5ad(args.adata_raw)

    if args.apply_ct_mapping:
        print(adata_raw.obs["celltype_semi_manual"])
        adata_raw.obs["celltype_semi_manual"] = (
            adata_raw.obs["celltype_semi_manual"].map(DEFAULT_CT_MAPPING)
        )

    var = pd.read_csv(op.join(output_dir, f"{experiment_name}_var.tsv"))
    obs = pd.read_csv(op.join(output_dir, f"{experiment_name}_obs.tsv"))

    #Start after correctly labelling the object
    obs = pd.read_csv(op.join(analysis_output_dir, 'leiden_remap.tsv'), sep = '\t')
    X_tnse=pd.read_csv(op.join(analysis_output_dir, 'tsne_embedding.tsv'), sep = '\t')
    X_pca=pd.read_csv(op.join(analysis_output_dir, 'pca_embedding.tsv'), sep = '\t')


    # ---- candidate edges based on cluster expression ----------------------
    grn_adata2 = ad.AnnData(shape=(obs.shape[0], var.shape[0]), var=var, obs=obs)
    grn_adata2.var = grn_adata2.var.set_index("index")
    grn_adata2.obs = grn_adata.obs

    add_neighbourhood_expression_mask(adata_raw, grn_adata2, strict=False, layer = 'counts', mask_data=False )
    grn_adata2 = add_cluster_based_candidate_edges(grn_adata2, threshold=0.5)

    index_list = np.where(grn_adata2.var["candidate_edge"])[0]
    grn_adata3 = retrieve_edges_by_index(grn_adata2, output_dir_grn, index_list)
    grn_adata3.obs = grn_adata2.obs

    grn_adata3.obsm["X_pca"] = X_pca.values
    grn_adata3.obsm["X_tsne"] = X_tnse.values

    cross_tab = pd.crosstab(grn_adata3.obs["celltype_semi_manual"], grn_adata3.obs["leiden_remap"])
    cross_tab.to_csv(op.join(analysis_output_dir, "crosstab.tsv"), sep="\t")

    # ---- external GRN (collectri) -----------------------------------------
    grn_adata3 = add_external_grn(grn_ad=grn_adata3, external_grn=collectri, name_grn="collectri")

        # ---- optionally persist processed AnnData objects ----------------------
    if args.save_objects:
        grn_adata3.write_h5ad(op.join(analysis_output_dir, f"{experiment_name}_processed.h5ad"))


    if args.apply_ct_mapping:
        for old, new in DEFAULT_CT_MAPPING.items():
            grn_adata3.var.columns = grn_adata3.var.columns.str.replace(old, new, regex=False)

    cluster_column = args.cluster_column
    add_neighbourhood_expression_mask(adata_raw, grn_adata3, strict=False, layer="counts")
    grn_adata3.layers["masked"] = np.multiply(grn_adata3.X, grn_adata3.layers["mask"])

    # ---- regulon sweeps: all sources, then TF-restricted -------------------
    for tp in args.top_per_source:
        keep_edges = select_top_edges(
            gene_inter_adata=grn_adata3, adata=adata_raw, top_per_source=tp,
            col_cluster=cluster_column, min_reg_size=args.min_reg_size, verbose=True,
        )
        all_regulons = make_cluster_regulon_dataframe(keep_edges)
        compute_regulon_tables(
            keep_edges, all_regulons, grn_adata3, adata_raw, cluster_column,
            analysis_output_dir, tp,
        )

    for tp in args.top_per_source:
        keep_edges = select_top_edges(
            gene_inter_adata=grn_adata3, adata=adata_raw, top_per_source=tp,
            col_cluster=cluster_column, min_reg_size=args.min_reg_size, verbose=True,
            tf_column="is_source_collectri",
        )
        all_regulons = make_cluster_regulon_dataframe(keep_edges)
        compute_regulon_tables(
            keep_edges, all_regulons, grn_adata3, adata_raw, cluster_column,
            analysis_output_dir, tp, suffix="_TF",
        )

    # ---- optional diagnostic plots (from the exploratory notebook) --------
    if args.extra_plots:
        # Regulon dotplot on the largest (last) sweep's regulons
        keep_edges = select_top_edges(
            gene_inter_adata=grn_adata3, adata=adata_raw, top_per_source=max(args.top_per_source),
            col_cluster=cluster_column, min_reg_size=args.min_reg_size, verbose=False,
        )
        all_regulons = make_cluster_regulon_dataframe(keep_edges)
        regus = aggregate_edges(keep_edges, grn_adata3, key="unique")
        all_signatures = ad.AnnData(X=regus, obs=adata_raw.obs, obsm=adata_raw.obsm)
        all_signatures.obs[cluster_column] = grn_adata3.obs[cluster_column].values
        sc.tl.rank_genes_groups(all_signatures, cluster_column, method="wilcoxon", key_added="wilcoxon")

        pp = rank_regulon_groups_dotplot(
            grn_adata3, adata_regl=all_signatures, return_fig=True,
            original_cluster_column=cluster_column, new_cluster_column=cluster_column,
            n_genes=10, figsize=(20, 2),
        )
        pp.savefig(op.join(analysis_output_dir, "regulon_dotplot.pdf"), bbox_inches="tight")
        plt.close("all")

        # Marker-based cell-type validation via random walk
        proba, celltypes, tgpc = run_celltype_random_walk(
            grn_adata3, DEFAULT_MARKER_DICT, restart_prob=0.2, n_jobs=20, max_iter=100,
            edge_idf=False, scale_by_marker_count=True, idf_log_scale=True, top_genes_k=20,
        )
        grn_adata3.obsm["celltype_rw_proba"] = proba
        grn_adata3.uns["celltype_rw_labels"] = celltypes
        grn_adata3.obs["celltype_rw_pred"] = np.array(celltypes)[proba.argmax(axis=1)]

        rw_cross_tab = pd.crosstab(
            grn_adata3.obs["celltype_semi_manual"], grn_adata3.obs["celltype_rw_pred"]
        )
        rw_cross_tab.to_csv(op.join(analysis_output_dir, "celltype_rw_crosstab.tsv"), sep="\t")

        with plt.rc_context():
            sc.pl.tsne(
                grn_adata3,
                color=["celltype_semi_manual", "celltype_rw_pred"],
                show=False,
            )
            plt.savefig(op.join(analysis_output_dir, "celltype_rw_tsne.pdf"), bbox_inches="tight")
            plt.close()

        # Combined cell-type co-expression network
        network = compute_combined_celltype_network(
            grn_adata3, tgpc, celltype_col=cluster_column, top_n=30,
            max_nodes_per_celltype=20, top_edges_per_celltype=100,
        )
        fig, _ = plot_combined_celltype_network(
            network, k=0.5, seed=1, node_radius_frac=(0.005, 0.025),
            palette=pal, figsize=(15, 15), layout="kk", drop_isolates=True,
        )


        importance = pd.DataFrame({"importance": network['importance']})
        fractions = pd.DataFrame(network['gene_fractions']).T
        fractions = pd.concat([importance, fractions], axis = 1)
        fractions.sort_values('importance', ascending=False)
        fractions.to_csv(op.join(analysis_output_dir, f"{experiment_name}_fractions.tsv"), sep = '\t')

        fig.savefig(op.join(analysis_output_dir, "celltype_network.pdf"), bbox_inches="tight")
        plt.close("all")

    # ---- optionally persist processed AnnData objects ----------------------
    if args.save_objects:
        grn_adata3.write_h5ad(op.join(analysis_output_dir, f"{experiment_name}_processed.h5ad"))


    print(f"\nDone. Outputs written to: {analysis_output_dir}")


if __name__ == "__main__":
    main()