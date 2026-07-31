# %%
import anndata as ad
import pandas as pd
from itertools import combinations
import numpy as np

ncores = 64

import sys

import anndata
from netmap.downstream.edge_selection import *
from netmap.downstream.clustering import *
from netmap.downstream.plotting import *

from netmap.downstream.regulon import *
import warnings

from netmap.utils.data_utils import *
from netmap.utils.tf_utils import *
from netmap.utils.netmap_config import NetmapConfig

from netmap.model.train_model import create_model_zoo
from netmap.grn.inferrence import inferrence, inferrence_model_wise
from netmap.masking.internal import *
from netmap.masking.external import *


import scipy.sparse as scs


import numpy as np
from pathlib import Path
from matplotlib import pyplot as plt

import argparse
import yaml


import pandas as pd
from itertools import combinations

def get_overlapping_signatures(all_regulons, focus = 'target'):
    genes_by_ct = all_regulons.groupby('cluster')[focus].apply(set).to_dict()
    all_cts = sorted(genes_by_ct.keys())
    pairwise_results = []

    # 2. Iterate through all pairwise combinations
    for ct1, ct2 in combinations(all_cts, 2):
        genes1 = genes_by_ct[ct1]
        genes2 = genes_by_ct[ct2]
        
        # Union of all genes in either cell type for this specific pair
        all_genes_in_pair = genes1.union(genes2)
        
        for gene in all_genes_in_pair:
            # Logic to determine status
            if gene in genes1 and gene in genes2:
                status = 'both'
            elif gene in genes1:
                status = f'only {ct1}'
            else:
                status = f'only {ct2}'
                
            pairwise_results.append({
                'celltype_1': ct1,
                'celltype_2': ct2,
                'status': status,
                'gene': gene
            })

    # 3. Create the final DataFrame
    pairwise_df = pd.DataFrame(pairwise_results)

    return pairwise_df
    

def save_clustering_score_to_yaml(analysis_output_dir, score):
    data = {'clustering_score': float(score) }
    with open(op.join(analysis_output_dir, 'clustering_score.yaml'), 'w') as file:
        yaml.dump(data, file, default_flow_style=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Netmap Analysis Script")
    parser.add_argument("--adata_raw", required=True, help="Path to raw h5ad")
    parser.add_argument("--collectri", required=True, help="Path to collectri.tsv")
    parser.add_argument("--results_dir", required=True, help="Path to the experiment folder containing grn/ and tsvs")
    parser.add_argument("--experiment_name", required=True, help="Experiment name")
    parser.add_argument("--min_score", type=float, default=0.85, help="Score threshold")
    parser.add_argument("--cluster_column", type=str, default='leiden_remap', help="Cluster column to use")
    parser.add_argument("--analysis_dir", required=True, help="Path to the experiment folder containing grn/ and tsvs")
    parser.add_argument("--save_objects", action='store_true', help= 'save the output object of the analysis to file')
    parser.add_argument("--rerun", action='store_true', help= 'save the output object of the analysis to file')


    args = parser.parse_args()

    # Paths setup
    experiment_name = args.experiment_name
    output_dir = op.join(args.results_dir, experiment_name)
    output_dir_grn = Path(op.join(output_dir, "grn"))
    
    analysis_output_dir = op.join(args.analysis_dir, "analysis", experiment_name)
    os.makedirs(analysis_output_dir, exist_ok=True)

    if not args.rerun and op.exists(op.join(analysis_output_dir, 'clustering_score.yaml')): 
        sys.exit("Analysis already exists. Skipping...")

    collectri = pd.read_csv(args.collectri, sep='\t')
    adata_raw = sc.read_h5ad(args.adata_raw)
    var = pd.read_csv(op.join(output_dir, f'{experiment_name}_var.tsv'))
    obs = pd.read_csv(op.join(output_dir, f'{experiment_name}_obs.tsv'))
    grn_adata = ad.AnnData(shape=(obs.shape[0], var.shape[0]), var=var, obs=obs)

    grn_adata.obs['celltype_semi_manual'] = adata_raw.obs['celltype_semi_manual'].values


    grn_adata = retrieve_top_edges(grn_adata, output_dir_grn, percentage=0.1 )
    grn_adata.var = grn_adata.var.set_index('index')


    with plt.rc_context(): 
        sc.pl.umap(adata_raw, color = ['celltype_semi_manual', 'leiden'], show=False)
        plt.savefig(op.join(analysis_output_dir, 'gex_clustering_annotation_umap.pdf'), bbox_inches="tight")



    sc.tl.pca(grn_adata, svd_solver = 'randomized', zero_center = False)
    sc.pp.neighbors(grn_adata, n_neighbors = 50)
    sc.tl.umap(grn_adata, n_components =50)
    sc.tl.tsne(grn_adata)
    sc.tl.leiden(grn_adata, resolution=0.29)


    with plt.rc_context(): 
        sc.pl.tsne(grn_adata, color = ['leiden', 'celltype_semi_manual'], show=False)
        plt.savefig(op.join(analysis_output_dir, 'grn_clustering_comparison_tsne.pdf'), bbox_inches="tight")



    print("\nClustering Score Evaluation")
    score, maping = unify_group_labelling(  adata_raw,grn_adata, 'celltype_semi_manual', 'leiden', True)
    print(score)
    save_clustering_score_to_yaml(analysis_output_dir, score)

    if score < args.min_score:
        raise ValueError(f"Clustering score {score} is below the required threshold of {args.min_score}")

    adata_raw.obs['celltype_remapped'] = grn_adata.obs['leiden_remap'].values


    grn_adata2 = ad.AnnData(shape = (obs.shape[0], var.shape[0]), var = var, obs = obs)
    grn_adata2.var = grn_adata2.var.set_index('index')
    grn_adata2.obs = grn_adata.obs


    add_neighbourhood_expression_mask(adata_raw, grn_adata2, strict=False, layer = 'counts' )
    grn_adata2 = add_cluster_based_candidate_edges(grn_adata2, threshold=0.5)


    index_list = np.where(grn_adata2.var['candidate_edge'])[0]
    grn_adata3 = retrieve_edges_by_index(grn_adata2, output_dir_grn, index_list)
    grn_adata3.obs = grn_adata2.obs
    grn_adata3.obsm['X_pca'] = grn_adata.obsm['X_pca']
    grn_adata3.obsm['X_umap'] = grn_adata.obsm['X_umap']

    cross_tab = pd.crosstab(grn_adata3.obs['celltype_semi_manual'], grn_adata3.obs['leiden_remap'])
    cross_tab.to_csv(op.join(analysis_output_dir, f'crosstab.tsv'), sep = '\t')


    grn_adata3 = add_external_grn(grn_ad=grn_adata3,external_grn=collectri,name_grn='collectri')

    cluster_column = 'leiden_remap'
    min_reg_size = 1
    add_neighbourhood_expression_mask(adata_raw, grn_adata3, strict=False, layer = 'counts' )




    for tp in [10, 20, 30, 40, 50]:
        keep_edges = select_top_edges(gene_inter_adata=grn_adata3, adata=adata_raw, top_per_source=tp, col_cluster=cluster_column, min_reg_size=min_reg_size, verbose=True)
        all_regulons = make_cluster_regulon_dataframe(keep_edges)
        overlapping_signatures = get_overlapping_signatures(all_regulons, focus = 'target')
        overlapping_signatures.to_csv(op.join(analysis_output_dir, f'overlapping_signatures_{tp}_target.tsv'), sep = '\t')

        overlapping_signatures = get_overlapping_signatures(all_regulons, focus = 'source')
        overlapping_signatures.to_csv(op.join(analysis_output_dir, f'overlapping_signatures_{tp}_source.tsv'), sep = '\t')

        for focus in ['source', 'target']:
            regus = aggregate_edges(keep_edges, grn_adata3, key='unique')
            all_signatures = ad.AnnData(X = regus, obs = adata_raw.obs, obsm = adata_raw.obsm)
            all_signatures.obs['leiden_remap'] = grn_adata3.obs['leiden_remap'].values

            sc.tl.rank_genes_groups(all_signatures, cluster_column, method='wilcoxon', key_added = "wilcoxon")
            
            
            regulons_collector = []
            bc_rank_collector = []
            for group in all_signatures.obs[cluster_column].unique():
                print(group)
                bcrank = sc.get.rank_genes_groups_df(all_signatures, group = group, key  = 'wilcoxon')
                print(all_signatures)
                bcrank[['celltype', 'gene']] = bcrank['names'].str.rsplit('_', n=1, expand=True)
                bcrank['group'] = group

                print(all_regulons)
                all_regulonsm  = all_regulons.merge(bcrank, left_on = [focus, 'cluster'], right_on=['gene', 'celltype'])
                all_regulonsm['group'] = group
                regulons_collector.append(all_regulonsm)
                bc_rank_collector.append(bcrank)
            
            regulons_collector = pd.concat(regulons_collector)
            bc_rank_collector = pd.concat(bc_rank_collector)
            
            regulons_collector.to_csv(op.join(analysis_output_dir, f'regulons_{tp}_{focus}.tsv'), sep = '\t')
            bc_rank_collector.to_csv(op.join(analysis_output_dir, f'regulon_rank_{tp}_{focus}.tsv'), sep = '\t')



    


    for tp in [10, 20, 30, 40, 50]:
        keep_edges = select_top_edges(gene_inter_adata=grn_adata3, adata=adata_raw, top_per_source=tp, col_cluster=cluster_column, min_reg_size=min_reg_size, verbose=True, tf_column = 'is_source_collectri')
        all_regulons = make_cluster_regulon_dataframe(keep_edges)
        
        overlapping_signatures = get_overlapping_signatures(all_regulons, focus = 'target')
        overlapping_signatures.to_csv(op.join(analysis_output_dir, f'overlapping_signatures_{tp}_target_TF.tsv'), sep = '\t')

        overlapping_signatures = get_overlapping_signatures(all_regulons, focus = 'source')
        overlapping_signatures.to_csv(op.join(analysis_output_dir, f'overlapping_signatures_{tp}_source_TF.tsv'), sep = '\t')

        for focus in ['source', 'target']:
            regus = aggregate_edges(keep_edges, grn_adata3, key='unique')
            all_signatures = ad.AnnData(X = regus, obs = adata_raw.obs, obsm = adata_raw.obsm)
            all_signatures.obs['leiden_remap'] = grn_adata3.obs['leiden_remap'].values

            sc.tl.rank_genes_groups(all_signatures, cluster_column, method='wilcoxon', key_added = "wilcoxon")
            
            regulons_collector = []
            bc_rank_collector = []
            for group in all_signatures.obs[cluster_column].unique():
                bcrank = sc.get.rank_genes_groups_df(all_signatures, group = group, key  = 'wilcoxon')
                bcrank[['celltype', 'gene']] = bcrank['names'].str.rsplit('_', n=1, expand=True)
                bcrank['group'] = group

                all_regulonsm  = all_regulons.merge(bcrank, left_on = [focus, 'cluster'], right_on=['gene', 'celltype'])
                all_regulonsm['group'] = group
                regulons_collector.append(all_regulonsm)
                bc_rank_collector.append(bcrank)
            
            regulons_collector = pd.concat(regulons_collector)
            bc_rank_collector = pd.concat(bc_rank_collector)
            
            regulons_collector.to_csv(op.join(analysis_output_dir, f'regulons_{tp}_{focus}_TF.tsv'), sep = '\t')
            bc_rank_collector.to_csv(op.join(analysis_output_dir, f'regulon_rank_{tp}_{focus}_TF.tsv'), sep = '\t')


    


    if args.save_objects:
        grn_adata3.write_h5ad(op.join(analysis_output_dir, f'{experiment_name}_processed.h5ad'))
