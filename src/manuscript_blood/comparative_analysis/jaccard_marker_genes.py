

import anndata as ad
import pandas as pd
from itertools import combinations
import numpy as np

# Number of cores to use
ncores = 64

import sys

import anndata
#from netmap.downstream import final_downstream
from netmap.downstream.edge_selection import *
from netmap.downstream.clustering import *
#from netmap.downstream.final_downstream import *

from netmap.downstream.regulon import *
import warnings

from netmap.utils.data_utils import *
from netmap.utils.tf_utils import *
from netmap.utils.netmap_config import NetmapConfig

from netmap.model.train_model import create_model_zoo
from netmap.grn.inferrence import inferrence, inferrence_model_wise
from netmap.masking.internal import *
from netmap.masking.external import *

from netmap.downstream.markers import *
import scipy.sparse as scs


import numpy as np
from pathlib import Path
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Netmap Analysis Script")
    parser.add_argument("--adata_raw", required=True, help="Path to raw h5ad")
    parser.add_argument("--results_dir", required=True, help="Path to the experiment folder containing grn/ and tsvs")
    parser.add_argument("--experiment_name", required=True, help="Experiment name")
    parser.add_argument("--cluster_column", type=str, default='leiden_remap', help="Cluster column to use")
    parser.add_argument("--analysis_dir", required=True, help="Path to the experiment folder containing grn/ and tsvs")


    args = parser.parse_args()

    # Paths setup
    experiment_name = args.experiment_name
    output_dir = op.join(args.results_dir, experiment_name)
    output_dir_grn = Path(op.join(output_dir, "grn"))
    
    analysis_output_dir = op.join(args.analysis_dir, "analysis", experiment_name)
    os.makedirs(analysis_output_dir, exist_ok=True)

    cluster_column = args.cluster_column

    pangalao_markers = pd.read_csv('/data_nfs/og86asub/netmap/netmap-evaluation/data/markers/PanglaoDB_markers_27_Mar_2020.tsv', sep  = '\t')
    cellmarker = pd.read_csv('/data_nfs/og86asub/netmap/netmap-evaluation/data/markers/Cell_marker_Seq.csv', sep ='\t')
    cluster_assignment = pd.read_csv('/data_nfs/og86asub/netmap/netmap-evaluation/data/markers/celltype_clustering.tsv', sep = '\t')


    adata_raw = sc.read_h5ad(args.adata_raw)


    grn_adata3 = sc.read_h5ad(op.join(analysis_output_dir, f'{experiment_name}_processed.h5ad'))
    
    top_per_source = 10
    min_reg_size = 1


    keep_edges = select_top_edges(gene_inter_adata=grn_adata3, adata=adata_raw, top_per_source=top_per_source, col_cluster=cluster_column, min_reg_size=min_reg_size, verbose=True)
    all_regulons = make_cluster_regulon_dataframe(keep_edges)
    regus = aggregate_edges(keep_edges, grn_adata3, key='unique')
    all_signatures = ad.AnnData(X = regus, obs = adata_raw.obs, obsm = adata_raw.obsm)
    all_signatures.obs['leiden_remap'] = grn_adata3.obs['leiden_remap'].values

     
    # Compute differential regulon activity using the nonparametric Wilcoxon test.
    # This identifies regulons whose activity significantly differs between clusters.
    sc.tl.rank_genes_groups(all_signatures, cluster_column, method='wilcoxon', key_added = "wilcoxon")

    # Visualize the top differentially active regulons per cluster.
    # The plot shows effect sizes and significance for the highest-ranked regulons.
    sc.pl.rank_genes_groups(all_signatures, n_genes=25, sharey=False, key="wilcoxon")



    # Handmade cluster mapping
    cluster_mapper = {'b_cells': ['B cells', 'B cells/ Dendritic cells'], 'cd14+_monocytes': ['Monocytes'], 'cd14-_monocytes': ['Monocytes'], 'cd4+_tcells': ['T cells', 'T cells/ NK cells'], 'nk_cells': ['T cells', 'T cells/ NK cells'],
    'DC':['B cells/ Dendritic cells'], 'cd8+_tcells': ['T cells', 'T cells/ NK cells'], 'pDC':['B cells/ Dendritic cells']}

    cm_filtered = cellmarker[
        (cellmarker.tissue_type.isin(['Peripheral blood', 'Blood'])) & 
        (cellmarker.species == 'Human')
    ]
    marker_sets = cm_filtered.groupby('cell_name')['marker'].apply(set).to_dict()

    def check_match(row):
        allowed = cluster_mapper.get(row['ct'], [])
        return 'on target' if row['cluster_high_level'] in allowed else 'off target'



    jaccard_data = prepare_jaccard_analysis_df(grn_adata3, all_signatures, keep_edges, marker_sets,  cluster_mapper)
    jaccard_data = jaccard_data.merge(cluster_assignment, left_on='celltype', right_on='index')
    jaccard_data['is_mapped'] = jaccard_data.apply(check_match, axis=1)
    jaccard_data = jaccard_data[jaccard_data['ct'].isin(cluster_mapper.keys())].copy()
    jaccard_data['database'] = 'cellmarker'


    #plot_jaccard_comparison(jaccard_data)

    ## Do the same for PangalaoDb
    immune =pangalao_markers[pangalao_markers['cell type'].str.contains('B cells|T cells|Monocytes|Dendritic|NK|Natural killer|Macrophages')]
    marker_sets = immune.groupby('cell type')['official gene symbol'].apply(set).to_dict()

    cluster_mapper = {'B cell': ['B cells', 'B cells memory', 'B cells naive'], 'CD14+ Mono': ['Monocytes', 'Macrophages'], 'CD14- Mono': ['Monocytes', 'Macrophages'], 'CD4+ T cell': ['T cells','Gamma delta T cells', 'T cells', 'T cells naive' ], 'NK cell': ['NK cells'],
    'DC':['Dendritic cells'], 'CD8+ T cell':['T cells','Gamma delta T cells', 'T cells', 'T cells naive' ],    'pDC':['Plasmacytoid dendritic cells'],}

    def check_match(row):
        allowed = cluster_mapper.get(row['ct'], [])
        return 'on target' if row['celltype'] in allowed else 'off target'


    full_df = prepare_jaccard_analysis_df(grn_adata3, all_signatures, keep_edges, marker_sets,  cluster_mapper)
    full_df = pd.DataFrame(full_df)
    full_df['is_mapped'] = full_df.apply(check_match, axis=1)
    full_df['database'] = 'pangalaodn'

    #plot_jaccard_comparison(full_df)


    # Concat and write results
    final_re  = pd.concat([jaccard_data, full_df])
    final_re.to_csv(op.join(analysis_output_dir, 'marker_gene_jaccard_overlaps.tsv'), sep = '\t')


