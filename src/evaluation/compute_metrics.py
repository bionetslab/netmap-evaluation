import sys
import os.path as op
import os
sys.path.append('/data_nfs/og86asub/netmap/netmap-evaluation/')

from netmap.src.utils.netmap_config import NetmapConfig
from src.methods.csnet.csnet_config import CsNetConfig
from src.methods.scgenerai.scgenerai_config import ScGeneRAIConfig
from src.data_simulation.data_simulation_config import DataSimulationConfig
from src.pipelines.utils import PipelineConfig
from captum.attr import GradientShap


import anndata
import numpy as np
import scanpy as sc
import pandas as pd
from sklearn.metrics.cluster import contingency_matrix
from scipy.optimize import linear_sum_assignment
import scipy.sparse as scs
from sklearn.cluster import SpectralClustering

from src.utils import write_config

import numpy as np
import pandas as pd
from typing import Optional


def downstream_recipe(adata)-> anndata.AnnData:
    """
    Downstream reciepe for an LRP anndata object:
    TODO: replace the config dict, to pass values
    """
    config = {'min_cells':1, 'n_neighbors': 30, 'leiden_resolution': 0.1, 'n_components': 30, 'knn_neighbors': 30}
    sc.pp.filter_genes(adata, min_cells=config['min_cells'])
    #sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
    #sc.pp.normalize_total(adata)
    #sc.pp.log1p(adata)

    sc.tl.pca(adata, svd_solver = 'randomized', zero_center = False)

    sc.pp.neighbors(adata, n_neighbors=config['knn_neighbors'])
    sc.tl.leiden(adata, resolution=config['leiden_resolution'])

    sc.tl.umap(adata, n_components = config['n_components'])
    return adata


def spectral_clustering(adata, n_clu = 2, key_added = 'spectral'):
    """
    Run sklearn spectral clustering on the neighbour matrix in the anndata object.

    Args:
    adata: Anndata object
    n_clu: Number of clusters to compute
    key_added: The key to add the new labelling to [Default: spectral]
    """
    sc.pp.neighbors(adata)
    ssc = SpectralClustering(n_clusters=n_clu,assign_labels='discretize',random_state=0, affinity= 'precomputed_nearest_neighbors').fit(adata.obsp['distances'])
    counter = 0
    key_added_t = key_added
    while key_added_t in adata.obs.columns:
        counter = counter + 1
        key_added_t = f'{key_added}_{counter}'
    adata.obs[key_added_t] = ssc.labels_
    adata.obs[key_added_t] = pd.Categorical(adata.obs[key_added_t])
    return adata


def unify_group_labelling(adata, grn_adata, col_adata, col_grn_adata, return_mapping=False):
    """
    Adjust group labelling such that grn_adata has the same group label than the 
    corresponding column in adata based on the grn column.

    Returns the data objects and a score of the matching as the cost of the matching
    divided by the number of cells.
    
    """

    cm = contingency_matrix(adata.obs[col_adata], grn_adata.obs[col_grn_adata])
    row_ind, col_ind = linear_sum_assignment(cm, maximize = True)
    
    names_ad = np.unique(adata.obs[col_adata])
    names_grn = np.unique(grn_adata.obs[col_grn_adata])
    mapping = {}
    reverse_mapping = {}
    for i in range(len(row_ind)):
        reverse_mapping[names_grn[col_ind[i]]] = names_ad[row_ind[i]]

    col_grn_adata_remapped = col_grn_adata + '_remap'
    if isinstance(np.unique(grn_adata.obs[col_grn_adata])[0], str):
        grn_adata.obs[col_grn_adata_remapped] = [reverse_mapping[a] for a in grn_adata.obs[col_grn_adata]]
    else:
        grn_adata.obs[col_grn_adata_remapped] = [reverse_mapping[int(a)] for a in grn_adata.obs[col_grn_adata]]

    grn_adata.obs[col_grn_adata_remapped] = pd.Categorical(grn_adata.obs[col_grn_adata_remapped])

    score = cm[row_ind, col_ind].sum()/adata.obs.shape[0]
    if return_mapping:
        return adata, grn_adata, score, reverse_mapping
    else:
        return adata, grn_adata, score

def compute_egde_overlaps_simple(grn_adata, net_list):
    net_recovery = {}
    for net_name, net in net_list:
        overlap_count = (pd.merge(grn_adata.var, net, on=['source', 'target'], how='inner').shape[0])
        reverse_overlap = (pd.merge(grn_adata.var, net, on=['target', 'source'], how='inner').shape[0])
        net_size = net.shape[0]
        overlap_perc = overlap_count/net_size
        reverse_overlap_perc = reverse_overlap/net_size
        net_recovery[net_name] = {'edge_overlap':overlap_count, 'net_size': net_size, 'overlap_percent':overlap_perc, 'reverse_overlap_percent': reverse_overlap_perc}

        net_recovery['total_number_edges'] = grn_adata.var.shape[0]
    return net_recovery



def process(grn_adata):
    if not scs.issparse(grn_adata.X):
        grn_adata.X[np.isnan(grn_adata.X) ] = 0
    grn_adata = downstream_recipe(grn_adata)
    grn_adata = spectral_clustering(grn_adata)
    return grn_adata





def build_augmented_network(net):
    """
    Insert an edge between two genes if they are regulated by the same gene.
    """
    augmented_net = []
    for s in net.source.unique():
        targets = net[net.source == s].target.unique()
        cro = np.transpose([np.tile(targets, len(targets)), np.repeat(targets, len(targets))])
        augmented_net.append(cro)
    augmented_net = pd.DataFrame(np.vstack(augmented_net))
    augmented_net.columns = ['source', 'target']
    augmented_net['edge' ] = augmented_net['source']+'_'+augmented_net['target']
    return augmented_net



def get_top_edges_per_cell(grn_adata, top_edges):
    partition_index = grn_adata.shape[1] - top_edges
    b = np.argpartition(grn_adata.X, partition_index, axis=1)    
    top_idx = b[:, partition_index:]
    
    edge_metadata_np = np.array(grn_adata.var.loc[:, ['source', 'target']])
    top_edges_metadata = edge_metadata_np[top_idx.ravel()]

    top_edges_per_cell = pd.DataFrame(top_edges_metadata, columns=['source', 'target'])
    top_edges_per_cell['cell_barcode'] = np.repeat(grn_adata.obs.index, repeats=top_edges)
    top_edges_per_cell['top_edges'] = top_edges
    
    return top_edges_per_cell

def get_top_edges_global(grn_adata, top_edges: int):
    top = []

    for t in top_edges:
        top.append(get_top_edges_per_cell(grn_adata, t))

    top = pd.concat(top)
    return top


def get_top_edges_by_target(grn_adata, top_edges: int):

    top = []
    for tar in np.unique(grn_adata.var.target):
        
        sub = grn_adata[:, grn_adata.var.target == tar]

        for t in top_edges:
            if t>sub.shape[1]:
                continue
            top.append(get_top_edges_per_cell(sub, t))

    top = pd.concat(top)
    return top


def get_top_edges_per_cell_per_cluster(grn_adata, nets, cluster_var = 'spectral_remap', top_edges=[100], group_by_target = False):
    
    grns = grn_adata.obs[cluster_var].unique()

    top_edges_per_cell_collector = []
    
    for i in range(len(nets)):
        for j in range(len(grns)):
            clu = grns[j]
            grn = nets[i]

            print(f"Cluster {clu} -- GRN {grn[0]}")

            grn_adata_sub = grn_adata[grn_adata.obs[cluster_var] == clu]
            # identify the genes with the highest mean
            if group_by_target:
                top_edges_per_cell = get_top_edges_by_target(grn_adata_sub,  top_edges)
            else:
                top_edges_per_cell = get_top_edges_global(grn_adata_sub,  top_edges)


            for te in top_edges:
                agg = top_edges_per_cell[top_edges_per_cell.top_edges==te].merge(grn[1], left_on=['source', 'target'], right_on = ['source', 'target']).loc[:,['cell_barcode']].groupby(['cell_barcode']).value_counts()
                
                agg = agg.reset_index()
                agg.columns = ['cell_barcode', 'egde_count']
                mean_edges_recovered = agg.loc[ :, "egde_count"].median()
                number_of_cells = agg.shape[0]

                #print(agg)
                tt = top_edges_per_cell[top_edges_per_cell.top_edges==te].loc[:,['source', 'target']].groupby(['source', 'target']).value_counts()

                #print(tt)
                agg['n_top'] = te
                agg['net'] = grn[0]
                # n_top, grn, TP = overlap, FP = top_k-TP
                # Number edges 
                n_edges = top_edges_per_cell[top_edges_per_cell.top_edges==te].shape[0]
                TP = agg.shape[0]
                FP = n_edges - TP
                FN = grn[1].shape[0]-TP
                gold_standard_edges = grn[1].shape[0]
                current_collector =  [te, grn[0], mean_edges_recovered, number_of_cells,  gold_standard_edges, grn_adata.var.shape[0]]

                if i == j:
                    current_collector.append('on_target')
                else: 
                    current_collector.append('off_target')
                
                #print(current_collector)
                top_edges_per_cell_collector.append(current_collector)


    top_edges_per_cell_collector = pd.DataFrame(top_edges_per_cell_collector)
    top_edges_per_cell_collector.columns = ['n_top', 'net', 'avg_edges_recovered', 'n_cells', 'gold_standard_edges', 'max_possible_edges', 'target']

    return top_edges_per_cell_collector


def reformat_dataframe(recovery_rates, config_name):
        # Assuming your DataFrame is named df
    recovery_rates = recovery_rates.pivot_table(
        index=['n_top', 'net'],
        columns='type',
        values='percentage_recovered'
    ).reset_index()
    recovery_rates['method'] = config_name
    recovery_rates = recovery_rates.loc[:, ['method', 'n_top', 'on_target', 'off_target']]
    return recovery_rates




def compute_metrics(grn_ads, nets, augmented_nets, global_nets, group_key='grn', group_by_target = False):
    collect_results = {}
    results = []
    for method in grn_ads:
        print(f'Running for {method}')
        
        collect_results[method] = compute_egde_overlaps_simple(grn_ads[method], nets)
        grn_ads[method] = process(grn_ads[method])
        print(f'Unify group for {method}')

        grn_ads[method], grn_ads[method], score = unify_group_labelling(grn_ads[method], grn_ads[method], group_key, 'spectral')
        collect_results[method]['clustering_score'] = float(score)

        if group_by_target:
            top_edges = [1,2,3,4,5,10, 15, 20, 25, 50]
        else:
            top_edges = [10, 50, 100, 500, 1000, 2500, 5000, 10000, 25000]
 

        print(f'Getting top edges for {method}')

        start = time.monotonic()
        top_edges_per_cell_collector = get_top_edges_per_cell_per_cluster(grn_ads[method], nets, cluster_var = 'spectral_remap', top_edges=top_edges, group_by_target = group_by_target)
        top_edges_per_cell_collector_augmented = get_top_edges_per_cell_per_cluster(grn_ads[method], augmented_nets, cluster_var = 'spectral_remap', top_edges=top_edges, group_by_target = group_by_target)
        top_edges_per_cell_collector_global = get_top_edges_per_cell_per_cluster(grn_ads[method], global_nets, cluster_var = 'spectral_remap', top_edges=top_edges,group_by_target = group_by_target)

        
        top_edges_per_cell_collector['method'] = method
        top_edges_per_cell_collector['net_type'] = 'strict'
        top_edges_per_cell_collector_augmented['method'] = method
        top_edges_per_cell_collector_augmented['net_type'] = 'extended'
        top_edges_per_cell_collector_global['method'] = method
        top_edges_per_cell_collector_global['net_type'] = 'unspecific'


        results.append(top_edges_per_cell_collector)
        results.append(top_edges_per_cell_collector_augmented)
        results.append(top_edges_per_cell_collector_global)


    results = pd.concat(results)
    return results

def parse_args():
    parser = argparse.ArgumentParser(description="Argument parser with two configuration files.")
    
    parser.add_argument(
        '--netmap_config',
        type=str,
        default='configuration.json',
        help='Path to the main configuration file (default: configuration.json)'
    )
    
    parser.add_argument(
        '--csnet_config',
        type=str,
        default='configuration.json',
        help='Path to the main configuration file (default: configuration.json)'
    )

    parser.add_argument(
        '--scgenerai_config',
        type=str,
        default='configuration.json',
        help='Path to the main configuration file (default: configuration.json)'
    )

    parser.add_argument(
        '--dataset_config',
        type=str,
        default='dataset_configuration.json',
        help='Path to the dataset configuration file (default: dataset_configuration.json)'
    )

    parser.add_argument(
        '--pipeline_config',
        type=str,
        default='dataset_configuration.json',
        help='Path to the dataset configuration file (default: dataset_configuration.json)'
    )

    parser.add_argument(
    "--config_list",  # name on the CLI - drop the `--` for positional/required parameters
    nargs="*",  # 0 or more values expected => creates a list
    type=str,
    default=None,  # default if nothing is provided
    )

    args = parser.parse_args()
    config_dict = {}
    for aa in args.config_list:
        n = aa.split('=')
        config_dict[n[0].strip()] = n[1].strip() 
    return args, config_dict


def config_reader(configs):
    config_dict = {}
    for c in configs:
        if  c.startswith('scgenerai'):
            scgenerai_config = ScGeneRAIConfig.read_yaml(configs[c])
            config_dict[c] = scgenerai_config
        elif c.startswith('csnet'):
            csnet_config = CsNetConfig.read_yaml(configs[c])
            config_dict[c] = csnet_config
        else:
            netmap_config = NetmapConfig.read_yaml(configs[c])
            config_dict[c] = netmap_config
    return config_dict

if  __name__ == '__main__':

    import argparse
    import time
    import cProfile

    args, config_dict = parse_args()
    # print(f"Netmap Configuration File: {args.netmap_config}")
    # print(f"CSnet Configuration File: {args.csnet_config}")
    # print(f"ScGeneRAI Configuration File: {args.scgenerai_config}")
    print(f"Dataset Configuration File: {args.dataset_config}")
    print(f"Pipeline Configuration File: {args.pipeline_config}")


    config_dict = config_reader(config_dict)
    dataset_config = DataSimulationConfig.read_yaml(args.dataset_config)
    pipeline_config = PipelineConfig.read_yaml(args.pipeline_config)


    # read network files
    
    nets = [(op.basename(op.dirname(filename)), pd.read_csv(op.join(pipeline_config.clustered_network_dir, filename), sep=dataset_config.separator)) for filename in dataset_config.edgelist]
    off_net = [(op.basename(op.dirname(filename)), pd.read_csv(op.join(pipeline_config.clustered_network_dir, filename), sep=dataset_config.separator)) for filename in dataset_config.common_edges]

    print(nets)
    # Augmented net contains edges between genes which are controlled by the same transcription factor
    augmented_nets = [(net[0], build_augmented_network(net[1])) for net in nets]

    grn_ads = {}

    for c in config_dict:
        try:
            print('reading file')
            if c.startswith('csnet'):
                scgenerai = sc.read_h5ad(op.join(config_dict[c].output_directory, config_dict[c].filename))
            else:
                scgenerai = sc.read_h5ad(op.join(config_dict[c].output_directory, config_dict[c].adata_filename))
                print(scgenerai.var)
            grn_ads[c] = scgenerai
        except:
            continue

    
    # READ all anndata objects
    # try:k_t
    #     print(op.join(scgenerai_config.output_directory, scgenerai_config.adata_filename))
    #     scgenerai = sc.read_h5ad(op.join(scgenerai_config.output_directory, scgenerai_config.adata_filename))
    #     grn_ads['scgenerai'] = scgenerai
    #     print(scgenerai)
    # except FileNotFoundError:
    #     print('ScGeneRAI not found')
    
    # try:
    #     csnet = sc.read_h5ad(op.join(csnet_config.output_directory,csnet_config.filename+".csn.h5ad"))
    #     grn_ads['csnet'] = csnet
    # except FileNotFoundError:
    #     print('ScGeneRAI not found')
    
    # try:
    #     netmap = sc.read_h5ad(op.join(netmap_config.output_directory,netmap_config.adata_filename))
    #     grn_ads['netmap'] = netmap
    #     grn_ads['netmap'].obs['grn'] = pd.Categorical(grn_ads['csnet'].obs['grn'])

    # except FileNotFoundError:
    #     print('Ntmap not found')




    overlaps = compute_metrics(grn_ads=grn_ads, nets= nets, augmented_nets=augmented_nets, global_nets=off_net, group_key=dataset_config.group_key, group_by_target=True)

    overlaps_ungrouped = compute_metrics(grn_ads=grn_ads, nets= nets, augmented_nets=augmented_nets, global_nets=off_net, group_key=dataset_config.group_key, group_by_target=False)

    outdir = op.join(pipeline_config.summary_output_dir, dataset_config.dataset_id)
    os.makedirs(outdir, exist_ok = True)
    overlaps.to_csv(op.join(outdir, 'overlaps_per_target.tsv'), sep='\t')
    overlaps_ungrouped.to_csv(op.join(outdir, 'overlaps_global_top_k.tsv'), sep='\t')



    
