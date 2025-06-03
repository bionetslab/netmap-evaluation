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
        if ((isinstance(names_ad[0], str)) & (isinstance(names_grn[0], str))):
            reverse_mapping[names_grn[col_ind[i]]] = names_ad[row_ind[i]]
        elif (isinstance(names_ad[0], str)):
            reverse_mapping[col_ind[i]] = names_ad[row_ind[i]]
        elif (isinstance(names_grn[0], str)):
            reverse_mapping[names_grn[col_ind[i]]] = row_ind[i]
        else:
            reverse_mapping[col_ind[i]] = row_ind[i]
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
        net_size = net.shape[0]
        overlap_perc = overlap_count/net_size
        net_recovery[net_name] = {'edge_overlap':overlap_count, 'net_size': net_size, 'overlap_percent':overlap_perc}

        net_recovery['total_number_edges'] = grn_adata.var.shape[0]
    return net_recovery


def process(grn_adata):
    if not scs.issparse(grn_adata.X):
        grn_adata.X[np.isnan(grn_adata.X) ] = 0
    grn_adata = downstream_recipe(grn_adata)
    grn_adata = spectral_clustering(grn_adata)
    return grn_adata

def get_top_edges_per_cell(grn_adata, nets, top_edges):
    b = np.argpartition(np.abs(grn_adata.X), top_edges)    # top 3 values from each row
    top_idx = b[:,-top_edges:]

    scgenerai_var_index_np = np.array(grn_adata.var)
    top_edges_per_cell = pd.DataFrame(np.concat(scgenerai_var_index_np[top_idx]))
    top_edges_per_cell.columns = ['source', 'target', 'n_cells']
    top_edges_per_cell['cell_barcode'] = np.repeat(grn_adata.obs.index, repeats=top_edges)
    for net_name, net in nets:
        print(top_edges_per_cell.merge(net, left_on=['source', 'target'], right_on = ['source', 'target']).groupby(['source', 'target']).count().shape[0])
    return top_edges_per_cell



def get_top_edges_per_cell_per_cluster(adata, grn, nets, cluster_var = 'spectral_remap', top_edges=500):
    
    adata = adata[adata.obs[cluster_var] == grn]
    # identify the genes with the highest mean
    top_edges_per_cell = get_top_edges_per_cell(adata, nets, top_edges)
    return top_edges_per_cell
    


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


def compute_metrics(grn_ads, nets):
    collect_results = {}
    for method in grn_ads:
        print(method)
        try:
            collect_results[method] = compute_egde_overlaps_simple(grn_ads[method], nets)
            grn_ads[method] = process(grn_ads[method])
            grn_ads[method], grn_ads[method], score = unify_group_labelling(grn_ads[method], grn_ads[method], 'grn', 'spectral')
            collect_results[method]['clustering_score'] = float(score)
        except FileNotFoundError:
            collect_results[method] = 'no data'
    return collect_results


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



    
    return parser.parse_args()

if  __name__ == '__main__':

    import argparse

    args = parse_args()
    print(f"Netmap Configuration File: {args.netmap_config}")
    print(f"CSnet Configuration File: {args.csnet_config}")
    print(f"ScGeneRAI Configuration File: {args.scgenerai_config}")
    print(f"Dataset Configuration File: {args.dataset_config}")
    print(f"Pipeline Configuration File: {args.pipeline_config}")


    # READ all configurations
    netmap_config = NetmapConfig.read_yaml(args.netmap_config)
    scgenerai_config = ScGeneRAIConfig.read_yaml(args.scgenerai_config)
    csnet_config = CsNetConfig.read_yaml(args.csnet_config)
    dataset_config = DataSimulationConfig.read_yaml(args.dataset_config)
    pipeline_config = PipelineConfig.read_yaml(args.pipeline_config)


    # read network files
    nets = [(op.basename(op.dirname(filename)), pd.read_csv(op.join(pipeline_config.clustered_network_dir, filename), sep='\t')) for filename in dataset_config.edgelist]

    # Augmented net contains edges between genes which are controlled by the same transcription factor
    augmented_nets = [(net[0], build_augmented_network(net[1])) for net in nets]

    grn_ads = {}
    # READ all anndata objects
    try:
        scgenerai = sc.read_h5ad(op.join(scgenerai_config.output_directory, scgenerai_config.adata_filename))
        grn_ads['scgenerai'] = scgenerai
    except FileNotFoundError:
        print('ScGeneRAI not found')
    
    try:
        csnet = sc.read_h5ad(op.join(csnet_config.output_directory,csnet_config.filename+".csn.h5ad"))
        grn_ads['csnet'] = csnet
    except FileNotFoundError:
        print('ScGeneRAI not found')
    
    try:
        netmap = sc.read_h5ad(op.join(netmap_config.output_directory,netmap_config.adata_filename))
        grn_ads['netmap'] = netmap
    except FileNotFoundError:
        print('Ntmap not found')


    collect_results = compute_metrics(grn_ads=grn_ads, nets= nets)

    outdir = op.join(pipeline_config.summary_output_dir, dataset_config.dataset_id)
    os.makedirs(outdir, exist_ok = True)
    write_config(collect_results, file=op.join(outdir, 'results.yaml'))


    
