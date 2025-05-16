import os
import sys
os.environ["CUDA_VISIBLE_DEVICES"]="0"


import torch
torch.cuda.is_available()

## External imports
import os
import pandas as pd
import numpy as np
import scanpy as sc
import time
import torch

import warnings
warnings.filterwarnings("ignore")
import os.path as op


from sklearn.model_selection import train_test_split

import pandas as pd

sys.path.append('/data_nfs/og86asub/netmap/netmap-evaluation/')
from netmap.src.utils.data_utils import *
from netmap.src.utils.tf_utils import *

### Internal imports
sys.path.append('/data_nfs/og86asub/netmap/NetMap_LRP')



from src.methods.scgenerai.scgenerai_config import ScGeneRAIConfig
from src.data_simulation.data_simulation_config import DataSimulationConfig

#Import cloned scgenerai/ can't be installed via pip, see issue on github
sys.path.append('/data_nfs/og86asub/netmap/scGeneRAI')
from scGeneRAI import scGeneRAI

def run_scgenerai(config, dataset_config):

    #setup temp dir for scGeneRAI to save results
    temp_dir = config.temp_dir
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    os.environ['scGeneRAI_TEMPDIR'] = temp_dir

    ## Load config and setup outputs
    os.makedirs(config.output_directory, exist_ok=True)
    sc.settings.figdir = config.output_directory

    rerun = config.rerun
    split = config.split

    
    ## load data
    adata = sc.read_h5ad(config.input_data)
    sc.pp.scale(adata)
    
    tf_genes = adata.var.index
    tf_indices, tf_gene_names = filter_tf_names(tf_genes, adata, config.tf_only)
    nr_tfs = len(tf_indices)


    ## Get the data matrix from the CustumAnndata obeject

    data = adata.X

    if split:
        data_train, data_test = train_test_split(data,test_size=config.test_size)
        row_names = adata.obs_names
        column_names = adata.var_names
        #scgenerai needs pandas df
        data_train_df = pd.DataFrame(data_train, index=row_names[:len(data_train)], columns=column_names)
        data_test_df = pd.DataFrame(data_test, index=row_names[len(data_train):], columns=column_names)
    else:
        #train == test when no split
        data_train_df = pd.DataFrame(data, index=adata.obs_names, columns=adata.var_names)
        data_test_df = data_train_df

    if rerun:
        start = time.monotonic()
        model = scGeneRAI()

        if torch.cuda.is_available():
            device = torch.device('cuda')
        else:
            device = torch.device('cpu')
    
        model.fit(data_train_df, model_depth=2, nepochs=100, lr=2e-2, batch_size=5, lr_decay = 0.99, descriptors = None, early_stopping = True, device_name = device)
        
        model.predict_networks(data_test_df, descriptors = None, LRPau = True, remove_descriptors = True, device_name = device, PATH = temp_dir)

        print(f'Elapsed time: {time.monotonic()-start}')

    res_files = os.listdir(op.join(temp_dir))
    #res_data = pd.concat([pd.read_csv(temp_dir + '/results/' + res_file) for res_file in res_files])
    res_data = pd.concat([pd.read_csv(temp_dir + '/results/' + res_file).assign(cell_name=res_file.split('_')[-1].replace('.csv', '')) for res_file in res_files])
    res_data = res_data.drop(res_data.columns[0], axis=1)
    #mean_resda = res_data[['LRP', 'source_gene', 'target_gene']].groupby(['source_gene', 'target_gene']).mean().reset_index()
    #get aa and nl like netmap
    unique_edges = res_data[['source_gene', 'target_gene']].drop_duplicates()

    nl = [(row['source_gene'], row['target_gene']) for _, row in unique_edges.iterrows()]

    pivot_table = res_data.pivot_table(index='cell_name', columns=['source_gene', 'target_gene'], values='LRP', fill_value=0)
    aa = pivot_table.values


    varnames =  [f'{str(x[0])}_{str(x[1])}' for x in np.array(nl)]
    grn_adata = attribution_to_anndata(aa, obs = adata.obs)
    grn_adata.var.index = np.array(varnames)
    grn_adata.write_h5ad(op.join(config.output_directory,config.adata_filename))

    # ## Run one round of 
    # grn_adata = d.downstream_recipe(grn_adata)
    # grn_adata = d.spectral_clustering(grn_adata, n_clu=dataset_config.n_celltypes)

    # ## Compute background for edges
    # rel = post.random_score_distributions(grn_adata)
    # list_of_edges = []
    # for i in rel.keys():
    #     list_of_edges = list_of_edges + rel[i].tolist()

    # #Subset GRN object to relevant edges
    # grn_adata_sub = grn_adata[:, list_of_edges]
    # sc.pp.filter_genes(grn_adata_sub, min_cells=50)


    # ## Run one round of 
    # adata = d.downstream_recipe(adata)
    # adata = d.spectral_clustering(adata, n_clu=dataset_config.n_celltypes)

    # # Compute the cluster matching and relable the spectral columns
    # adata, grn_adata_sub, score = post.unify_group_labelling(adata, grn_adata_sub, col_adata = 'spectral', col_grn_adata = 'spectral')

    # ## Initialize result dictionary
    # trial_results = {'cluster_matching_score': score}

    
    # ## Plot the contingency matrix
    # cm = plottu.compute_contingency(adata, grn_adata_sub, config, col_adata = 'spectral', col_grn_adata = 'spectral_remap')

    
    # ## Save plots
    # sc.pl.umap(grn_adata_sub, color = ['spectral_remap'], title = 'GRN based embedding', save = '_grn.pdf')
    # sc.pl.umap(adata, color = ['spectral'], title = 'GEX based embedding', save = '_gex.pdf')
    # plottu.plot_differential_expression(grn_adata_sub, column = 'spectral_remap', suffix = 'grn.pdf')
    # plottu.plot_differential_expression(adata, column='spectral', suffix='gex.pdf')
    
    
    # ## Save edgelist
    # #summary= post.create_all_summaries(grn_adata, grn_adata_sub, rel, cluster_col = 'spectral')
    # #summary.to_csv(op.join(config['results']['output_directory'],config['results']['grn']), sep='\t', index = False)

    # ## Benchmark
    # net = pd.read_csv(op.join(op.dirname(config['data']['input_data']), 'net.tsv'), sep='\t')
    # trial_results = trial_results | eu.get_edge_overlap(grn_adata_sub, net,top_n=200)


    # with open(op.join(op.join(config['results']['output_directory']), 'performance.json'), 'w') as f:
    #     json.dump(trial_results, f)


def parse_args():
    parser = argparse.ArgumentParser(description="Argument parser with two configuration files.")
    
    parser.add_argument(
        '--config',
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
    
    return parser.parse_args()

if __name__ == "__main__":
    import argparse

    args = parse_args()
    print(f"Main Configuration File: {args.config}")
    print(f"Dataset Configuration File: {args.dataset_config}")

    
    config = ScGeneRAIConfig.read_yaml(args.config)
    dataset_config = DataSimulationConfig.read_yaml(args.dataset_config)

    
    run_scgenerai(config, dataset_config)