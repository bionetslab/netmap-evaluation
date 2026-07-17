import os
import sys
os.environ["CUDA_VISIBLE_DEVICES"]="1"


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
from netmap.utils.data_utils import *
from netmap.utils.tf_utils import *


from src.methods.scgenerai.scgenerai_config import ScGeneRAIConfig
from src.data_simulation.data_simulation_config import DataSimulationConfig

#Import cloned scgenerai/ can't be installed via pip, see issue on github
sys.path.append('/data_nfs/og86asub/netmap/scGeneRAI')
from scGeneRAI import scGeneRAI

from src.utils import write_config, split_index


def run_scgenerai(input_data_file, output_directory):

    start_total = time.monotonic()
    

    ## Load config and setup outputs
    os.makedirs(output_directory, exist_ok=True)
    
    #setup temp dir for scGeneRAI to save results
    temp_dir = op.join(output_directory, 'tmp')
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    os.environ['scGeneRAI_TEMPDIR'] = temp_dir

    
    ## load data
    adata = sc.read_h5ad(input_data_file)
    sc.pp.scale(adata)

    data = adata.X

    data_train_df = pd.DataFrame(data, index=adata.obs_names, columns=adata.var_names)
    data_test_df = data_train_df

    grn_col = [str(x) for x in adata.obs.grn]
    descriptors = pd.DataFrame({'grn': grn_col})
    print(descriptors)


 
    start = time.monotonic()
    model = scGeneRAI()

    if torch.cuda.is_available():
        device = torch.device('cuda')
    else:
        device = torch.device('cpu')

    print(data_test_df.shape)
    model.fit(data_train_df, model_depth=2, nepochs=100, lr=2e-2, batch_size=5, lr_decay = 0.99,  early_stopping = True, device_name = device)
    
    #Remove interactions between descriptors and variables
    # The documentation is wrong and to disable averaging we have to set LRPau false.
    model.predict_networks(data_test_df,  device_name = device, PATH = temp_dir, LRPau=False )
    

    model_time = time.monotonic()-start
    print(f'Elapsed time: {model_time}')

    res_files = os.listdir(op.join(temp_dir, 'results'))
    #res_data = pd.concat([pd.read_csv(temp_dir + '/results/' + res_file) for res_file in res_files])
    res_data = pd.concat([pd.read_csv(op.join(temp_dir,'results', res_file)).assign(cell_name=res_file.split('_')[1].replace('.csv', '')) for res_file in res_files])
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
    grn_adata = split_index(grn_adata)
    grn_adata.write_h5ad(op.join(output_directory,'scgenerai.h5ad'))


    time_elapsed_total = time.monotonic()-start_total
    write_config({'time_elapsed_total': time_elapsed_total, 'time_elapsed_scgenerai': model_time}, file=op.join(output_directory, 'results.yaml'))

  

def parse_args():
    parser = argparse.ArgumentParser(description="Argument parser with two configuration files.")
    
    parser.add_argument(
        '--input_data',
        type=str,
        default='configuration.json',
        help='Path to the main configuration file (default: configuration.json)'
    )
    parser.add_argument(
        '--output_dir',
        type=str,
        default='dataset_configuration.json',
        help='Path to the dataset configuration file (default: dataset_configuration.json)'
    )
    
    return parser.parse_args()

if __name__ == "__main__":
    import argparse

    args = parse_args()

    
    run_scgenerai(args.input_data, args.output_dir)