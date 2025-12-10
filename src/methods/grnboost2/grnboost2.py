import os
import sys

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

import pandas as pd
sys.path.append('/data_nfs/og86asub/netmap/netmap-evaluation/')

from src.data_simulation.data_simulation_config import DataSimulationConfig
from grnboost2_config import GRNBoost2Config
import signifikante

import pandas as pd
import numpy as np
from signifikante.algo import grnboost2
from signifikante.utils import load_tf_names
import warnings
from src.utils import write_config

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



def run_grnboost_per_cluster(config):

    start = time.monotonic()
    os.makedirs(config.output_directory, exist_ok = True)
    # copy configuration file to the folder
    config.write_yaml(yaml_file=op.join(config.output_directory, 'config.yaml'))
    

    adata = sc.read_h5ad(config.input_data)
    if config.tf_only:
        tf_names = pd.read_csv(config.transcription_factors, header=None)
        tf_names.columns = ['tf_names']
        tf_names = list(tf_names.tf_names)
    else:
        tf_names = list(adata.var.index)
        print(tf_names)
        

    if 'grn' not in adata.obs.columns:
        raise ValueError("The 'grn' key is not in adata.obs.")
    
    # Initialize an empty list to store the GRNs for each cluster
    grn_list = []
    
    clusters = adata.obs['grn'].unique()
    

    start_method = time.monotonic()
    for cluster in clusters:
        
        adata_cluster = adata[adata.obs['grn'] == cluster, :].copy()
        
        expr_matrix = pd.DataFrame(adata_cluster.X,
                                   index=adata_cluster.obs_names,
                                   columns=adata_cluster.var_names)
        
        # 2. Run GRNBoost2
        # `verbose=False` to keep the output clean
        grn_df = grnboost2(expression_data=expr_matrix,
                           tf_names=tf_names,
                           verbose=False,
                           target_names = 'all')
        
        # 3. Add the cluster variable and save the GRN
        grn_df['grn'] = cluster
        grn_list.append(grn_df)
    
    # Concatenate all the individual GRN dataframes into a single dataframe
    final_grn_df = pd.concat(grn_list, ignore_index=True)
    time_method = time.monotonic()-start_method

    final_grn_df.to_csv(op.join(config.output_directory, config.grn))
    
    time_elapsed = time.monotonic()-start
    write_config({'time_elapsed_total': time_elapsed, 'time_elapsed_grnboost2': time_method}, file=op.join(config.output_directory, 'results.yaml'))
    return final_grn_df


if __name__ == "__main__":
    import argparse

    args = parse_args()
    print(f"Main Configuration File: {args.config}")
    print(f"Dataset Configuration File: {args.dataset_config}")

    
    config = GRNBoost2Config.read_yaml(args.config)
    #dataset_config = DataSimulationConfig.read_yaml(args.dataset_config)

    
    run_grnboost_per_cluster(config)