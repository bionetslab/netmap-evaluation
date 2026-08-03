import os
import sys
import pandas as pd
import numpy as np
import anndata as ad
import time


import warnings
warnings.filterwarnings("ignore")
import os.path as op

sys.path.append('/data_nfs/og86asub/netmap/netmap-evaluation/')

from src.data_simulation.data_simulation_config import DataSimulationConfig
import signifikante
from signifikante.algo import grnboost2
from signifikante.utils import load_tf_names
from src.utils import write_config



def run_grnboost_per_cluster(input_data_file, transcription_factors, output_directory, tf_only=False, outfile = 'grn.tsv'):

    start = time.monotonic()
    os.makedirs(output_directory, exist_ok = True)
    # copy configuration file to the folder
    #config.write_yaml(yaml_file=op.join(output_directory, 'config.yaml'))
    

    adata = ad.read_h5ad(input_data_file)

    if tf_only:
        tf_names = pd.read_csv(transcription_factors, header=None)
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

    final_grn_df.to_csv(op.join(output_directory, outfile))
    
    time_elapsed = time.monotonic()-start
    write_config({'time_elapsed_total': time_elapsed, 'time_elapsed_grnboost2': time_method}, file=op.join(output_directory, 'results.yaml'))
    return final_grn_df


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run data generation with specified configuration")
    parser.add_argument("-d", "--data_file", type=str, help="dataset id", required=True)
    parser.add_argument("-t", "--tf_file", type=str, help="dataset id", required=True)
    parser.add_argument("-o", "--output_dir", type=str, help="dataset id", required=True)

    args = parser.parse_args()

    run_grnboost_per_cluster(input_data_file=args.data_file, transcription_factors = args.tf_file, output_directory = args.output_dir)