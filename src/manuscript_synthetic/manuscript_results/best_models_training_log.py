import sys
sys.path.append('/data_nfs/og86asub/netmap/netmap-evaluation/')

import scanpy as sc
import time 
import os.path as op
import os
import numpy as np
import pandas as pd
pd.set_option('display.float_format', lambda x: '%.3f' % x)
import scipy.sparse as scs

import torch
import yaml


from netmap.utils.misc import write_config


from netmap.utils.data_utils import *
from netmap.utils.tf_utils import *

from netmap.model.train_model import create_model_zoo
import netmap.grn.inferrence as inferrence
from src.data_simulation.data_simulation_config import DataSimulationConfig

from compute_metrics import build_augmented_network, create_forward_reverse
import compute_metrics
import scipy.integrate as si



def process_results(all_overlaps):
    all_overlaps['configuration'] = all_overlaps['method']+'_'+all_overlaps['layer']
    # this should be changed if we use different sized popultaitons
    all_overlaps['factor'] = all_overlaps['mean_cell_count']/500
    all_overlaps['factor'] = all_overlaps['factor'].apply(lambda x: min(1, x))
    all_overlaps['weighted_overlap'] = all_overlaps['avg_edges_recovered']*(all_overlaps['factor'])

    all_overlaps = all_overlaps[~all_overlaps.avg_edges_recovered.isna()]
    all_overlaps = all_overlaps[~all_overlaps.mean_cell_count.isna()]
    all_overlaps['percentage_overlap'] = all_overlaps['weighted_overlap']/all_overlaps['gold_standard_edges']

    all_overlaps['precision'] = all_overlaps['weighted_overlap']/(all_overlaps['avg_edges_recovered'] + (all_overlaps['max_possible_edges']*all_overlaps['n_top']))
    return all_overlaps

def compute_area_over_diagonal(group):
    # Ensure data is sorted by x-axis ('n_top') for correct trapezoidal integration
    group = group.sort_values(by='n_top')
    
    # Calculate the vertical distance from the diagonal: f(x) = y - x
    difference = group['percentage_overlap'] - group['n_top']
    
    # We only want the area where the line is above the diagonal (difference > 0).
    # Set negative differences to 0.
    positive_difference = np.maximum(difference, 0)
    
    # Use numpy's trapezoidal rule (numerical integration) to compute the area: 
    # Area = Integral of (y - x) dx over the positive region
    area = si.trapezoid(y=difference, x=group['n_top'])
    return area


def read_config(file):
    with open(file, "r") as f:
        config = yaml.safe_load(f)
    return config


if __name__== '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Hello.")

    parser.add_argument("-n", "--network", type=str, help="dataset id", required=True)
    parser.add_argument("-d", "--data_config", type=str, help="dataset id", required=True)
    args = parser.parse_args()


    #net = 'net_84_10865_net_88_10937_net_90_11013'
    net = args.network
    dat_conf = args.data_config


    
    dataset_config = f"/data_nfs/og86asub/netmap/netmap-evaluation/results/configurations/data_simulation/{dat_conf}/{net}.config.yaml"
    data_path = f"/data_nfs/og86asub/netmap/netmap-evaluation/data/simulated_data/{dat_conf}/{net}/data.h5ad"
    outdir = f'/data_nfs/og86asub/netmap/netmap-evaluation/results/summaries_best_models_log3/{net}'

    os.makedirs(outdir, exist_ok=True)

    if op.isfile(op.join(outdir, f'{net}.yaml')):
        print('Processing has begun')
        if op.isfile(op.join(outdir, f'runtime_all.yaml')):
            print('Processing has finished')
            exit(0)
        
    

    write_config(c = {'network': net}, file= op.join(outdir, f'{net}.yaml'))


    dataset_config =DataSimulationConfig.read_yaml(dataset_config)

    # read network files
    # Augmented net contains edges between genes which are controlled by the same transcription factor
    nets = [(op.basename(op.dirname(filename)), pd.read_csv(op.join(filename), sep=dataset_config.separator)) for filename in dataset_config.edgelist]
    off_net = [(op.basename(op.dirname(filename)), pd.read_csv(op.join(filename), sep=dataset_config.separator)) for filename in dataset_config.common_edges]
    augmented_nets = [(net[0], build_augmented_network(net[1])) for net in nets]

    print(nets)

    #Create the networks for forward and reverse nets
    forward_reverse_nets = [(net[0], create_forward_reverse(net[1])) for net in nets]
    forward_reverse_nets_augmented = [(net[0], build_augmented_network(net[1])) for net in forward_reverse_nets]
    forward_reverse_off = [(net[0], create_forward_reverse(net[1])) for net in off_net]

    combined_net = [n[1] for n in nets] + [n[1] for n in off_net]
    combined_net = np.concatenate(combined_net)
    combined_net = pd.DataFrame(combined_net)
    combined_net.columns = ['source', 'target']
    combined_net['weight'] = 1

    # load data
    adata = sc.read_h5ad(data_path)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    print(adata.shape)
    ## Get the data matrix from the CustumAnndata obeject

    gene_names = np.array(adata.var.index)
    model_start = time.monotonic()

    # if config.layer == 'counts':
    #     data_tensor = adata.layers['counts']
    # else:
    data_tensor = adata.X

    if scs.issparse(data_tensor):
        data_tensor = torch.tensor(data_tensor.todense(), dtype=torch.float32)
    else:
        data_tensor = torch.tensor(data_tensor, dtype=torch.float32)

    start_time_all = time.monotonic()
    hyper = {
        '64_1':[64]
    }    

    dropout_percentage = [0.1]
    model_type = ['NegativeBinomialAutoencoder']
    background_type = ['zeros']
    attribution_raw = [False]
    
    time_collector = {}
    for hd in hyper.keys():   
        for dp in dropout_percentage:
            for mt in model_type:
                if op.isfile(op.join(outdir, f'{hd}_{dp}_{mt}_training_time.json')):
                    print(op.join(outdir, f'{hd}_{dp}_{mt}_training_time.json'))
                    continue
                start_training = time.monotonic() 
                model_zoo = create_model_zoo(data_tensor,  n_models=10, n_epochs=10000, model_type=mt, dropout_rate=dp, hidden_dim = hyper[hd] )
                end_training = time.monotonic()

                if len(model_zoo)==0:
                    print(f'Model failed to train {hd}_{dp}_{mt}, passing on')
                    continue
                
                write_config(c = {'training_time': end_training-start_training}, file= op.join(outdir, f'{hd}_{dp}_{mt}_training_time.json'))

                for xai in ['GradientShap']:
                    
                    if xai == 'GradientShap':
                        bg_modes = background_type
                    else:
                        bg_modes = ['zeros'] #Dummy (will not be used)

                    for bg in bg_modes:
                        if op.isfile(op.join(outdir, f'{xai}_{hd}_{dp}_{mt}_{bg}_clustering_score.json')):
                            # Check if output for this configuration exists and if so continue.
                            print(op.join(outdir, f'{xai}_{hd}_{dp}_{mt}_{bg}_clustering_score.json'))
                            break
                        for raw in attribution_raw:
                            start_inference = time.monotonic()
                            grn_adata2 = inferrence.inferrence(model_zoo, data_tensor.cuda(), gene_names, xai_method=xai, background_type = 'zeros', backing_file=None)
                            end_inference = time.monotonic()
                            if grn_adata2 is None:
                                continue
                            grn_adata2.obs['grn'] = pd.Categorical(adata.obs['grn'])
                            grn_adata2.write_h5ad(op.join(outdir, f'{xai}_{hd}_{dp}_{mt}_{bg}_{raw}_grn.h5ad'))
                            grn_ads = {f'{xai}_{hd}_{dp}_{mt}_{bg}' : grn_adata2}



                            overlaps_ungrouped, collect_results = compute_metrics.compute_metrics(grn_ads=grn_ads, nets= nets, augmented_nets=augmented_nets, global_nets=off_net, group_key=dataset_config.group_key, group_by_target=False, aggregate=False)
                            
                            # Add percentage and precision
                            overlaps_ungrouped = process_results(overlaps_ungrouped)
                            
                            area_results = overlaps_ungrouped.groupby(['net', 'target', 'direction', 'method', 'net_type', 'layer']).apply(compute_area_over_diagonal).rename('Area Over Diagonal')
                            overlaps_ungrouped.to_csv(op.join(outdir, f'{xai}_{hd}_{dp}_{mt}_{bg}_{raw}_overlaps_global_top_k.tsv'), sep='\t')

                                                                                    
                            aggregated_performance = compute_metrics.compute_aggregated_grn_result(grn_adata2, nets, cluster_col = 'spectral_remap' )
                            aggregated_performance.to_csv(op.join(outdir, f'{xai}_{hd}_{dp}_{mt}_{bg}_{raw}_aggregated_performance.tsv'), sep = '\t')
                            
                            aggregated_performance_scaled = compute_metrics.compute_aggregated_grn_result(grn_adata2, nets, cluster_col = 'grn' )
                            aggregated_performance_scaled.to_csv(op.join(outdir, f'{xai}_{hd}_{dp}_{mt}_{bg}_{raw}_aggregated_performance_grn.tsv'), sep = '\t')
                            
                            
                            area_results = area_results.reset_index()
                            area_results.to_csv(op.join(outdir, f'{xai}_{hd}_{dp}_{mt}_{bg}_{raw}_area_over_diagnoal.tsv'), sep='\t')

                            collect_results[f'{xai}_{hd}_{dp}_{mt}_{bg}']['training_time'] = end_training - start_training
                            collect_results[f'{xai}_{hd}_{dp}_{mt}_{bg}']['inference_time'] = end_inference - start_inference
                            write_config(c = collect_results, file= op.join(outdir, f'{xai}_{hd}_{dp}_{mt}_{bg}_{raw}_clustering_score.json'))

    end_time_all = time.monotonic()
    write_config(c = {'total_time': end_time_all-start_time_all}, file= op.join(outdir, f'runtime_all.yaml'))
