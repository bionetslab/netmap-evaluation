
import sys
import os.path as op
import os
sys.path.append('/data_nfs/og86asub/netmap/netmap-evaluation/')

from netmap.utils.netmap_config import NetmapConfig
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
import itertools
import warnings


import numpy as np
from collections import Counter
from netmap.masking.external import *
import time

from compute_metrics import build_augmented_network, create_forward_reverse
import compute_metrics
import scipy.integrate as si
import src.evaluation.compute_metrics as compute_metrics

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
    outdir = f'/data_nfs/og86asub/netmap/netmap-evaluation/results/summaries_scgenerai/{net}'
    grn_path = f"/data_nfs/og86asub/netmap/netmap-evaluation/results/scgenerai/config/{dat_conf}/{net}/grn_lrp.h5ad"

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

    # read grn file
    adata = sc.read_h5ad(data_path)
    grn_adata2 = sc.read_h5ad(grn_path)
                                
    if grn_adata2 is None:
        exit

    scgenerai_name = f'{dat_conf}_{net}'


    grn_adata2.obs['grn'] = pd.Categorical(adata.obs['grn'])

    grn_ads = {'scgenerai':grn_adata2}
    overlaps_ungrouped, collect_results = compute_metrics.compute_metrics(grn_ads=grn_ads, nets= nets, augmented_nets=augmented_nets, global_nets=off_net, group_key=dataset_config.group_key, group_by_target=False, aggregate=False)
    
    # Add percentage and precision
    overlaps_ungrouped = process_results(overlaps_ungrouped)
    
    overlaps_ungrouped.to_csv(op.join(outdir, f'{scgenerai_name}_overlaps_global_top_k.tsv'), sep='\t')



    overlaps_averaged, collect_results_avg = compute_metrics.compute_metrics(grn_ads=grn_ads, nets= forward_reverse_nets, augmented_nets=forward_reverse_nets_augmented, global_nets=forward_reverse_off, group_key=dataset_config.group_key, group_by_target=False)
    overlaps_averaged = process_results(overlaps_averaged)

    overlaps_averaged.to_csv(op.join(outdir, f'{scgenerai_name}_overlaps_global_top_k_fr.tsv'), sep='\t')



                                                            
    aggregated_performance = compute_metrics.compute_aggregated_grn_result(grn_adata2, nets, cluster_col = 'spectral_remap' )
    aggregated_performance.to_csv(op.join(outdir, f'{scgenerai_name}_aggregated_performance.tsv'), sep = '\t')
    
    sc.pp.normalize_total(grn_adata2)
    aggregated_performance_scaled = compute_metrics.compute_aggregated_grn_result(grn_adata2, nets, cluster_col = 'spectral_remap' )
    aggregated_performance_scaled.to_csv(op.join(outdir, f'{scgenerai_name}_aggregated_performance_scaled.tsv'), sep = '\t')
    
    
    write_config(c = collect_results, file= op.join(outdir, f'{scgenerai_name}_clustering_score.json'))

