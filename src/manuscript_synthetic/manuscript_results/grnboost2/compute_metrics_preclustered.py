import pandas as pd
import numpy as np

import sys
sys.path.append('/data_nfs/og86asub/netmap/netmap-evaluation/')
from src.data_simulation.data_simulation_config import DataSimulationConfig
import os.path as op
import os
from src.evaluation.compute_metrics import build_augmented_network, create_forward_reverse
import src.evaluation.compute_metrics
from src.utils import write_config


def calculate_recovered_edges(inferred_grn, gold_standard_grn, k_values):
    """
    Computes the percentage of recovered edges for increasing k top edges in an inferred GRN.

    Args:
        inferred_grn (str): Path to a CSV file for the inferred GRN.
                                 The file should have columns: 'regulator', 'target', 'score'.
        gold_standard_grn (str): Path to a CSV file for the gold standard GRN.
                                      The file should have columns: 'regulator', 'target'.
        k_values (list): A list of integers representing the number of top edges to consider.

    Returns:
        pd.DataFrame: A DataFrame with columns 'k' and 'percentage_recovered',
                      showing the recovery percentage for each k value.
    """
    # Sort inferred edges by score in descending order
    inferred_grn = inferred_grn.sort_values(by='importance', ascending=False).reset_index(drop=True)

    # Create a set of gold standard edges for efficient lookup
    gold_standard_edges = set(zip(gold_standard_grn['source'], gold_standard_grn['target']))
    total_gold_standard_edges = len(gold_standard_edges)

    # Filter k_values to not exceed the total number of inferred edges
    max_k = len(inferred_grn)

    # Initialize a list to store results
    results = []

    for k in k_values:
        top_k_inferred = inferred_grn.head(k)

        # Create a set of the top k inferred edges
        top_k_edges = set(zip(top_k_inferred['TF'], top_k_inferred['target']))

        # Find the intersection (recovered edges)
        recovered_edges = len(top_k_edges.intersection(gold_standard_edges))

        # Calculate percentage of recovered edges
        if total_gold_standard_edges > 0:
            percentage = (recovered_edges / total_gold_standard_edges)
        else:
            percentage = 0  # Avoid division by zero if gold standard is empty

        results.append({'n_top': k, 'percentage_recovered': percentage, 'tp': recovered_edges, 'gs_count': total_gold_standard_edges, 'pp': len(top_k_edges) })

    results = pd.DataFrame(results)
    return results


def compute_metric(net, nets, k_thresholds, per_target=True):
    recovery_rates = []
    for grn in net.grn.unique():
        for n in range(len(nets)):
            if per_target:
                recovery_rate = calculate_recovered_edges_per_target( net[net.grn == grn],nets[n], k_thresholds)
            else:
                recovery_rate = calculate_recovered_edges(net[net.grn == grn], nets[n], k_thresholds)
            if grn-1 == n:
                recovery_rate['type'] = 'on_target'
            else:
                recovery_rate['type'] = 'off_target'
            
            recovery_rate['grn'] = int(grn-1)
            recovery_rate['net'] = n
            recovery_rates.append(recovery_rate)
    recovery_rates = pd.concat(recovery_rates)

    return recovery_rates




def reformat_dataframe(recovery_rates, config_name):
        # Assuming your DataFrame is named df
    recovery_rates = recovery_rates.pivot_table(
        index=['n_top', 'net', 'grn', 'pp', 'gs_count'],
        columns='type',
        values=['percentage_recovered', 'tp' ]
    ).reset_index()
    #recovery_rates.columns = ['_'.join(col).strip() for col in recovery_rates.columns.values]
    recovery_rates['method'] = config_name
    return recovery_rates



def calculate_recovered_edges_per_target(inferred_grn, gold_standard_grn, k_values):
    # Sort inferred edges by score in descending ordergrn_adata2
    inferred_grn = inferred_grn.sort_values(by='importance', ascending=False).reset_index(drop=True)

    # Create a set of gold standard edges for efficient lookup
    gold_standard_edges = set(zip(gold_standard_grn['source'], gold_standard_grn['target']))
    total_gold_standard_edges = len(gold_standard_edges)

    # Filter k_values to not exceed the total number of inferred edges
    max_k = len(inferred_grn)
    effective_k_values = [k for k in k_values if k <= max_k]
    effective_k_values.append(inferred_grn.shape[0])

    # Initialize a list to store results
    results = []

    for k in effective_k_values:
        
        top_k_inferred = inferred_grn.groupby('target').apply(lambda x: x.nlargest(k, 'importance')).reset_index(drop=True)
        # Create a set of the top k inferred edges
        top_k_edges = set(zip(top_k_inferred['TF'], top_k_inferred['target']))

        # Find the intersection (recovered edges)
        recovered_edges = len(top_k_edges.intersection(gold_standard_edges))

        # Calculate percentage of recovered edges
        if total_gold_standard_edges > 0:
            percentage = (recovered_edges / total_gold_standard_edges)
        else:
            percentage = 0  # Avoid division by zero if gold standard is empty

        results.append({'n_top': k, 'percentage_recovered': percentage})

    results = pd.DataFrame(results)

    return results




if  __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description="Argument parser with two configuration files.")
    
    parser.add_argument(
        '--dataset_config',
        type=str,
        default='dataset_configuration.json',
        help='Path to the dataset configuration file (default: dataset_configuration.json)'
    )


    parser.add_argument(
    "--clustered_network_dir",  # name on the CLI - drop the `--` for positional/required parameters
    nargs="*",  # 0 or more values expected => creates a list
    type=str,
    default=None,  # default if nothing is provided
    )

    parser.add_argument(
    "--summary_output_dir",  # name on the CLI - drop the `--` for positional/required parameters
    type=str,
    default=None,  # default if nothing is provided
    )

    parser.add_argument(
    "--grn",  # name on the CLI - drop the `--` for positional/required parameters
    type=str,
    default=None,  # default if nothing is provided
    )

    args = parser.parse_args()


    dataset_config = DataSimulationConfig.read_yaml(args.dataset_config)

    # read network files
    nets = [pd.read_csv(filename, sep=dataset_config.separator) for filename in dataset_config.edgelist]

    off_net = [(op.basename(op.dirname(filename)), pd.read_csv( filename, sep=dataset_config.separator)) for filename in dataset_config.common_edges]

    # Augmented net contains edges between genes which are controlled by the same transcription factor
    #augmented_nets = [(net[0], build_augmented_network(net[1])) for net in nets]



    results_per_target_collect = []
    results_global_collect = []

    net = pd.read_csv(args.grn)
    # k_thresholds = [10, 20, 50, 100, 200, 500, 5000, 7500, 10000]
    
    top_edges = [0.001, 0.01, 0.05, 0.1, 0.2, 0.25, 0.5, 0.75, 1.0]

    ## TODO: is it inflated?? 
    k_thresholds = [int(np.round(net.shape[0] * t)) for t in top_edges]



    results_global = compute_metric(net, nets, k_thresholds, per_target=False)
    #results_global = reformat_dataframe(results_global, c)
    #print(results_global)
    #results_global['config'] = c
    results_global_collect.append(results_global)


    # k_thresholds = [1,2, 3, 4 , 5, 10, 15, 20, 25, 50 , 75, 100]

    # results_per_target = compute_metric(net, nets, k_thresholds, per_target=True)
    # results_per_target = reformat_dataframe(results_per_target, c)

    # results_per_target_collect.append(results_per_target)

    results_global = pd.concat(results_global_collect)
    # results_per_target = pd.concat(results_per_target_collect)
    dataset_type = op.basename(op.dirname(args.dataset_config))
    outdir = op.join(args.summary_output_dir, dataset_type, dataset_config.dataset_id)

    os.makedirs(outdir, exist_ok = True)
    #write_config(collect_results, file=op.join(outdir, 'results.yaml'))
    # results_per_target.to_csv(op.join(outdir, 'preclustered_overlaps_per_target.tsv'), sep='\t')
    results_global.to_csv(op.join(outdir, 'preclustered_overlaps_global_top_k.tsv'), sep='\t')
        




    
