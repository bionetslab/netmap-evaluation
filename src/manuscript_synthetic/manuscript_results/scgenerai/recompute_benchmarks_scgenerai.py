import sys
sys.path.append('/data_nfs/og86asub/netmap/netmap-evaluation/')

import scanpy as sc
import os.path as op
import os
import pandas as pd
pd.set_option('display.float_format', lambda x: '%.3f' % x)

from netmap.utils.data_utils import *
from netmap.utils.tf_utils import *

from src.data_simulation.data_simulation_config import DataSimulationConfig

from src.evaluation.compute_metrics import build_augmented_network, create_forward_reverse
import src.evaluation.compute_metrics as compute_metrics

from src.utils import write_config


def load_networks(dataset_config):
    nets = [(op.basename(op.dirname(filename)), pd.read_csv(op.join(filename), sep=dataset_config.separator)) for filename in dataset_config.edgelist]
    off_net = [(op.basename(op.dirname(filename)), pd.read_csv(op.join(filename), sep=dataset_config.separator)) for filename in dataset_config.common_edges]
    augmented_nets = [(net[0], build_augmented_network(net[1])) for net in nets]
    return nets, off_net, augmented_nets

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
    
def run(grn_path, data_path, dataset_config_path, summary_output_dir):
    os.makedirs(summary_output_dir, exist_ok=True)

    # if op.isfile(op.join(summary_output_dir, 'scgenerai_clustering_score.json')):
    #     print(op.join(summary_output_dir, 'scgenerai_clustering_score.json'))
    #     return
    print('hop')
    dataset_config = DataSimulationConfig.read_yaml(dataset_config_path)
    nets, off_net, augmented_nets = load_networks(dataset_config)

    adata = sc.read_h5ad(data_path)
    grn_adata = sc.read_h5ad(grn_path)
    grn_adata.obs['grn'] = pd.Categorical(adata.obs['grn'])

    grn_ads = {'scgenerai': grn_adata}
    overlaps_ungrouped, collect_results = compute_metrics.compute_metrics(
        grn_ads=grn_ads, nets=nets, augmented_nets=augmented_nets, global_nets=off_net,
        group_key=dataset_config.group_key, group_by_target=False, aggregate=False
    )
    overlaps_ungrouped = process_results(overlaps_ungrouped)

    overlaps_ungrouped.to_csv(op.join(summary_output_dir, f'scgenerai_overlaps_global_top_k.tsv'), sep='\t')

    aggregated_performance = compute_metrics.compute_aggregated_grn_result(grn_adata, nets, cluster_col='grn')
    aggregated_performance.to_csv(op.join(summary_output_dir, 'scgenerai_aggregated_performance_grn.tsv'), sep='\t')


    write_config(c=collect_results, file=op.join(summary_output_dir, 'scgenerai_clustering_score.json'))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Recompute benchmark metrics for an already-generated scGeneRAI GRN (no training).")

    parser.add_argument("-g", "--grn_path", type=str, help="Path to the previously generated scgenerai.h5ad", required=True)
    parser.add_argument("--dataset_config", type=str, required=True, help="Path to the dataset configuration yaml")
    parser.add_argument("--summary_output_dir", type=str, required=True, help="Directory to write the evaluation summary files to")
    parser.add_argument("--data_path", type=str, required=True, help="Path to the data.h5ad the GRN was inferred from")

    args = parser.parse_args()
    print('FLOP')
    run(
        grn_path=args.grn_path,
        data_path=args.data_path,
        dataset_config_path=args.dataset_config,
        summary_output_dir=args.summary_output_dir,
    )
