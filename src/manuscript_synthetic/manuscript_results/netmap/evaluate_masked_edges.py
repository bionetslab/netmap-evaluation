import sys
sys.path.append('/data_nfs/og86asub/netmap/netmap-evaluation/')

import scanpy as sc
import os.path as op
import os
import glob
import numpy as np
import pandas as pd
pd.set_option('display.float_format', lambda x: '%.3f' % x)

import yaml

from netmap.utils.misc import write_config
from netmap.utils.data_utils import *
from netmap.utils.tf_utils import *

from netmap.downstream.clustering import unify_group_labelling
from netmap.masking.internal import add_neighbourhood_expression_mask, add_cluster_based_candidate_edges

from src.data_simulation.data_simulation_config import DataSimulationConfig

from src.evaluation.compute_metrics import build_augmented_network, create_forward_reverse
import src.evaluation.compute_metrics as compute_metrics
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

    # Use numpy's trapezoidal rule (numerical integration) to compute the area:
    # Area = Integral of (y - x) dx over the positive region
    area = si.trapezoid(y=difference, x=group['n_top'])
    return area


def load_networks(dataset_config):
    nets = [(op.basename(op.dirname(filename)), pd.read_csv(op.join(filename), sep=dataset_config.separator)) for filename in dataset_config.edgelist]
    off_net = [(op.basename(op.dirname(filename)), pd.read_csv(op.join(filename), sep=dataset_config.separator)) for filename in dataset_config.common_edges]
    augmented_nets = [(net[0], build_augmented_network(net[1])) for net in nets]
    return nets, off_net, augmented_nets


def select_cluster_specific_edges(adata, grn_adata, group_key, leiden_resolution, mask_threshold):
    """Cluster the GRN attributions, mask edges by neighbourhood co-expression,
    and keep only edges that are candidates in at least one cluster.
    """
    sc.tl.pca(grn_adata, svd_solver='randomized', zero_center=False)
    sc.pp.neighbors(grn_adata, n_neighbors=50)
    sc.tl.leiden(grn_adata, resolution=leiden_resolution)

    score, mapping = unify_group_labelling(adata, grn_adata, group_key, 'leiden', True)

    add_neighbourhood_expression_mask(adata, grn_adata, strict=False, layer='X')
    grn_adata = add_cluster_based_candidate_edges(grn_adata, cluster_column='leiden_remap', threshold=mask_threshold)

    index_list = np.where(grn_adata.var['candidate_edge'])[0]
    grn_adata_masked = grn_adata[:, index_list].copy()

    return grn_adata_masked, score


def evaluate_one(grn_path, adata, nets, off_net, augmented_nets, group_key, sumout, leiden_resolution, mask_threshold):
    prefix = op.basename(grn_path).replace('_grn.h5ad', '')

    if op.isfile(op.join(sumout, f'{prefix}_masked_clustering_score.json')):
        print(op.join(sumout, f'{prefix}_masked_clustering_score.json'))
        return

    grn_adata = sc.read_h5ad(grn_path)

    n_edges_total = grn_adata.shape[1]
    grn_adata_masked, clustering_score = select_cluster_specific_edges(
        adata, grn_adata, group_key=group_key, leiden_resolution=leiden_resolution, mask_threshold=mask_threshold
    )
    n_edges_kept = grn_adata_masked.shape[1]

    grn_ads = {prefix: grn_adata_masked}

    overlaps_ungrouped, collect_results = compute_metrics.compute_metrics(
        grn_ads=grn_ads, nets=nets, augmented_nets=augmented_nets, global_nets=off_net,
        group_key=group_key, group_by_target=False, aggregate=False
    )

    overlaps_ungrouped = process_results(overlaps_ungrouped)
    area_results = overlaps_ungrouped.groupby(['net', 'target', 'direction', 'method', 'net_type', 'layer']).apply(compute_area_over_diagonal).rename('Area Over Diagonal')
    overlaps_ungrouped.to_csv(op.join(sumout, f'{prefix}_masked_overlaps_global_top_k.tsv'), sep='\t')

    aggregated_performance = compute_metrics.compute_aggregated_grn_result(grn_adata_masked, nets, cluster_col='leiden_remap')
    aggregated_performance.to_csv(op.join(sumout, f'{prefix}_masked_aggregated_performance.tsv'), sep='\t')

    aggregated_performance_grn = compute_metrics.compute_aggregated_grn_result(grn_adata_masked, nets, cluster_col='grn')
    aggregated_performance_grn.to_csv(op.join(sumout, f'{prefix}_masked_aggregated_performance_grn.tsv'), sep='\t')

    area_results = area_results.reset_index()
    area_results.to_csv(op.join(sumout, f'{prefix}_masked_area_over_diagnoal.tsv'), sep='\t')

    collect_results[prefix]['leiden_clustering_score'] = float(clustering_score)
    collect_results[prefix]['n_edges_total'] = int(n_edges_total)
    collect_results[prefix]['n_edges_kept'] = int(n_edges_kept)
    write_config(c=collect_results, file=op.join(sumout, f'{prefix}_masked_clustering_score.json'))


def run(output_dir, dataset_config_path, data_path, summary_output_dir, leiden_resolution=0.29, mask_threshold=0.5):
    os.makedirs(summary_output_dir, exist_ok=True)

    dataset_config = DataSimulationConfig.read_yaml(dataset_config_path)
    nets, off_net, augmented_nets = load_networks(dataset_config)

    adata = sc.read_h5ad(data_path)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.tl.pca(adata, svd_solver='randomized', zero_center=False)

    grn_paths = sorted(glob.glob(op.join(output_dir, '*_grn.h5ad')))
    if not grn_paths:
        print(f'No *_grn.h5ad files found in {output_dir}')
        return

    for grn_path in grn_paths:
        print(f'Evaluating {grn_path}')
        evaluate_one(
            grn_path=grn_path, adata=adata, nets=nets, off_net=off_net, augmented_nets=augmented_nets,
            group_key=dataset_config.group_key, sumout=summary_output_dir,
            leiden_resolution=leiden_resolution, mask_threshold=mask_threshold
        )


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Evaluate already-trained netmap GRNs on cluster-specific, mask-filtered edges.")

    parser.add_argument("-o", "--output_dir", type=str, help="Directory containing the previously generated *_grn.h5ad files", required=True)

    parser.add_argument(
        '--dataset_config',
        type=str,
        default='dataset_configuration.json',
        help='Path to the dataset configuration file (default: dataset_configuration.json)'
    )

    parser.add_argument(
        "--summary_output_dir",
        type=str,
        default=None,
        required=True,
        help="Directory to write the evaluation summary files to",
    )

    parser.add_argument(
        "--data_path",
        type=str,
        default=None,
        required=True,
        help="Path to the data.h5ad the GRNs were inferred from",
    )

    parser.add_argument("--leiden_resolution", type=float, default=0.29, help="Leiden resolution used to cluster the GRN attributions")
    parser.add_argument("--mask_threshold", type=float, default=0.5, help="Minimum per-cluster mask support fraction for an edge to be kept")

    args = parser.parse_args()

    run(
        output_dir=args.output_dir,
        dataset_config_path=args.dataset_config,
        data_path=args.data_path,
        summary_output_dir=args.summary_output_dir,
        leiden_resolution=args.leiden_resolution,
        mask_threshold=args.mask_threshold,
    )
