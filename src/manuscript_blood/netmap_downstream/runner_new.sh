#!/usr/bin/env bash
set -euo pipefail

# ---------------------------------------------------------------------------
# Runner for run_netmap_analysis.py
# Mirrors the existing batch_analyse_output.py invocation style.
# ---------------------------------------------------------------------------
SCRIPT="src/manuscript_blood/netmap_downstream/run_netmap_analysis.py"
RESULTS_DIR="/data_nfs/og86asub/netmap/netmap-evaluation/netmap/case_studies/blood/"
ANALYSIS_DIR="/data_nfs/og86asub/netmap/netmap-evaluation/results/case_studies/blood/"
COLLECTRI="/data_nfs/og86asub/netmap/netmap-evaluation/data/input_network/collectri.tsv"
DATA_ROOT="/data_nfs/og86asub/netmap/netmap-evaluation/data/blood/reprocessed"

EXP_NAME="bd-rhap-rep2-markers"

python src/manuscript_blood/netmap_downstream/run_netmap_analysis.py \
    --adata_raw '/data_nfs/og86asub/netmap/netmap-evaluation/data/blood/reprocessed/bd-rhap-rep2/bd-rhap-rep2_with_markers.h5ad' \
    --collectri '/data_nfs/og86asub/netmap/netmap-evaluation/data/input_network/collectri.tsv' \
    --results_dir '/data_nfs/og86asub/netmap/netmap-evaluation/netmap/case_studies/blood/' \
    --experiment_name "$EXP_NAME" \
    --min_score 0.85 \
    --cluster_column 'leiden_remap' \
    --analysis_dir '/data_nfs/og86asub/netmap/netmap-evaluation/results/case_studies/blood/' \
    --top_per_source 10 20 30 40 50 \
    --min_reg_size 1 \
    --top_edge_percentage 0.1 \
    --n_neighbors 50 \
    --leiden_resolution 0.29 \
    --apply_ct_mapping \
    --extra_plots \
    --save_objects \
    --rerun




# EXP_NAME="10x-rep1-kallisto-cellbender-markers"
# python src/manuscript_blood/netmap_downstream/run_netmap_analysis.py \
#   --adata_raw "$DATA_ROOT/10x-rep1-kallisto-cellbender/10x-rep1-kallisto-cellbender_with_markers.h5ad" \
#   --collectri "$COLLECTRI" \
#   --results_dir "$RESULTS_DIR" \
#   --experiment_name "$EXP_NAME" \
#     --min_score 0.85 \
#     --cluster_column 'leiden_remap' \
#     --analysis_dir '/data_nfs/og86asub/netmap/netmap-evaluation/results/case_studies/blood/' \
#     --top_per_source 10 20 30 40 50 \
#     --min_reg_size 1 \
#     --top_edge_percentage 0.1 \
#     --n_neighbors 50 \
#     --leiden_resolution 0.29 \
#     --apply_ct_mapping \
#     --extra_plots \
#     --save_objects \
#     --rerun

# EXP_NAME="bd-rhap-rep1-markers"
# python src/manuscript_blood/netmap_downstream/run_netmap_analysis.py \
#   --adata_raw "$DATA_ROOT/bd-rhap-rep1/bd-rhap-rep1_with_markers.h5ad" \
#   --collectri "$COLLECTRI" \
#   --results_dir "$RESULTS_DIR" \
#   --experiment_name "$EXP_NAME" \
#     --min_score 0.85 \
#     --cluster_column 'leiden_remap' \
#     --analysis_dir '/data_nfs/og86asub/netmap/netmap-evaluation/results/case_studies/blood/' \
#     --top_per_source 10 20 30 40 50 \
#     --min_reg_size 1 \
#     --top_edge_percentage 0.1 \
#     --n_neighbors 50 \
#     --leiden_resolution 0.29 \
#     --apply_ct_mapping \
#     --extra_plots \
#     --save_objects \
#     --rerun

# EXP_NAME="10x-rep2-kallisto-cellbender-markers"
# echo $EXP_NAME
# python src/manuscript_blood/netmap_downstream/run_netmap_analysis.py \
#   --adata_raw "$DATA_ROOT/10x-rep2-kallisto-cellbender/10x-rep2-kallisto-cellbender_with_markers.h5ad" \
#   --collectri "$COLLECTRI" \
#   --results_dir "$RESULTS_DIR" \
#   --experiment_name "$EXP_NAME" \
#     --min_score 0.85 \
#     --cluster_column 'leiden_remap' \
#     --analysis_dir '/data_nfs/og86asub/netmap/netmap-evaluation/results/case_studies/blood/' \
#     --top_per_source 10 20 30 40 50 \
#     --min_reg_size 1 \
#     --top_edge_percentage 0.1 \
#     --n_neighbors 50 \
#     --leiden_resolution 0.29 \
#     --apply_ct_mapping \
#     --extra_plots \
#     --save_objects \
#     --rerun