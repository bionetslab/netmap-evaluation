#!/bin/bash

# Define the base directories and paths
MANIFEST="src/manuscript_synthetic/manuscript_results/grnboost2/pixi.toml"
SCRIPT="src/manuscript_synthetic/manuscript_results/grnboost2/compute_metrics_preclustered.py"

CLUSTERED_NET_DIR="/data_nfs/og86asub/netmap/netmap-evaluation/data/clustered_network"
GRN_BASE_DIR="/data_nfs/og86asub/netmap/netmap-evaluation/results/final_benchmark/grnboost2"
CONFIG_BASE_DIR="/data_nfs/og86asub/netmap/netmap-evaluation/results/configurations/data_simulation"
SUMMARY_BASE_DIR="results/summaries_final_benchmark/grnboost2"

# Target configuration folders to loop over
CONFIGS=("config_easy" "config_noise" "config_three" "config_three_noise" "config_five" "config_ten")

for config in "${CONFIGS[@]}"; do
    # We look inside the GRN results directory to see which subfolders successfully finished
    config_grn_dir="$GRN_BASE_DIR/$config"
    
    if [ -d "$config_grn_dir" ]; then
        echo "Computing metrics for configuration: $config"
        
        for subfolder_path in "$config_grn_dir"/*; do
            if [ -d "$subfolder_path" ]; then
                # Extract the subfolder name (e.g., net_53_11196_...)
                subfolder=$(basename "$subfolder_path")
                
                grn_file="$subfolder_path/grn.tsv"
                config_yaml="$CONFIG_BASE_DIR/$config/${subfolder}.config.yaml"
                summary_out="$SUMMARY_BASE_DIR/$config/$subfolder"
                
                # Check that both the computed GRN network and its simulation config exist before evaluating
                if [ -f "$grn_file" ] && [ -f "$config_yaml" ]; then
                    echo "  -> Processing $subfolder"
                    
                    pixi run python "$SCRIPT" \
                        --clustered_network_dir "$CLUSTERED_NET_DIR" \
                        --summary_output_dir "$summary_out" \
                        --grn "$grn_file" \
                        --dataset_config "$config_yaml"
                else
                    [ ! -f "$grn_file" ] && echo "  [Warning] Missing grn.tsv for $subfolder"
                    [ ! -f "$config_yaml" ] && echo "  [Warning] Missing expected YAML configuration: $config_yaml"
                fi
            fi
        done
    else
        echo "[Warning] GRN directory not found for config: $config"
    fi
done

echo "All metrics computation completed!"