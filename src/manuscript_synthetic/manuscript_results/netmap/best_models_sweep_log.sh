#!/bin/bash

# Define the base directories and paths
SCRIPT="src/manuscript_synthetic/manuscript_results/netmap/best_models_training_log.py"

DATA_BASE_DIR="/data_nfs/og86asub/netmap/netmap-evaluation/data/simulated_data_final_benchmark"
CONFIG_BASE_DIR="/data_nfs/og86asub/netmap/netmap-evaluation/results/configurations/data_simulation"
CLUSTERED_NET_DIR="/data_nfs/og86asub/netmap/netmap-evaluation/data/clustered_network"

OUTPUT_BASE_DIR="results/final_benchmark/netmap"
SUMMARY_BASE_DIR="/data_nfs/og86asub/netmap/netmap-evaluation/results/summaries_final_benchmark/netmap"

# Target configuration folders to loop over
CONFIGS=("config_easy"  "config_five" "config_ten")
CONFIGS=("config_noise" "config_three" "config_three_noise")


for config in "${CONFIGS[@]}"; do
    config_data_dir="$DATA_BASE_DIR/$config"
    
    if [ -d "$config_data_dir" ]; then
        echo "Processing training logs for configuration: $config"
        
        for subfolder_path in "$config_data_dir"/*; do
            if [ -d "$subfolder_path" ]; then
                # Extract the subfolder name (e.g., net_53_11196_...)
                subfolder=$(basename "$subfolder_path")
                
                data_file="$subfolder_path/data.h5ad"
                config_yaml="$CONFIG_BASE_DIR/$config/${subfolder}.config.yaml"
                
                # Dynamically construct output and summary directories matching the structure
                out_dir="$OUTPUT_BASE_DIR/$config/$subfolder"
                summary_out="$SUMMARY_BASE_DIR/$config/$subfolder"
                
                # Verify required inputs exist before executing
                if [ -f "$data_file" ] && [ -f "$config_yaml" ]; then
                    echo "  -> Processing $subfolder"
                    
                    pixi run python "$SCRIPT" \
                        -o "$out_dir" \
                        --dataset_config "$config_yaml" \
                        --clustered_network_dir "$CLUSTERED_NET_DIR" \
                        --summary_output_dir "$summary_out" \
                        --data_path "$data_file"
                else
                    [ ! -f "$data_file" ] && echo "  [Warning] Missing data.h5ad for $subfolder"
                    [ ! -f "$config_yaml" ] && echo "  [Warning] Missing expected YAML configuration: $config_yaml"
                fi
            fi
        done
    else
        echo "[Warning] Data directory not found for config: $config"
    fi
done

echo "All training log processing completed!"