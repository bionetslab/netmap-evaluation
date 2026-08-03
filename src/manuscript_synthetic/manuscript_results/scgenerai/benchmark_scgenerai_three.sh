#!/bin/bash

# Define the base directories and paths
SCRIPT="src/manuscript_synthetic/manuscript_results/scgenerai/run_scgenerai.py"

DATA_BASE_DIR="/data_nfs/og86asub/netmap/netmap-evaluation/data/simulated_data_final_benchmark"
CONFIG_BASE_DIR="/data_nfs/og86asub/netmap/netmap-evaluation/results/configurations/data_simulation"
CLUSTERED_NET_DIR="/data_nfs/og86asub/netmap/netmap-evaluation/data/clustered_network"

OUTPUT_BASE_DIR="results/final_benchmark/scgenerai"
SUMMARY_BASE_DIR="/data_nfs/og86asub/netmap/netmap-evaluation/results/summaries_final_benchmark/scgenerai"

# Target configuration folders to loop over
CONFIGS=("config_easy")
CONFIGS=("config_noise" "config_three")
CONFIGS=("config_easy" "config_three_noise")



for config in "${CONFIGS[@]}"; do
    config_data_dir="$DATA_BASE_DIR/$config"

    if [ -d "$config_data_dir" ]; then
        echo "Running scgenerai for configuration: $config"

        for subfolder_path in "$config_data_dir"/*; do
            if [ -d "$subfolder_path" ]; then
                # Extract the subfolder name (e.g., net_53_11196_...)
                subfolder=$(basename "$subfolder_path")

                input_data="$subfolder_path/data.h5ad"
                config_yaml="$CONFIG_BASE_DIR/$config/${subfolder}.config.yaml"

                # Dynamically construct output paths preserving the correct hierarchy
                out_dir="$OUTPUT_BASE_DIR/$config/$subfolder"
                summary_out="$SUMMARY_BASE_DIR/$config/$subfolder"
                output_h5ad="$out_dir/scgenerai.h5ad"

                # Skip if this dataset has already been processed
                if [ -f "$output_h5ad" ]; then
                    echo "  -> Skipping $subfolder (scgenerai.h5ad already exists)"
                    continue
                fi

                # Verify required inputs exist before executing
                if [ -f "$input_data" ] && [ -f "$config_yaml" ]; then
                    echo "  -> Processing $subfolder"

                    python "$SCRIPT" \
                        --output_dir "$out_dir" \
                        --dataset_config "$config_yaml" \
                        --clustered_network_dir "$CLUSTERED_NET_DIR" \
                        --summary_output_dir "$summary_out" \
                        --input_data "$input_data"
                else
                    [ ! -f "$input_data" ] && echo "  [Warning] Missing data.h5ad for $subfolder"
                    [ ! -f "$config_yaml" ] && echo "  [Warning] Missing expected YAML configuration: $config_yaml"
                fi
            fi
        done
    else
        echo "[Warning] Data directory not found for config: $config"
    fi
done

echo "All scgenerai processing completed!"
