#!/bin/bash

# Define the base directories and paths
SCRIPT="src/manuscript_synthetic/manuscript_results/scgenerai/recompute_benchmarks_scgenerai.py"

DATA_BASE_DIR="/data_nfs/og86asub/netmap/netmap-evaluation/data/simulated_data_final_benchmark"
CONFIG_BASE_DIR="/data_nfs/og86asub/netmap/netmap-evaluation/results/configurations/data_simulation"

# Already-generated scgenerai.h5ad files from benchmark_scgenerai.sh
OUTPUT_BASE_DIR="results/final_benchmark/scgenerai"
# Where the recomputed benchmark summaries are written
SUMMARY_BASE_DIR="/data_nfs/og86asub/netmap/netmap-evaluation/results/summaries_final_benchmark/scgenerai"

# Target configuration folders to loop over
CONFIGS=("config_noise" "config_three" "config_three_noise")
CONFIGS=("config_easy")


for config in "${CONFIGS[@]}"; do
    config_data_dir="$DATA_BASE_DIR/$config"

    if [ -d "$config_data_dir" ]; then
        echo "Recomputing scgenerai benchmarks for configuration: $config"

        for subfolder_path in "$config_data_dir"/*; do
            if [ -d "$subfolder_path" ]; then
                # Extract the subfolder name (e.g., net_53_11196_...)
                subfolder=$(basename "$subfolder_path")

                data_file="$subfolder_path/data.h5ad"
                config_yaml="$CONFIG_BASE_DIR/$config/${subfolder}.config.yaml"
                grn_file="$OUTPUT_BASE_DIR/$config/$subfolder/scgenerai.h5ad"

                # Dynamically construct the summary directory matching the structure
                summary_out="$SUMMARY_BASE_DIR/$config/$subfolder"

                # Verify required inputs exist before executing
                if [ -f "$data_file" ] && [ -f "$config_yaml" ] && [ -f "$grn_file" ]; then
                    echo "  -> Processing $subfolder"
                    echo $SCRIPT
                    python "$SCRIPT" \
                        -g "$grn_file" \
                        --dataset_config "$config_yaml" \
                        --summary_output_dir "$summary_out" \
                        --data_path "$data_file"
                else
                    [ ! -f "$data_file" ] && echo "  [Warning] Missing data.h5ad for $subfolder"
                    [ ! -f "$config_yaml" ] && echo "  [Warning] Missing expected YAML configuration: $config_yaml"
                    [ ! -f "$grn_file" ] && echo "  [Warning] Missing scgenerai.h5ad output: $grn_file"
                fi
            fi
        done
    else
        echo "[Warning] Data directory not found for config: $config"
    fi
done

echo "All scgenerai benchmark recomputation completed!"
