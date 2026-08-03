#!/bin/bash

# Define the base directories
DATA_BASE="data/simulated_data_final_benchmark"
OUTPUT_BASE="results/final_benchmark/grnboost2"
TF_FILE="data/annotation/allTFs_hg38.txt"
MANIFEST="src/manuscript_synthetic/manuscript_results/grnboost2/pixi.toml"
SCRIPT="src/manuscript_synthetic/manuscript_results/grnboost2/run_grnboost2.py"

# List of top-level configuration folders to iterate over
CONFIGS=("config_easy" "config_noise" "config_three" "config_three_noise" "config_five" "config_ten")

for config in "${CONFIGS[@]}"; do
    config_dir="$DATA_BASE/$config"
    
    # Check if the configuration folder actually exists before proceeding
    if [ -d "$config_dir" ]; then
        echo "Processing configuration: $config"
        
        # Loop through every subfolder inside this config folder
        for subfolder_path in "$config_dir"/*; do
            if [ -d "$subfolder_path" ]; then
                # Extract just the subfolder name (e.g., net_53_11196_...)
                subfolder=$(basename "$subfolder_path")
                
                # Check if the expected h5ad data file exists
                data_file="$subfolder_path/data.h5ad"
                if [ -f "$data_file" ]; then
                    echo "  -> Running subfolder: $subfolder"
                    
                    # Construct the matching output directory
                    out_dir="$OUTPUT_BASE/$config/$subfolder"
                    
                    # Execute the Pixi task
                    pixi run --manifest-path "$MANIFEST" python "$SCRIPT" \
                        --data_file "$data_file" \
                        --tf_file "$TF_FILE" \
                        --output_dir "$out_dir"
                else
                    echo "  [Warning] Missing data.h5ad in $subfolder_path"
                fi
            fi
        done
    else
        echo "[Warning] Configuration folder not found: $config_dir"
    fi
done

echo "All configurations completed!"