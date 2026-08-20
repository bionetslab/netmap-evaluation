#!/bin/bash
# Download the 10x PBMC multiome dataset and LINGER's matching reference data
# (PBMC.md tutorial: https://github.com/Durenlab/LINGER/blob/main/docs/PBMC.md).
#
# Usage:
#   ./download_pbmc_multiome_for_linger.sh -o <data_outdir> [-g <grn_bulk_dir>] [-s]
#
#   -o  directory to download/extract the 10x feature matrix + cell-type labels into
#   -g  directory to download/extract LINGER's bulk reference GRN data into
#       (large, shared across runs; omit -g to skip this download)
#   -s  skip the bulk reference download even if -g is given (e.g. already downloaded)
set -euo pipefail

DATA_OUTDIR=""
GRN_DIR=""
SKIP_BULK_REF=0

while getopts "o:g:sh" opt; do
    case "$opt" in
        o) DATA_OUTDIR="$OPTARG" ;;
        g) GRN_DIR="$OPTARG" ;;
        s) SKIP_BULK_REF=1 ;;
        h)
            echo "Usage: $0 -o <data_outdir> [-g <grn_bulk_dir>] [-s]"
            exit 0
            ;;
        *)
            echo "Usage: $0 -o <data_outdir> [-g <grn_bulk_dir>] [-s]"
            exit 1
            ;;
    esac
done

if [[ -z "$DATA_OUTDIR" ]]; then
    echo "ERROR: -o <data_outdir> is required" >&2
    exit 1
fi

mkdir -p "$DATA_OUTDIR"

echo "Downloading 10x PBMC multiome feature matrix..."
wget -O "$DATA_OUTDIR/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.tar.gz" \
    https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.tar.gz
tar -xzf "$DATA_OUTDIR/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.tar.gz" -C "$DATA_OUTDIR"
gzip -df "$DATA_OUTDIR"/filtered_feature_bc_matrix/*.gz

echo "Downloading LINGER's PBMC cell-type label file..."
LABEL_COOKIES="$(mktemp)"
wget --load-cookies "$LABEL_COOKIES" \
    "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies "$LABEL_COOKIES" --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=17PXkQJr8fk0h90dCkTi3RGPmFNtDqHO_' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=17PXkQJr8fk0h90dCkTi3RGPmFNtDqHO_" \
    -O "$DATA_OUTDIR/PBMC_label.txt"
rm -f "$LABEL_COOKIES"

if [[ -n "$GRN_DIR" && "$SKIP_BULK_REF" -eq 0 ]]; then
    mkdir -p "$GRN_DIR"
    echo "Downloading LINGER's bulk reference GRN data (this is large)..."
    BULK_COOKIES="$(mktemp)"
    wget --load-cookies "$BULK_COOKIES" \
        "https://drive.usercontent.google.com/download?export=download&confirm=$(wget --quiet --save-cookies "$BULK_COOKIES" --keep-session-cookies --no-check-certificate 'https://drive.usercontent.google.com/download?id=1jwRgRHPJrKABOk7wImKONTtUupV7yJ9b' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1jwRgRHPJrKABOk7wImKONTtUupV7yJ9b" \
        -O "$GRN_DIR/data_bulk.tar.gz"
    tar -xzf "$GRN_DIR/data_bulk.tar.gz" -C "$GRN_DIR"
    rm -f "$BULK_COOKIES"
fi

echo ""
echo "Done."
echo "  10x feature matrix: $DATA_OUTDIR/filtered_feature_bc_matrix"
echo "  Cell-type labels:   $DATA_OUTDIR/PBMC_label.txt"
if [[ -n "$GRN_DIR" ]]; then
    echo "  GRNdir:              $GRN_DIR/data_bulk/"
fi
