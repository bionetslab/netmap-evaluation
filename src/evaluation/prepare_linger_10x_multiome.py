"""Turn a raw 10x multiome feature matrix + cell-type label file into the
paired, QC'd, barcode-synced RNA/ATAC h5ad pair that run_linger.py expects.

Mirrors PBMC.md tutorial steps 1-4 (data ingestion):
https://github.com/Durenlab/LINGER/blob/main/docs/PBMC.md
i.e. LingerGRN.preprocess.get_adata() (splits the combined 10x multiome matrix
into RNA/ATAC by feature type and attaches the cell-type label) followed by
min_genes/min_cells filtering and RNA/ATAC barcode synchronization.

Expects the output of download_pbmc_multiome_for_linger.sh (or any 10x
multiome filtered_feature_bc_matrix + matching cell-type label file with
'barcode_use'/'label' columns).
"""
import os
import os.path as op

import scanpy as sc
import scipy.io
import pandas as pd


def run(matrix_dir, label_file, outdir, min_genes=200, min_cells=3):
    from LingerGRN.preprocess import get_adata

    matrix = scipy.io.mmread(op.join(matrix_dir, 'matrix.mtx'))
    features = pd.read_csv(op.join(matrix_dir, 'features.tsv'), sep='\t', header=None)
    barcodes = pd.read_csv(op.join(matrix_dir, 'barcodes.tsv'), sep='\t', header=None)
    label = pd.read_csv(label_file, sep='\t', header=0)

    adata_RNA, adata_ATAC = get_adata(matrix, features, barcodes, label)

    sc.pp.filter_cells(adata_RNA, min_genes=min_genes)
    sc.pp.filter_genes(adata_RNA, min_cells=min_cells)
    sc.pp.filter_cells(adata_ATAC, min_genes=min_genes)
    sc.pp.filter_genes(adata_ATAC, min_cells=min_cells)

    selected_barcode = list(set(adata_RNA.obs['barcode']) & set(adata_ATAC.obs['barcode']))
    rna_idx = pd.DataFrame(range(adata_RNA.shape[0]), index=adata_RNA.obs['barcode'].values)
    adata_RNA = adata_RNA[rna_idx.loc[selected_barcode][0]].copy()
    atac_idx = pd.DataFrame(range(adata_ATAC.shape[0]), index=adata_ATAC.obs['barcode'].values)
    adata_ATAC = adata_ATAC[atac_idx.loc[selected_barcode][0]].copy()

    os.makedirs(outdir, exist_ok=True)
    adata_RNA.write(op.join(outdir, 'adata_RNA.h5ad'))
    adata_ATAC.write(op.join(outdir, 'adata_ATAC.h5ad'))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Ingest a raw 10x multiome feature matrix into the paired h5ad run_linger.py expects."
    )
    parser.add_argument("-m", "--matrix_dir", type=str, required=True,
                         help="10x filtered_feature_bc_matrix dir (matrix.mtx, features.tsv, barcodes.tsv)")
    parser.add_argument("-l", "--label_file", type=str, required=True,
                         help="cell-type label file with 'barcode_use'/'label' columns")
    parser.add_argument("-o", "--outdir", type=str, required=True,
                         help="where to write adata_RNA.h5ad / adata_ATAC.h5ad")
    parser.add_argument("--min_genes", type=int, default=200)
    parser.add_argument("--min_cells", type=int, default=3)
    args = parser.parse_args()

    run(args.matrix_dir, args.label_file, args.outdir, min_genes=args.min_genes, min_cells=args.min_cells)
