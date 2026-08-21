"""Single combined preprocessing step: raw 10x multiome -> the shared RNA/ATAC/
gene-set inputs both run_netmap_pbmc.py and run_linger.py need.

Replaces the previous two-script pipeline (prepare_linger_10x_multiome.py +
prepare_common_pbmc_100genes.py). That split ran two independent scanpy QC
passes back to back (LINGER's own hardcoded pct_counts_mt<5 inside
LingerGRN.preprocess.get_adata(), then a second, differently-thresholded QC
pass at the gene-selection stage) and required LINGER's package just to
ingest the raw matrix. Here there is exactly one QC pass, one preprocessing
script, and no LINGER import at all.

get_adata()'s RNA/ATAC-split, barcode/label attachment, and mt-filter logic
is replicated by hand below (verified against the installed LingerGRN==1.110
source) with the mito threshold exposed as --max_pct_mt (default 10) instead
of LINGER's hardcoded 5.

Steps:
  1. load_raw_multiome(): split the combined 10x matrix into RNA/ATAC by
     feature type, attach obs['barcode']/obs['sample']/obs['label'],
     restrict to labeled barcodes.
  2. qc_filter(): pct_counts_mt < max_pct_mt on RNA, then min_genes/min_cells
     filtering on both modalities, then RNA/ATAC barcode intersection.
  3. select_common_genes(): normalize_total -> log1p -> rank by variability
     -> drop PanglaoDB-ubiquitous genes -> keep top n_genes -> force-include
     every TF from tf_names_file present in the data. Top-variability
     ranking alone can select a gene panel with few or no real TFs (they
     tend to be moderately expressed regulators, not top-dispersion genes),
     and LINGER's preprocess.TF_expression() raises IndexError downstream if
     its TF/gene-panel intersection is empty -- TFs are added unconditionally,
     bypassing the ubiquitousness/variability filters, since LINGER needs a
     real TF pool regardless of whether they'd otherwise qualify.

Writes, all sharing the exact same final cell set:
  - adata_netmap.h5ad: log-normalized, final gene set, with a 'counts' layer
    of raw counts (run_netmap_pbmc.py's masking step needs it), for
    run_netmap_pbmc.py.
  - adata_RNA_matched.h5ad: raw counts, same final gene set, with
    obs['sample']/obs['label'], for run_linger.py.
  - adata_ATAC_matched.h5ad: raw counts, same final cells, full peak set
    (LINGER's cis-regulatory step restricts REs to those near the given
    target genes itself), for run_linger.py.
  - preprocessing_summary.yaml: n_cells/n_genes/n_tfs_included/genes/timing.

No LINGER import means this can run in the main env (not linger_env).
Outputs are still passed through sanitize_nullable_dtypes() before writing
so LINGER's pinned anndata==0.9.2 can read them regardless.
"""
import os
import os.path as op
import sys
import time

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy.io

sys.path.append(op.abspath(op.join(op.dirname(__file__), '..', '..')))

from src.utils import write_config, sanitize_nullable_dtypes


def load_raw_multiome(matrix_dir, label_file):
    matrix = scipy.io.mmread(op.join(matrix_dir, 'matrix.mtx'))
    matrix.data = matrix.data.astype(np.float32)
    features = pd.read_csv(op.join(matrix_dir, 'features.tsv'), sep='\t', header=None)
    barcodes = pd.read_csv(op.join(matrix_dir, 'barcodes.tsv'), sep='\t', header=None)
    label = pd.read_csv(label_file, sep='\t', header=0)
    label.index = label['barcode_use']

    adata = ad.AnnData(X=matrix.T.tocsr())
    adata.var['gene_ids'] = features[1].values
    adata.var_names = features[1].values
    adata.obs['barcode'] = barcodes[0].values
    adata.obs_names = barcodes[0].values
    if len(barcodes[0].values[0].split('-')) == 2:
        adata.obs['sample'] = [b.split('-')[1] for b in barcodes[0].values]
    else:
        adata.obs['sample'] = '1'

    adata_RNA = adata[:, (features[2] == 'Gene Expression').values].copy()
    adata_ATAC = adata[:, (features[2] == 'Peaks').values].copy()

    idx = adata_RNA.obs['barcode'].isin(label['barcode_use'].values)
    adata_RNA = adata_RNA[idx].copy()
    adata_ATAC = adata_ATAC[idx].copy()
    adata_RNA.obs['label'] = label.loc[adata_RNA.obs['barcode']]['label'].values
    adata_ATAC.obs['label'] = label.loc[adata_ATAC.obs['barcode']]['label'].values

    return adata_RNA, adata_ATAC


def qc_filter(adata_RNA, adata_ATAC, max_pct_mt=10, min_genes=500, min_cells=50):
    adata_RNA.var['mt'] = adata_RNA.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata_RNA, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    adata_RNA = adata_RNA[adata_RNA.obs.pct_counts_mt < max_pct_mt, :].copy()
    adata_RNA.var.index = adata_RNA.var['gene_ids'].values
    adata_RNA.var_names_make_unique()
    adata_RNA.var['gene_ids'] = adata_RNA.var.index

    sc.pp.filter_cells(adata_RNA, min_genes=min_genes)
    sc.pp.filter_genes(adata_RNA, min_cells=min_cells)
    sc.pp.filter_cells(adata_ATAC, min_genes=min_genes)
    sc.pp.filter_genes(adata_ATAC, min_cells=min_cells)

    selected_barcode = list(set(adata_RNA.obs['barcode']) & set(adata_ATAC.obs['barcode']))
    rna_idx = pd.DataFrame(range(adata_RNA.shape[0]), index=adata_RNA.obs['barcode'].values)
    adata_RNA = adata_RNA[rna_idx.loc[selected_barcode][0]].copy()
    atac_idx = pd.DataFrame(range(adata_ATAC.shape[0]), index=adata_ATAC.obs['barcode'].values)
    adata_ATAC = adata_ATAC[atac_idx.loc[selected_barcode][0]].copy()
    return adata_RNA, adata_ATAC


def select_common_genes(adata_RNA, panglaodb_csv, tf_names_file, n_genes=1000, ubiquitousness_threshold=0.1):
    adata = adata_RNA.copy()
    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=10000)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata)

    panglaodb = pd.read_csv(panglaodb_csv)
    high_ui = panglaodb[panglaodb['UI'] > ubiquitousness_threshold]['Official gene symbol'].values
    candidates = adata.var[~adata.var.index.isin(high_ui)]
    top_genes = candidates.sort_values('dispersions_norm', ascending=False).index[:n_genes]

    tf_names = pd.read_csv(tf_names_file, header=None)[0].values
    tfs_present = adata.var_names[adata.var_names.isin(tf_names)]

    final_genes = pd.Index(top_genes).union(tfs_present)
    return adata[:, final_genes].copy(), tfs_present


def run(matrix_dir, label_file, panglaodb_csv, tf_names_file, outdir, max_pct_mt=10, min_genes=500,
        min_cells=50, n_genes=1000, ubiquitousness_threshold=0.1):
    os.makedirs(outdir, exist_ok=True)

    start = time.monotonic()
    adata_RNA, adata_ATAC = load_raw_multiome(matrix_dir, label_file)
    adata_RNA, adata_ATAC = qc_filter(adata_RNA, adata_ATAC, max_pct_mt=max_pct_mt,
                                       min_genes=min_genes, min_cells=min_cells)
    netmap_adata, tfs_present = select_common_genes(
        adata_RNA, panglaodb_csv, tf_names_file, n_genes=n_genes,
        ubiquitousness_threshold=ubiquitousness_threshold,
    )
    elapsed = time.monotonic() - start

    final_cells = netmap_adata.obs_names
    final_genes = netmap_adata.var_names

    rna_matched = adata_RNA[final_cells, final_genes].copy()
    atac_matched = adata_ATAC[final_cells].copy()

    sanitize_nullable_dtypes(netmap_adata).write_h5ad(op.join(outdir, 'adata_netmap.h5ad'))
    sanitize_nullable_dtypes(rna_matched).write_h5ad(op.join(outdir, 'adata_RNA_matched.h5ad'))
    sanitize_nullable_dtypes(atac_matched).write_h5ad(op.join(outdir, 'adata_ATAC_matched.h5ad'))

    write_config(
        {
            'n_cells': int(len(final_cells)),
            'n_genes': int(len(final_genes)),
            'n_tfs_included': int(len(tfs_present)),
            'genes': list(final_genes),
            'max_pct_mt': max_pct_mt,
            'preprocessing_time': elapsed,
        },
        file=op.join(outdir, 'preprocessing_summary.yaml'),
    )


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Combined preprocessing: raw 10x multiome -> shared RNA/ATAC/gene-set inputs for both '
                     'Netmap (run_netmap_pbmc.py) and LINGER (run_linger.py).'
    )
    parser.add_argument('-m', '--matrix_dir', type=str, required=True,
                         help='10x filtered_feature_bc_matrix dir (matrix.mtx, features.tsv, barcodes.tsv)')
    parser.add_argument('-l', '--label_file', type=str, required=True,
                         help="cell-type label file with 'barcode_use'/'label' columns")
    parser.add_argument('-p', '--panglaodb_csv', type=str, required=True,
                         help="PanglaoDB marker table CSV with 'UI'/'Official gene symbol' columns")
    parser.add_argument('-t', '--tf_names_file', type=str, required=True,
                         help='headerless TF-name list; every TF present in the data is force-included '
                              'regardless of variability/ubiquitousness (LINGER needs a real TF pool)')
    parser.add_argument('-o', '--outdir', type=str, required=True)
    parser.add_argument('--max_pct_mt', type=float, default=10,
                         help="max %% mitochondrial reads per cell (LINGER's own get_adata() hardcodes 5)")
    parser.add_argument('--min_genes', type=int, default=500)
    parser.add_argument('--min_cells', type=int, default=50)
    parser.add_argument('--n_genes', type=int, default=1000,
                         help='number of top-variability, non-ubiquitous genes to keep')
    parser.add_argument('--ubiquitousness_threshold', type=float, default=0.1)
    args = parser.parse_args()

    run(args.matrix_dir, args.label_file, args.panglaodb_csv, args.tf_names_file, args.outdir,
        max_pct_mt=args.max_pct_mt, min_genes=args.min_genes, min_cells=args.min_cells,
        n_genes=args.n_genes, ubiquitousness_threshold=args.ubiquitousness_threshold)
