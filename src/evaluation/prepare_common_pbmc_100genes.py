"""Common gene/cell selection so Netmap and LINGER run on identical PBMC data.

QC filter -> normalize_total(1e4) -> log1p -> rank all genes by variability
-> drop PanglaoDB-ubiquitous genes -> keep the top n_genes by variability.

n_genes defaults to 1000 rather than ~100: LINGER's own LL_net.py hardcodes
batchsize=50 in its per-chromosome SHAP batching (TF_RE_LINGER_chr) and
crashes with UnboundLocalError if a chromosome has fewer than 50 associated
RE-TG rows, which happens reliably at a ~100-gene scale spread over 23
chromosomes; ~1000 genes keeps most chromosomes above that floor.

Takes the RNA/ATAC h5ad pair produced by prepare_linger_10x_multiome.py
(already QC'd, barcode-synced, with obs['sample']/obs['label']) and writes
three outputs sharing the exact same final cell set and the exact same final
gene set, each in the representation its consumer expects:
  - adata_netmap_100genes.h5ad: log-normalized (matches
    src/manuscript_blood/model_training/train_models_blood.py's expected
    input), for run_netmap_pbmc.py.
  - adata_RNA_matched.h5ad: raw counts, same cells + genes, for run_linger.py.
  - adata_ATAC_matched.h5ad: raw counts, same cells, full peak set (ATAC
    doesn't need gene-restriction -- LINGER's cis-regulatory step already
    restricts REs to those near the given target genes), for run_linger.py.

Restricting LINGER's adata_RNA to the same gene set (rather than the full
transcriptome) also makes it a closed TF/target system over that panel, same
as Netmap's autoencoder -- required for the runtime/result comparison
between the two to be apples-to-apples rather than LINGER scoring a much
larger gene universe.
"""
import os
import os.path as op
import sys
import time

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.append(op.abspath(op.join(op.dirname(__file__), '..', '..')))

from src.utils import write_config, sanitize_nullable_dtypes


def select_top_variable_genes(adata_rna, panglaodb_csv, n_genes=1000, min_genes=500, min_cells=50,
                               max_pct_mt=20, ubiquitousness_threshold=0.1):
    adata = adata_rna.copy()

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    adata = adata[adata.obs.pct_counts_mt < max_pct_mt].copy()

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=10000)
    adata.layers["count_norm"] = adata.X.copy()
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata)

    panglaodb = pd.read_csv(panglaodb_csv)
    high_ui = panglaodb[panglaodb["UI"] > ubiquitousness_threshold]["Official gene symbol"].values

    candidates = adata.var[~adata.var.index.isin(high_ui)]
    top_genes = candidates.sort_values("dispersions_norm", ascending=False).index[:n_genes]

    return adata[:, top_genes].copy()


def run(rna_h5ad, atac_h5ad, panglaodb_csv, outdir, n_genes=1000, min_genes=500, min_cells=50,
        max_pct_mt=20, ubiquitousness_threshold=0.1):
    os.makedirs(outdir, exist_ok=True)

    adata_RNA = sc.read_h5ad(rna_h5ad)
    adata_ATAC = sc.read_h5ad(atac_h5ad)

    start = time.monotonic()
    netmap_adata = select_top_variable_genes(
        adata_RNA, panglaodb_csv, n_genes=n_genes, min_genes=min_genes, min_cells=min_cells,
        max_pct_mt=max_pct_mt, ubiquitousness_threshold=ubiquitousness_threshold,
    )
    elapsed = time.monotonic() - start

    final_cells = netmap_adata.obs_names
    final_genes = netmap_adata.var_names

    rna_matched = sanitize_nullable_dtypes(adata_RNA[final_cells, final_genes].copy())
    atac_matched = sanitize_nullable_dtypes(adata_ATAC[final_cells].copy())

    sanitize_nullable_dtypes(netmap_adata).write_h5ad(op.join(outdir, "adata_netmap_100genes.h5ad"))
    rna_matched.write_h5ad(op.join(outdir, "adata_RNA_matched.h5ad"))
    atac_matched.write_h5ad(op.join(outdir, "adata_ATAC_matched.h5ad"))

    write_config(
        {
            "n_cells": int(len(final_cells)),
            "n_genes": int(len(final_genes)),
            "genes": list(final_genes),
            "selection_time": elapsed,
        },
        file=op.join(outdir, "preprocessing_summary.yaml"),
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Select a shared top-variable-gene, shared-cell RNA/ATAC subset for a fair Netmap/LINGER comparison."
    )
    parser.add_argument("-r", "--rna_h5ad", type=str, required=True,
                         help="adata_RNA.h5ad from prepare_linger_10x_multiome.py")
    parser.add_argument("-a", "--atac_h5ad", type=str, required=True,
                         help="adata_ATAC.h5ad from prepare_linger_10x_multiome.py")
    parser.add_argument("-p", "--panglaodb_csv", type=str, required=True,
                         help="PanglaoDB marker table CSV with 'UI'/'Official gene symbol' columns")
    parser.add_argument("-o", "--outdir", type=str, required=True)
    parser.add_argument("--n_genes", type=int, default=1000,
                         help="number of top-variability, non-ubiquitous genes to keep")
    parser.add_argument("--min_genes", type=int, default=500)
    parser.add_argument("--min_cells", type=int, default=50)
    parser.add_argument("--max_pct_mt", type=float, default=20)
    parser.add_argument("--ubiquitousness_threshold", type=float, default=0.1)
    args = parser.parse_args()

    run(args.rna_h5ad, args.atac_h5ad, args.panglaodb_csv, args.outdir, n_genes=args.n_genes,
        min_genes=args.min_genes, min_cells=args.min_cells, max_pct_mt=args.max_pct_mt,
        ubiquitousness_threshold=args.ubiquitousness_threshold)
