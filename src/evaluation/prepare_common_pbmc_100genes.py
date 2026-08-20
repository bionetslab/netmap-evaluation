"""Common gene/cell selection so Netmap and LINGER run on identical PBMC data.

Replicates the ~100-gene selection used for the blood manuscript
(object_with_markers in
src/manuscript_blood/data_preprocessing/preprocess_blood_rhap_rep2_include_markers.ipynb):
QC filter -> normalize_total(1e4) -> log1p -> HVG(n_top_genes) -> keep genes
that are (HVG AND dropout<95% AND not PanglaoDB-ubiquitous) OR in the curated
immune MARKER_DICT below (same list as the blood pipeline; PBMC is the same
immune/blood cell population so it is reused verbatim rather than recreated).

Takes the RNA/ATAC h5ad pair produced by prepare_linger_10x_multiome.py
(already QC'd, barcode-synced, with obs['sample']/obs['label']) and writes
three outputs sharing the exact same final cell set and the exact same final
~100-gene set, each in the representation its consumer expects:
  - adata_netmap_100genes.h5ad: log-normalized (matches
    src/manuscript_blood/model_training/train_models_blood.py's expected
    input), for run_netmap_pbmc.py.
  - adata_RNA_matched.h5ad: raw counts, same cells + genes, for run_linger.py.
  - adata_ATAC_matched.h5ad: raw counts, same cells, full peak set (ATAC
    doesn't need gene-restriction -- LINGER's cis-regulatory step already
    restricts REs to those near the given target genes), for run_linger.py.

Restricting LINGER's adata_RNA to the same ~100 genes (rather than the full
transcriptome) also makes it a closed TF/target system over that gene panel,
same as Netmap's autoencoder -- required for the runtime/result comparison
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

MARKER_DICT = {
    # --- B cells and Plasma cells ---
    "B cell": ["CD19", "CD79A", "CD79B", "CD21", "IGHM", "IGHD", "IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHE", "IGKC", "IGLC2", "IGLC3", "MS4A1", "PAX5", "HLA-DR", "HLA-DQ"],
    "Transitional B cell": ["MME", "CD9", "CD38", "PAX5", "FCER2"],
    "Naive B cell": ["IL4R", "FCER2", "IGHM", "IGHD"],
    "Memory B cell": ["CD27", "AIM2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHE"],
    "Effector B cell": ["ITGAX", "FCRL4", "FCRL5", "TBX21", "ZEB2", "PDCD1"],
    "Plasma cell": ["CD19", "PRDM1", "XBP1", "MZB1", "SLAMF7", "CD27", "CD38", "IGHM", "IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGKC", "IGLC2", "IGLC3"],
    "Core naive B cell": ["IL4R", "FCER2", "IGHM", "IGHD"],
    "ISG+ naive B cell": ["STAT3", "STAT1", "IFI44L", "ISG15"],
    "Early memory B cell": ["CD27", "AIM2", "IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4"],
    "Core memory B cell": ["CD27", "AIM2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHE"],
    "Type 2 polarized memory B cell": ["IL4R", "FCER2", "COCH", "IGHG1", "IGHG4", "IGHE"],
    "CD95 memory B cell": ["FAS", "AIM2", "IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4"],
    "Activated memory B cell": ["FOS", "CD69", "JUN", "MCL1", "MYC"],
    "CD27+ effector B cell": ["CD27"],
    "CD27- effector B cell": ["ITGAX", "TBX21", "ZEB2"],

    # --- T cell (shared Level 1) ---
    "T cell": ["TRAC", "TRDC", "CD3D", "CD3E", "CD3G"],

    # --- CD4 T, DN T, and Treg cells ---
    "Naive CD4 T cell": ["CD27", "CCR7", "SELL", "TCF7", "LEF1"],
    "Memory CD4 T cell": ["ITGB1"],
    "Treg": ["FOXP3", "IL2RA", "IKZF2", "RTKN2"],
    "DN T cell": ["TRAC"],
    "Proliferating T cell": ["MKI67"],
    "Core naive CD4 T cell": ["CD27", "CCR7", "SELL", "TCF7", "LEF1"],
    "SOX4+ naive CD4 T cell": ["SOX4"],
    "ISG+ naive CD4 T cell": ["MX1", "IFI44"],
    "CM CD4 T cell": ["CCR7", "SELL", "LEF1"],
    "GZMB- CD27- EM CD4 T cell": [],
    "GZMB- CD27+ EM CD4 T cell": ["CD27", "GZMK"],
    "KLRF1- GZMB+ CD27- memory CD4 T cell": ["GZMB", "CCL5"],
    "ISG+ memory CD4 T cell": ["MX1", "IFI44"],
    "Naive CD4 Treg": ["CD27", "CCR7", "SELL", "TCF7", "LEF1"],
    "Memory CD4 Treg": ["ITGB1"],
    "KLRB1+ memory CD4 Treg": ["KLRB1"],
    "GZMK+ memory CD4 Treg": ["GZMK"],
    "Memory CD8 Treg": ["CD8A"],
    "KLRB1+ memory CD8 Treg": ["KLRB1"],

    # --- CD8 T, gdT, and MAIT cells ---
    "Naive CD8 T cell": ["CD27", "CCR7", "SELL", "TCF7", "LEF1"],
    "Memory CD8 T cell": ["ITGB1", "GZMA", "GZMB", "GZMK", "TRAC"],
    "CD8aa": ["CD8A", "KLRC2", "IKZF2", "IL21R"],
    "MAIT": ["SLC4A10", "KLRB1"],
    "gdT": ["TRDC", "TRGC1", "TRGC2"],
    "Core naive CD8 T cell": ["CD27", "CCR7", "SELL", "TCF7", "LEF1"],
    "SOX4+ naive CD8 T cell": ["SOX4"],
    "ISG+ naive CD8 T cell": ["MX1", "IFI44"],
    "CM CD8 T cell": ["CCR7", "SELL", "LEF1"],
    "GZMK- CD27+ EM CD8 T cell": ["CD27"],
    "GZMK+ CD27+ EM CD8 T cell": ["CD27", "GZMK"],
    "KLRF1- GZMB+ CD27- EM CD8 T cell": ["GZMB", "CCL5"],
    "KLRF1+ GZMB+ CD27- EM CD8 T cell": ["KLRF1", "GZMB", "CCL5"],
    "ISG+ memory CD8 T cell": ["MX1", "IFI44"],
    "Naive Vd1 gdT": ["CCR7", "SELL", "LEF1"],
    "SOX4+ Vd1 gdT": ["SOX4"],
    "KLRF1- effector Vd1 gdT": ["TRDC", "TRDV1"],
    "KLRF1+ effector Vd1 gdT": ["TRDC", "TRDV1", "KLRF1"],
    "GZMB+ Vd2 gdT": ["TRDC", "TRDV2", "GZMB"],
    "GZMK+ Vd2 gdT": ["TRDC", "TRDV2", "GZMK"],
    "CD8 MAIT": ["CD8A"],
    "CD4 MAIT": ["CD4"],
    "ISG+ MAIT": ["MX1", "IFI44"],

    # --- NK cells and ILCs ---
    "NK cell": ["CD3E"],
    "ILC": [],
    "CD56bright NK cell": ["NCAM1"],
    "CD56dim NK cell": ["FCGR3A"],
    "Proliferating NK cell": ["MKI67"],
    "GZMK- CD56dim NK cell": [],
    "GZMK+ CD56dim NK cell": ["GZMK"],
    "Adaptive NK cell": ["KLRC2", "FCGR3A"],
    "ISG+ CD56dim NK cell": ["ISG15", "MX1", "MX2"],

    # --- Monocytes ---
    "Monocyte": ["FCN1", "CTSS"],
    "CD14 monocyte": ["CD14", "VCAN", "S100A8", "S100A9"],
    "CD16 monocyte": ["FCGR3A", "CDKN1C", "LST1"],
    "Intermediate monocyte": ["CD14", "FCGR3A", "HLA-DPA1", "HLA-DOA", "HLA-DRA", "CD74"],
    "Core CD14 monocyte": [],
    "IL1B+ CD14 monocyte": ["IL1B", "CCL3", "CXCL8"],
    "ISG+ CD14 monocyte": ["MX1", "IFI44L", "IFI6"],
    "Core CD16 monocyte": [],
    "C1Q+ CD16 monocyte": ["C1QA", "C1QB"],
    "ISG+ CD16 monocyte": ["MX1", "IFI44L", "IFI6"],

    # --- Dendritic cells ---
    "DC": ["CST3", "FLT3", "HLA-DPA1", "HLA-DRA", "CD74"],
    "ASDC": ["AXL", "SIGLEC6", "HAMP"],
    "cDC1": ["CLEC9A", "XCR1", "IDO1", "C1orf54"],
    "cDC2": ["CD1C", "FCN1", "PILRA"],
    "pDC": ["PTCRA", "SMIM5", "LAMP5", "IL3RA", "JCHAIN"],
    "CD14+ cDC2": ["CST3", "CD74", "HLA-DRA", "HLA-DPA1", "CD14", "S100A8", "S100A9", "VCAN", "CD163"],
    "HLA-DRhi cDC2": ["HLA-DPA1", "HLA-DRA", "CD1C"],
    "ISG+ cDC2": ["MX1", "IFI44L", "IFI6"],

    # --- Progenitors, Erythrocytes, and Platelets ---
    "Erythrocyte": ["HBA1", "HBA2", "HBB"],
    "Progenitor cell": ["SMIM24", "CD34"],
    "Platelet": ["PPBP", "TUBB1"],
    "CLP cell": [],
    "CMP cell": [],
    "BaEoMaP cell": [],
}
MARKER_FLAT = sorted({gene for genes in MARKER_DICT.values() for gene in genes})


def select_marker_genes(adata_rna, panglaodb_csv, n_top_genes=2500, min_genes=500, min_cells=50,
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
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)

    panglaodb = pd.read_csv(panglaodb_csv)
    high_ui = panglaodb[panglaodb["UI"] > ubiquitousness_threshold]["Official gene symbol"].values

    flg = ((adata.var.pct_dropout_by_counts < 95) & (adata.var.highly_variable)
           & (~adata.var.index.isin(high_ui))) | (adata.var.index.isin(MARKER_FLAT))

    return adata[:, flg].copy()


def run(rna_h5ad, atac_h5ad, panglaodb_csv, outdir, n_top_genes=2500, min_genes=500, min_cells=50,
        max_pct_mt=20, ubiquitousness_threshold=0.1):
    os.makedirs(outdir, exist_ok=True)

    adata_RNA = sc.read_h5ad(rna_h5ad)
    adata_ATAC = sc.read_h5ad(atac_h5ad)

    start = time.monotonic()
    netmap_adata = select_marker_genes(
        adata_RNA, panglaodb_csv, n_top_genes=n_top_genes, min_genes=min_genes, min_cells=min_cells,
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
        description="Select a shared ~100-gene, shared-cell RNA/ATAC subset for a fair Netmap/LINGER comparison."
    )
    parser.add_argument("-r", "--rna_h5ad", type=str, required=True,
                         help="adata_RNA.h5ad from prepare_linger_10x_multiome.py")
    parser.add_argument("-a", "--atac_h5ad", type=str, required=True,
                         help="adata_ATAC.h5ad from prepare_linger_10x_multiome.py")
    parser.add_argument("-p", "--panglaodb_csv", type=str, required=True,
                         help="PanglaoDB marker table CSV with 'UI'/'Official gene symbol' columns")
    parser.add_argument("-o", "--outdir", type=str, required=True)
    parser.add_argument("--n_top_genes", type=int, default=2500)
    parser.add_argument("--min_genes", type=int, default=500)
    parser.add_argument("--min_cells", type=int, default=50)
    parser.add_argument("--max_pct_mt", type=float, default=20)
    parser.add_argument("--ubiquitousness_threshold", type=float, default=0.1)
    args = parser.parse_args()

    run(args.rna_h5ad, args.atac_h5ad, args.panglaodb_csv, args.outdir, n_top_genes=args.n_top_genes,
        min_genes=args.min_genes, min_cells=args.min_cells, max_pct_mt=args.max_pct_mt,
        ubiquitousness_threshold=args.ubiquitousness_threshold)
