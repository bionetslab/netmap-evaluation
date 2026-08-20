"""Per-cell AUC benchmark: LINGER cell-level output vs. sparsified (+
ATAC-filtered) Netmap output, reusing the exact 3 ChIP-Atlas cases already
defined in src/manuscript_comparisons/compare_benchmark.ipynb.

That notebook's bm_trans() computes ONE pooled AUC per (method, case) from a
cell-type-pooled TG-by-TF score matrix -- one number summarizing an entire
cell type. This script keeps the same cases/ground truth (top-N ChIP-Atlas
target genes per TF) but computes one AUC PER INDIVIDUAL CELL instead, from
a cells x edges score matrix: LINGER's genuinely single-cell
cell_level_trans_regulatory.h5ad (see run_linger.py), or Netmap's sparsified/
ATAC-masked grn h5ad (see run_netmap_pbmc.py / run_atac_grn_mask.py).

CASES (verbatim from compare_benchmark.ipynb):
  case1: TF=IRF4,  cell_type='naive B cells',       ground truth 40215_gene_score_5fold.txt
  case2: TF=SPI1,  cell_type='classical monocytes',  ground truth 85986_gene_score_5fold.txt
  case3: TF=FOXP3, cell_type='naive CD4 T cells',    ground truth 44098_gene_score_5fold.txt

Ground truth is gene-level and cell-type-specific, not single-cell-resolved,
so the same top-N gold target-gene set applies uniformly to every cell of
the matching cell type -- only the per-cell predicted score vector varies.

--cell_type_col's values must match these case cell-type strings exactly;
if your data's cell-type labels (e.g. from PBMC_label.txt/Azimuth) use
different naming, pass --cell_type_map to override per case rather than
relabeling your data.

run_linger.py's cell-level h5ad carries no obs metadata (just barcodes), so
--linger_reference_h5ad (e.g. adata_RNA_matched.h5ad) is used to attach
--cell_type_col onto it by barcode before filtering by cell type.
"""
import os
import os.path as op

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as scs
from sklearn.metrics import roc_auc_score

CASES = {
    "case1": {"tf": "IRF4", "cell_type": "naive B cells", "ground_truth_file": "40215_gene_score_5fold.txt"},
    "case2": {"tf": "SPI1", "cell_type": "classical monocytes", "ground_truth_file": "85986_gene_score_5fold.txt"},
    "case3": {"tf": "FOXP3", "cell_type": "naive CD4 T cells", "ground_truth_file": "44098_gene_score_5fold.txt"},
}


def load_gold_targets(ground_truth_file, n_top=500):
    data = pd.read_csv(ground_truth_file, sep="\t", skiprows=5, header=0)
    ranked = data.groupby("symbol")["score"].max().sort_values(ascending=False)
    return set(ranked.index[:n_top])


def attach_cell_type(grn_adata, reference_h5ad, cell_type_col):
    reference = ad.read_h5ad(reference_h5ad)
    common = grn_adata.obs_names.intersection(reference.obs_names)
    grn_adata = grn_adata[common].copy()
    grn_adata.obs[cell_type_col] = reference.obs.loc[common, cell_type_col].values
    return grn_adata


def per_cell_auc(grn_adata, tf, gold_targets, cell_type_col, cell_type_value,
                  source_col="source", target_col="target", require_col=None):
    mask_cells = (grn_adata.obs[cell_type_col] == cell_type_value).to_numpy()
    if not mask_cells.any():
        return pd.DataFrame(columns=["cell", "auc", "n_targets", "n_positive"])

    var = grn_adata.var
    edge_mask = (var[source_col] == tf).to_numpy()
    if require_col is not None and require_col in var.columns:
        edge_mask &= var[require_col].to_numpy(dtype=bool)

    sub = grn_adata[mask_cells, edge_mask]
    targets = sub.var[target_col].to_numpy()
    labels = np.isin(targets, list(gold_targets)).astype(int)

    X = sub.X
    if scs.issparse(X):
        X = X.toarray()

    rows = []
    for i, cell in enumerate(sub.obs_names):
        auc = np.nan
        if 0 < labels.sum() < len(labels):
            try:
                auc = roc_auc_score(labels, X[i, :])
            except ValueError:
                auc = np.nan
        rows.append({"cell": cell, "auc": auc, "n_targets": len(labels), "n_positive": int(labels.sum())})
    return pd.DataFrame(rows)


def run(ground_truth_dir, outdir, linger_trans_h5ad=None, linger_reference_h5ad=None,
        netmap_h5ad=None, netmap_atac_h5ad=None, cell_type_col="label", cell_type_map=None,
        n_top_genes=500, atac_regulator_col="regulator_atac"):
    os.makedirs(outdir, exist_ok=True)
    cell_type_map = cell_type_map or {}

    methods = {}
    if linger_trans_h5ad is not None:
        linger_adata = ad.read_h5ad(linger_trans_h5ad)
        if linger_reference_h5ad is not None:
            linger_adata = attach_cell_type(linger_adata, linger_reference_h5ad, cell_type_col)
        methods["LINGER"] = dict(adata=linger_adata, source_col="TF", target_col="TG", require_col=None)
    if netmap_h5ad is not None:
        methods["Netmap_sparsified"] = dict(adata=ad.read_h5ad(netmap_h5ad), source_col="source",
                                             target_col="target", require_col=None)
    if netmap_atac_h5ad is not None:
        methods["Netmap_sparsified_ATAC"] = dict(adata=ad.read_h5ad(netmap_atac_h5ad), source_col="source",
                                                  target_col="target", require_col=atac_regulator_col)

    all_results = []
    for case_name, case in CASES.items():
        gold_targets = load_gold_targets(op.join(ground_truth_dir, case["ground_truth_file"]), n_top=n_top_genes)
        cell_type_value = cell_type_map.get(case_name, case["cell_type"])

        for method_name, spec in methods.items():
            result = per_cell_auc(spec["adata"], case["tf"], gold_targets, cell_type_col, cell_type_value,
                                   source_col=spec["source_col"], target_col=spec["target_col"],
                                   require_col=spec["require_col"])
            result["case"] = case_name
            result["method"] = method_name
            result["tf"] = case["tf"]
            result["cell_type"] = cell_type_value
            all_results.append(result)

    all_results = pd.concat(all_results, ignore_index=True) if all_results else pd.DataFrame()
    all_results.to_csv(op.join(outdir, "per_cell_auc.tsv"), sep="\t", index=False)

    summary = (all_results.groupby(["case", "method"])["auc"]
               .agg(["mean", "median", "std", "count"]).reset_index())
    summary.to_csv(op.join(outdir, "per_cell_auc_summary.tsv"), sep="\t", index=False)
    return all_results, summary


if __name__ == "__main__":
    import argparse

    def _parse_map(value):
        mapping = {}
        for item in value.split(","):
            if not item:
                continue
            k, v = item.split("=", 1)
            mapping[k] = v
        return mapping

    parser = argparse.ArgumentParser(
        description="Per-cell AUC benchmark of LINGER vs. sparsified(+ATAC) Netmap GRNs, reusing "
                    "compare_benchmark.ipynb's 3 ChIP-Atlas cases."
    )
    parser.add_argument("-g", "--ground_truth_dir", type=str, required=True,
                         help="dir containing 40215_gene_score_5fold.txt / 85986_.../ 44098_...")
    parser.add_argument("-o", "--outdir", type=str, required=True)
    parser.add_argument("--linger_trans_h5ad", type=str, default=None,
                         help="cell_level_trans_regulatory.h5ad from run_linger.py")
    parser.add_argument("--linger_reference_h5ad", type=str, default=None,
                         help="adata_RNA_matched.h5ad (or similar), used only to attach --cell_type_col "
                              "onto the LINGER h5ad by barcode, since run_linger.py's cell-level output "
                              "carries no obs metadata")
    parser.add_argument("--netmap_h5ad", type=str, default=None,
                         help="{model_name}_sparsified.h5ad from run_netmap_pbmc.py")
    parser.add_argument("--netmap_atac_h5ad", type=str, default=None,
                         help="{model_name}_sparsified_{name}_masked.h5ad from run_atac_grn_mask.py "
                              "(run with --grn_h5ad pointed at the sparsified h5ad)")
    parser.add_argument("--cell_type_col", type=str, default="label")
    parser.add_argument("--cell_type_map", type=_parse_map, default={},
                         help="override a case's cell-type string, e.g. "
                              "'case1=Naive B,case2=CD14 Mono' if your labels don't match "
                              "the ChIP-Atlas case names verbatim")
    parser.add_argument("--n_top_genes", type=int, default=500)
    parser.add_argument("--atac_regulator_col", type=str, default="regulator_atac")
    args = parser.parse_args()

    run(args.ground_truth_dir, args.outdir, linger_trans_h5ad=args.linger_trans_h5ad,
        linger_reference_h5ad=args.linger_reference_h5ad, netmap_h5ad=args.netmap_h5ad,
        netmap_atac_h5ad=args.netmap_atac_h5ad, cell_type_col=args.cell_type_col,
        cell_type_map=args.cell_type_map, n_top_genes=args.n_top_genes,
        atac_regulator_col=args.atac_regulator_col)
