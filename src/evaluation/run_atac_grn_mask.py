"""Annotate a Netmap-predicted GRN with ATAC-supported TF->gene edges.

Runs netmap.atac's SnapATAC2-based pipeline (qc.run_qc_pipeline ->
network.run_atac_grn_pipeline) on a cells x peaks ATAC AnnData, then merges
the result into an existing Netmap `grn_adata` via
netmap.atac.network.add_atac_grn_support. Each `Source_Target` edge in
`grn_adata.var` gets a `regulator_{name}` flag -- True only if the TF has a
real CIS-BP motif match (SnapATAC2 `add_tf_binding`) inside a peak that is
both accessible in this dataset and within the TSS window of the target gene
-- plus an `{name}_cor_score` (TF/target activity correlation).

netmap.atac is not yet committed/pip-installable from the netmap repo, and
snapatac2 is not a netmap-evaluation dependency, so this script must be run
from netmap's `atac_env` pixi environment (sets PYTHONPATH to netmap/src),
not netmap-evaluation's own pixi env:
    cd <netmap_repo>/atac_env && pixi run python \
        <netmap-evaluation_repo>/src/evaluation/run_atac_grn_mask.py ...

The `grn_h5ad` input is a Netmap output h5ad (e.g. `{model_name}_grn.h5ad`
from src/manuscript_blood/model_training/train_models_blood.py), whose `.var`
has `source`/`target` columns.

Verified end to end against snapatac2==2.8.0 on real data (a subset of a real
PBMC ATAC dataset, real hg38 sequence, real CIS-BP motif scan) before being
wired up here; API calls match the installed snapatac2 version.
"""
import os
import os.path as op

import anndata as ad
import pandas as pd


def run(grn_h5ad, atac_h5ad, outdir, peak_col='gene_ids', genome='hg38',
        link_upstream=250000, link_downstream=250000, motif_pvalue=1e-5, min_abs_cor=0.1,
        qc_min_counts=1000, qc_n_features=200000, name='atac'):
    import snapatac2 as snap
    from netmap.atac import qc, network

    genome_obj = {'hg38': snap.genome.hg38, 'mm10': snap.genome.mm10}[genome]

    grn_adata = ad.read_h5ad(grn_h5ad)
    atac_adata = ad.read_h5ad(atac_h5ad)

    peak_adata = ad.AnnData(
        X=atac_adata.X,
        obs=pd.DataFrame(index=atac_adata.obs_names),
        var=pd.DataFrame(index=atac_adata.var[peak_col].astype(str).values),
    )

    qc.run_qc_pipeline(peak_adata, min_counts=qc_min_counts, n_features=qc_n_features)
    _pruned, edge_table = network.run_atac_grn_pipeline(
        peak_adata, genome_obj,
        link_upstream=link_upstream, link_downstream=link_downstream,
        motif_pvalue=motif_pvalue, min_abs_cor=min_abs_cor,
    )

    grn_adata = network.add_atac_grn_support(grn_adata, edge_table, name=name)

    os.makedirs(outdir, exist_ok=True)
    out_file = op.join(outdir, op.basename(grn_h5ad).replace('.h5ad', f'_{name}_masked.h5ad'))
    grn_adata.write(out_file)
    edge_table.to_csv(op.join(outdir, f'{name}_tf_gene_edges.tsv'), sep='\t', index=False)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Annotate a Netmap GRN with ATAC-supported TF-gene edges (SnapATAC2).")
    parser.add_argument("-r", "--grn_h5ad", type=str, required=True,
                         help="Netmap-predicted grn_adata h5ad (.var has 'source'/'target' columns)")
    parser.add_argument("-a", "--atac_h5ad", type=str, required=True,
                         help="ATAC AnnData h5ad; peak coordinates read from var[--peak_col]")
    parser.add_argument("-o", "--outdir", type=str, required=True)
    parser.add_argument("--peak_col", type=str, default="gene_ids",
                         help="var column holding 'chr:start-end' peak coordinates")
    parser.add_argument("--genome", type=str, default="hg38", choices=["hg38", "mm10"])
    parser.add_argument("--link_upstream", type=int, default=250000)
    parser.add_argument("--link_downstream", type=int, default=250000)
    parser.add_argument("--motif_pvalue", type=float, default=1e-5)
    parser.add_argument("--min_abs_cor", type=float, default=0.1)
    parser.add_argument("--qc_min_counts", type=int, default=1000,
                         help="min total peak counts per cell (SnapATAC2 default; lower for a small/subset ATAC input)")
    parser.add_argument("--qc_n_features", type=int, default=200000,
                         help="max peaks marked informative (SnapATAC2 default; lower for a small/subset ATAC input)")
    parser.add_argument("--name", type=str, default="atac")
    args = parser.parse_args()

    run(args.grn_h5ad, args.atac_h5ad, args.outdir, peak_col=args.peak_col, genome=args.genome,
        link_upstream=args.link_upstream, link_downstream=args.link_downstream,
        motif_pvalue=args.motif_pvalue, min_abs_cor=args.min_abs_cor,
        qc_min_counts=args.qc_min_counts, qc_n_features=args.qc_n_features, name=args.name)
