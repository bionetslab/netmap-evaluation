"""Run LINGER (https://github.com/Durenlab/LINGER) to produce cell-level GRNs.

LINGER's public API (LingerGRN.LL_net) only ships "cell population" and
"cell type specific" GRN functions. Those functions select cells via a
boolean mask on `adata_RNA.obs['label'] == celltype` and mean-pool RNA/ATAC
values across the matching cells before scoring edges. A group of exactly
one cell degenerates to that cell's own value, so passing a per-cell label
(e.g. the barcode) instead of a cell-type annotation and calling with
celltype='all' makes these functions emit genuinely single-cell-resolution
networks. This is not documented as a feature, but it follows directly from
the grouping math in LL_net.py (verified against the source) and requires no
modification of LINGER's internals.

Two cell-level network types are supported:
  - cis:   RE -> TG cis-regulatory network. Cheap: needs only pseudobulking
           + preprocess() (peak overlap with the GRNdir reference), no
           neural-network training, since the RE-TG prior and distance
           weights come straight from the GRNdir bulk reference.
  - trans: TF -> TG trans-regulatory network (the classic "GRN"). Expensive:
           needs the population-level LINGER model trained first, and its
           per-cell TF-RE binding matrix is a dense RE x TF table written to
           disk per cell before being multiplied into the trans matrix.

Both require paired, barcode-aligned scRNA + scATAC AnnData input with
`obs['sample']` (for pseudobulking) and `obs['label']` (cell type, only used
by the population/cell-type-specific steps, not by the cell-level ones). To
build that pair from a raw 10x multiome download, run
download_pbmc_multiome_for_linger.sh then prepare_linger_10x_multiome.py
first. To run LINGER on the exact same cells/genes as Netmap (for a runtime
comparison), run prepare_common_pbmc_100genes.py on top of that pair first
and pass its adata_RNA_matched.h5ad/adata_ATAC_matched.h5ad here instead of
prepare_linger_10x_multiome.py's raw output.

Cell-level extraction loops LINGER's per-group pipeline once per requested
cell, each iteration re-scanning all chromosomes, so it does not scale to a
full single-cell dataset. Use --n_cells/--cells_file to restrict the run to
a manageable subset (a few hundred cells for `cis`, far fewer for `trans`).
"""
import os
import os.path as op
import sys
import time

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.sparse import csr_matrix

sys.path.append(op.abspath(op.join(op.dirname(__file__), '..', '..')))

from src.utils import write_config


def _pseudobulk_and_preprocess(adata_RNA, adata_ATAC, grn_dir, genome, method, outdir):
    from LingerGRN.pseudo_bulk import pseudo_bulk
    from LingerGRN.preprocess import preprocess

    samplelist = list(adata_ATAC.obs['sample'].unique())
    n_samples = adata_RNA.obs['sample'].nunique()
    singlepseudobulk = n_samples * n_samples > 100

    TG_pseudobulk = pd.DataFrame([])
    RE_pseudobulk = pd.DataFrame([])
    for sample in samplelist:
        rna_s = adata_RNA[adata_RNA.obs['sample'] == sample]
        atac_s = adata_ATAC[adata_ATAC.obs['sample'] == sample]
        tg, re = pseudo_bulk(rna_s, atac_s, singlepseudobulk)
        TG_pseudobulk = pd.concat([TG_pseudobulk, tg], axis=1)
        RE_pseudobulk = pd.concat([RE_pseudobulk, re], axis=1)

    RE_pseudobulk[RE_pseudobulk > 100] = 100
    TG_pseudobulk = TG_pseudobulk.fillna(0)
    RE_pseudobulk = RE_pseudobulk.fillna(0)
    pd.DataFrame(adata_ATAC.var['gene_ids']).to_csv(op.join(outdir, 'Peaks.txt'), header=None, index=None)

    preprocess(TG_pseudobulk, RE_pseudobulk, grn_dir, genome, method, outdir)


def _train_population_model(adata_RNA, adata_ATAC, grn_dir, genome, method, activation, species, outdir):
    import LingerGRN.LINGER_tr as LINGER_tr
    from LingerGRN.LL_net import TF_RE_binding, cis_reg

    LINGER_tr.training(grn_dir, method, outdir, activation, species)
    TF_RE_binding(grn_dir, adata_RNA, adata_ATAC, genome, method, outdir)
    cis_reg(grn_dir, adata_RNA, adata_ATAC, genome, method, outdir)


def _select_cells(adata_RNA, adata_ATAC, n_cells, cells_file):
    if cells_file is not None:
        wanted = set(pd.read_csv(cells_file, header=None)[0].astype(str))
        keep = adata_RNA.obs_names.isin(wanted)
    elif n_cells is not None and n_cells < adata_RNA.n_obs:
        rng = np.random.default_rng(0)
        idx = rng.choice(adata_RNA.n_obs, size=n_cells, replace=False)
        keep = np.zeros(adata_RNA.n_obs, dtype=bool)
        keep[idx] = True
    else:
        keep = np.ones(adata_RNA.n_obs, dtype=bool)
    return adata_RNA[keep].copy(), adata_ATAC[keep].copy()


def _run_cell_level(adata_RNA, adata_ATAC, grn_dir, genome, method, outdir, networks):
    from LingerGRN.LL_net import (
        cell_type_specific_TF_RE_binding,
        cell_type_specific_cis_reg,
        cell_type_specific_trans_reg,
    )

    barcodes = adata_RNA.obs_names.astype(str).tolist()
    adata_RNA = adata_RNA.copy()
    adata_RNA.obs['label'] = barcodes

    need_binding = 'trans' in networks
    need_cis = 'cis' in networks or 'trans' in networks

    if need_binding:
        cell_type_specific_TF_RE_binding(grn_dir, adata_RNA, adata_ATAC, genome, 'all', outdir, method)
    if need_cis:
        cell_type_specific_cis_reg(grn_dir, adata_RNA, adata_ATAC, genome, 'all', outdir, method)
    if 'trans' in networks:
        cell_type_specific_trans_reg(grn_dir, adata_RNA, 'all', outdir)

    return barcodes


def _assemble_cis_adata(outdir, barcodes, keep_intermediate):
    records = []
    for bc in barcodes:
        f = op.join(outdir, f'cell_type_specific_cis_regulatory_{bc}.txt')
        df = pd.read_csv(f, sep='\t', header=None, names=['RE', 'TG', 'score'])
        df['cell'] = bc
        records.append(df)

    all_edges = pd.concat(records, ignore_index=True)
    all_edges['edge'] = all_edges['RE'] + '|' + all_edges['TG']
    edge_ids = pd.Index(all_edges['edge'].unique())
    cell_ids = pd.Index(barcodes)

    row = pd.Categorical(all_edges['cell'], categories=cell_ids).codes
    col = pd.Categorical(all_edges['edge'], categories=edge_ids).codes
    mat = csr_matrix((all_edges['score'].values, (row, col)), shape=(len(cell_ids), len(edge_ids)))

    var = all_edges.drop_duplicates('edge').set_index('edge').loc[edge_ids, ['RE', 'TG']]
    adata = ad.AnnData(X=mat, obs=pd.DataFrame(index=cell_ids), var=var)

    if not keep_intermediate:
        for bc in barcodes:
            os.remove(op.join(outdir, f'cell_type_specific_cis_regulatory_{bc}.txt'))
    return adata


def _assemble_trans_adata(outdir, barcodes, keep_intermediate):
    records = []
    for bc in barcodes:
        f = op.join(outdir, f'cell_type_specific_trans_regulatory_{bc}.txt')
        dense = pd.read_csv(f, sep='\t', index_col=0)
        dense.index.name = 'TG'
        long = dense.reset_index().melt(id_vars='TG', var_name='TF', value_name='score')
        long = long[long['score'] != 0]
        long['cell'] = bc
        records.append(long)

    all_edges = pd.concat(records, ignore_index=True)
    all_edges['edge'] = all_edges['TF'] + '|' + all_edges['TG']
    edge_ids = pd.Index(all_edges['edge'].unique())
    cell_ids = pd.Index(barcodes)

    row = pd.Categorical(all_edges['cell'], categories=cell_ids).codes
    col = pd.Categorical(all_edges['edge'], categories=edge_ids).codes
    mat = csr_matrix((all_edges['score'].values, (row, col)), shape=(len(cell_ids), len(edge_ids)))

    var = all_edges.drop_duplicates('edge').set_index('edge').loc[edge_ids, ['TF', 'TG']]
    adata = ad.AnnData(X=mat, obs=pd.DataFrame(index=cell_ids), var=var)

    if not keep_intermediate:
        for bc in barcodes:
            os.remove(op.join(outdir, f'cell_type_specific_trans_regulatory_{bc}.txt'))
            binding_f = op.join(outdir, f'cell_type_specific_TF_RE_binding_{bc}.txt')
            if op.isfile(binding_f):
                os.remove(binding_f)
    return adata


def run(rna_h5ad, atac_h5ad, grn_dir, outdir, genome='hg38', species='Human', method='LINGER',
        activation='ReLU', networks=('cis',), n_cells=None, cells_file=None,
        keep_intermediate=False, skip_population=False):
    os.makedirs(outdir, exist_ok=True)
    if not outdir.endswith('/'):
        outdir = outdir + '/'
    if not grn_dir.endswith('/'):
        grn_dir = grn_dir + '/'

    adata_RNA = sc.read_h5ad(rna_h5ad)
    adata_ATAC = sc.read_h5ad(atac_h5ad)

    for col in ('sample', 'label'):
        if col not in adata_RNA.obs.columns:
            raise ValueError(f"adata_RNA.obs is missing the '{col}' column required by LINGER.")
    if adata_RNA.n_obs != adata_ATAC.n_obs or not (adata_RNA.obs_names == adata_ATAC.obs_names).all():
        raise ValueError("adata_RNA and adata_ATAC must have identical, barcode-aligned obs_names.")

    if 'trans' in networks and n_cells is None and cells_file is None:
        print('WARNING: --networks includes "trans" with no cell subset given. '
              'Per-cell trans-regulatory extraction writes a dense RE x TF matrix to disk for every '
              'cell and will not scale to a full dataset. Consider --n_cells/--cells_file.')

    start_all = time.monotonic()

    if not skip_population:
        _pseudobulk_and_preprocess(adata_RNA, adata_ATAC, grn_dir, genome, method, outdir)
        if 'trans' in networks:
            _train_population_model(adata_RNA, adata_ATAC, grn_dir, genome, method, activation, species, outdir)

    cell_RNA, cell_ATAC = _select_cells(adata_RNA, adata_ATAC, n_cells, cells_file)
    barcodes = _run_cell_level(cell_RNA, cell_ATAC, grn_dir, genome, method, outdir, networks)

    if 'cis' in networks:
        cis_adata = _assemble_cis_adata(outdir, barcodes, keep_intermediate)
        cis_adata.write(op.join(outdir, 'cell_level_cis_regulatory.h5ad'))
    if 'trans' in networks:
        trans_adata = _assemble_trans_adata(outdir, barcodes, keep_intermediate)
        trans_adata.write(op.join(outdir, 'cell_level_trans_regulatory.h5ad'))

    write_config(
        {'total_time': time.monotonic() - start_all, 'n_cells': len(barcodes), 'networks': list(networks)},
        file=op.join(outdir, 'runtime_all.yaml'),
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run LINGER and assemble cell-level GRN(s) into h5ad.")
    parser.add_argument("-r", "--rna_h5ad", type=str, required=True,
                         help="paired scRNA AnnData h5ad, obs must have 'sample' and 'label'")
    parser.add_argument("-a", "--atac_h5ad", type=str, required=True,
                         help="paired scATAC AnnData h5ad, barcode-aligned with rna_h5ad")
    parser.add_argument("-g", "--grn_dir", type=str, required=True,
                         help="path to LINGER's downloaded bulk reference data (GRNdir)")
    parser.add_argument("-o", "--outdir", type=str, required=True,
                         help="output directory for LINGER intermediates and final h5ad(s)")
    parser.add_argument("--genome", type=str, default="hg38", choices=["hg19", "hg38"])
    parser.add_argument("--species", type=str, default="Human")
    parser.add_argument("--method", type=str, default="LINGER", choices=["LINGER", "baseline"])
    parser.add_argument("--activation", type=str, default="ReLU", choices=["ReLU", "sigmoid", "tanh"])
    parser.add_argument("--networks", type=str, nargs="+", default=["cis"], choices=["cis", "trans"],
                         help="which cell-level network(s) to assemble into h5ad")
    parser.add_argument("--n_cells", type=int, default=None,
                         help="randomly subsample this many cells for cell-level extraction")
    parser.add_argument("--cells_file", type=str, default=None,
                         help="headerless text file of barcodes to restrict cell-level extraction to")
    parser.add_argument("--keep_intermediate", action="store_true",
                         help="keep LINGER's per-cell txt outputs instead of deleting them after assembly")
    parser.add_argument("--skip_population", action="store_true",
                         help="skip pseudobulk/preprocess/training steps (reuse a prior run's outdir)")
    args = parser.parse_args()

    run(args.rna_h5ad, args.atac_h5ad, args.grn_dir, args.outdir, genome=args.genome, species=args.species,
        method=args.method, activation=args.activation, networks=args.networks, n_cells=args.n_cells,
        cells_file=args.cells_file, keep_intermediate=args.keep_intermediate, skip_population=args.skip_population)
