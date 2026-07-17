# %%
import scanpy as sc
import sys
import pandas as pd
import numpy as np
sys.path.append('/data/bionets/og86asub/netmap/scGeneRAI')
from scGeneRAI import scGeneRAI
import os
import networkx as nx
import os.path as op

if __name__ == "__main__":

    ex_data_us = pd.read_csv('/data/bionets/og86asub/netmap/netmap-evaluation/src/manuscript_synthetic/manuscript_results/figures/epi_top2000.csv')
    ex_data_descriptors = ex_data_us.iloc[:, 2002:2012].copy()
    #ex_data_us = pd.read_csv('/data_nfs/og86asub/netmap/netmap-evaluation/src/manuscript_synthetic/manuscript_results/figures/epi_top2000.csv')

    # ex_data_us = ex_data_us.iloc[:,2:802].copy()

    # means = ex_data_us.mean(axis=0)
    # sds = ex_data_us.std(axis=0)
    # ex_data_us = (ex_data_us-means)/sds


    #model =scGeneRAI()

    #model.fit(ex_data_us, nepochs = 100, model_depth =2, descriptors = ex_data_descriptors, early_stopping=True, device_name = 'cuda')
    #model.predict_networks(ex_data_us, descriptors = ex_data_descriptors, PATH = '.')



    res_files = os.listdir('./results')
    chunk_size = 300

    pivot_chunks = []
    edge_chunks = []

    for i in range(0, len(res_files), chunk_size):
        j = min(i + chunk_size, len(res_files))
        batch = res_files[i:j]
        batch_df = pd.concat([pd.read_csv(op.join('./results', res_file)).assign(cell_name=res_file.split('_')[-1].replace('.csv', '')) for res_file in batch], ignore_index=True)

        batch_df = batch_df.drop(batch_df.columns[0], axis=1)

        # collect unique edges for this chunk
        edge_chunks.append(batch_df[['source_gene', 'target_gene']].drop_duplicates())

        # pivot this chunk and store it
        chunk_pivot = batch_df.pivot_table(
            index='cell_name', columns=['source_gene', 'target_gene'], values='LRP', fill_value=0
        )
        pivot_chunks.append(chunk_pivot)

        print(chunk_pivot)

    # combine unique edges across all chunks
    unique_edges = pd.concat(edge_chunks, ignore_index=True).drop_duplicates()
    nl = [(row['source_gene'], row['target_gene']) for _, row in unique_edges.iterrows()]

    # combine pivoted chunks: rows (cell_name) are disjoint across chunks,
    # columns (gene pairs) may differ, so align via outer concat and fill 0
    pivot_table = pd.concat(pivot_chunks, axis=0).fillna(0)
    aa = pivot_table.values

    pivot_table.to_csv('./scgenerai_pandas.tsv', sep ='\t')


    ada2 = sc.AnnData(X = aa, obs = ex_data_descriptors.iloc[0:aa.shape[0]])
    ada2.write_h5ad('./scgenerai.h5ad')
    sc.pp.scale(ada2)
    sc.pp.normalize_total(ada2)
    sc.tl.pca(ada2)
    sc.pp.neighbors(ada2)
    sc.tl.umap(ada2)

    sc.pl.umap(ada2, save='_with_descriptors.pdf', color =['cell_type_epi', 'patient_id'])


    # ada = sc.AnnData(X = ex_data_us, obs = ex_data_descriptors)
    # sc.pp.scale(ada)
    # sc.pp.normalize_total(ada)
    # sc.tl.pca(ada)
    # sc.pp.neighbors(ada)
    # sc.tl.umap(ada)


    # sc.pl.umap(ada, save=True, color =['cell_type_epi', 'patient_id'])


