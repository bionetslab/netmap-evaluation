import anndata as ad
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import combinations
from matplotlib.pyplot import rc_context
import numpy as np

# Number of cores to use
ncores = 64

import sys

import anndata
from netmap.downstream import final_downstream


import warnings

from netmap.utils.data_utils import *
from netmap.utils.tf_utils import *
from netmap.utils.netmap_config import NetmapConfig

from netmap.model.train_model import create_model_zoo
from netmap.grn.inferrence import inferrence, inferrence_model_wise
from netmap.masking.internal import *
from netmap.masking.external import *

from netmap.downstream.edge_selection import *
from netmap.downstream.clustering import *
from netmap.downstream.final_downstream import *

import scipy.sparse as scs
import torch
import pandas as pd

import os.path as op
import time

def create_anndata_object(directory, sample_id, condition_name):

    adata1 = sc.read_mtx(op.join(directory, f'{sample_id}_matrix.mtx.gz'))
    adata1_bc  = pd.read_csv(op.join(directory, f'{sample_id}_barcodes.tsv.gz'), header = None)
    adata1_f = pd.read_csv(op.join(directory, f'{sample_id}_features.tsv.gz'), header = None)

    adata1  = adata1.T

    adata1_f.columns = ['gene_name']
    adata1_bc.columns = ['cell_id']

    adata1.var['gene_name'] = adata1_f['gene_name'].values
    adata1.obs['cell_id'] = adata1_bc['cell_id'].values

    adata1.obs['condition'] = condition_name
    adata1.obs['barcode'] = adata1.obs['condition']+adata1.obs['cell_id']

    adata1.var = adata1.var.set_index('gene_name')

    return adata1



if __name__ == '__main__':

    # set output directory
    output_dir = "/data_nfs/og86asub/netmap/netmap-evaluation/results/cd8_usecase_low_ui"

    # Create both directories regardless of whether they already exist
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.join(output_dir, "networks"), exist_ok=True)

    directory_name = '/data_nfs/og86asub/netmap/netmap-evaluation/data/cd8_mouse/'
    sample_id = '72h/GSM8286681_10XSC009-02-RNA'
    condition_name = 'exhausted_d3_'

    adata_72 = create_anndata_object(directory_name, sample_id, condition_name)

    directory_name = '/data_nfs/og86asub/netmap/netmap-evaluation/data/cd8_mouse/96h'
    sample_id = 'GSM8286680_10XSC009-01-RNA'
    condition_name = 'exhausted_d4_'

    adata_96 = create_anndata_object(directory_name, sample_id, condition_name)

    directory_name = '/data_nfs/og86asub/netmap/netmap-evaluation/data/cd8_mouse/active'
    sample_id = 'GSM8286683_10XSC011-02-RNA'
    condition_name = 'activated_d1_'

    adata_active = create_anndata_object(directory_name, sample_id, condition_name)

    directory_name = '/data_nfs/og86asub/netmap/netmap-evaluation/data/cd8_mouse/naive'
    sample_id = 'GSM8286682_10XSC011-01-RNA'
    condition_name = 'naive_d4_'

    adata_naive = create_anndata_object(directory_name, sample_id, condition_name)

    barcode = pd.read_csv('/data_nfs/og86asub/netmap/netmap-evaluation/data/cd8_mouse/ga_an0682_10x_gex_exvivo_p14_gp33_stim_int_filtered_cell_barcodes.txt')
    barcode.columns = ['barcode']
    

    adata = ad.concat([adata_72, adata_96, adata_active, adata_naive])
    adata = adata[adata.obs['barcode'].isin(barcode.barcode)]

    #print(adata.var_names)
    #adata.var = adata.var.set_index('gene_name')


    adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    # mitochondrial genes, "MT-" for human, "Mt-" for mouse
    adata.var["mt"] = adata.var_names.str.startswith("Mt-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
    # hemoglobin genes
    adata.var["hb"] = adata.var_names.str.contains("^Hb[^(P)]")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

    sc.pp.highly_variable_genes(adata, n_top_genes=4000, flavor='seurat')

    
    # pangalao = pd.read_csv('/data_nfs/og86asub/netmap/netmap-evaluation/netmap/resources/panglaodb_mouse.csv')
    # high_ui = pangalao[pangalao['UI']>0.3]['Official gene symbol'].values
    # print(high_ui)
    # adata.var['gene_name_upper'] = adata.var.index.str.upper()
    # adata = adata[:, ~adata.var['gene_name_upper'].isin(high_ui)].copy()
    # print(adata.var.gene_name_upper)
    # print(f'Anndata after filtering: {adata.shape}')

    # adata = adata[:, adata.var.highly_variable == True].copy()
    # adata = adata[:, adata.var.pct_dropout_by_counts<97]
    
    # print(adata.var.index)

    # Genes/markers that should always survive filtering, no matter what
    predefined_markers = ['Cd3e', 'Cd19', 'Ptprc', 'Epcam', 'Pdcd1', 'Lag3', 'Havcr2', 'Tigit', 
                            'Tox', 'Eomes', 'Maf', 'Bach1', 'Cd200r1', 'Ikzf2', 'Tigit', 'Fasl', 'Ccl3',
                            'Ccl4', 'Ccl5','Gzma','Gzmb', 'Gzmc', 'Gzmk', 'Gzmm', 'Ifng', 'Prf1', 'Ctsd', 
                            'Ctsw', 'Tnf', 'Nkg7', 'Klrg1', 'Cd40lg', 'Tgfb1', 'Tgfb2',
                            'Tgfb3', 'Csf1', 'Csf2', 'Lif', 'Osm', 'Lta', 'Cd69', 'Il2ra',
                            'Hla-Dr1', 'Cd44', 'Prdm1', 'Tbx21', 'Tnfrsf9', 'Nfatc1', 'Nfatc2',
                            'Sell', 'Ccr7', 'Tcf7', 'Lef1', 'Il7r', 'Slamf6','Klf2','Klf3']
    predefined_markers_upper = pd.Index(predefined_markers).str.upper()

    pangalao = pd.read_csv('/data_nfs/og86asub/netmap/netmap-evaluation/netmap/resources/panglaodb_mouse.csv')
    high_ui = pangalao[pangalao['UI'] > 0.3]['Official gene symbol'].values
    print(high_ui)

    adata.var['gene_name_upper'] = adata.var.index.str.upper()
    is_marker = adata.var['gene_name_upper'].isin(predefined_markers_upper)

    # Step 1: high-UI filter — skip for markers
    keep = is_marker | ~adata.var['gene_name_upper'].isin(high_ui)
    adata = adata[:, keep].copy()
    print(adata.var.gene_name_upper)
    print(f'Anndata after filtering: {adata.shape}')

    # Recompute marker mask after subsetting
    is_marker = adata.var['gene_name_upper'].isin(predefined_markers_upper)

    # Step 2: highly_variable filter — skip for markers
    keep = is_marker | (adata.var.highly_variable == True)
    adata = adata[:, keep].copy()

    # Recompute again
    is_marker = adata.var['gene_name_upper'].isin(predefined_markers_upper)

    # Step 3: pct_dropout filter — skip for markers
    keep = is_marker | (adata.var.pct_dropout_by_counts < 97)
    adata = adata[:, keep].copy()

    print(adata.var.index)

    
    adata.write_h5ad('/data_nfs/og86asub/netmap/netmap-evaluation/data/blood/reprocessed/cd8tcells_extra_genes.h5ad')
    
    time_tracker  = []
    start = time.monotonic()
    
    print(f'Anndata after filtering: {adata.shape}')
    gene_names = np.array(adata.var.index)
    data_tensor = adata.X


    if scs.issparse(data_tensor):
        data_tensor = torch.tensor(data_tensor.todense(), dtype=torch.float32)
    else:
        data_tensor = torch.tensor(data_tensor, dtype=torch.float32)

    model_zoo = create_model_zoo(data_tensor,  n_models=10, n_epochs=10000, model_type='NBAutoencoder', latent_dim= 8, dropout_rate=0.1, hidden_dim = [64] )
    grn_adata = inferrence(model_zoo, data_tensor.cuda(), gene_names, xai_method='GradientShap', background_type = 'zeros', backing_file=op.join(output_dir, 'grn_adata.h5'))
    
    model_name = 'cd8_tcells_extra_genes'
    #Save anndata obs to grn obs for reference
    grn_adata.obs = adata.obs
    grn_adata.write_h5ad( op.join(output_dir, f'{model_name}_grn.h5ad'))
    
    grn_adata.var.to_csv(op.join(output_dir, f'{model_name}_var.tsv'), header = '\t')
    # save the original obs
    adata.obs.to_csv(op.join(output_dir, f'{model_name}_obs.tsv'), header = '\t')

    stop = time.monotonic()
    time_tracker.append([model_name, stop-start])

    time_tracker_df = pd.DataFrame(time_tracker)
    time_tracker_df.to_csv(op.join(output_dir, 'time_tracker.tsv'))


