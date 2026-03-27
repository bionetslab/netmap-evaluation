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

    directory_name = '/data_nfs/og86asub/netmap/netmap-evaluation/netmap/data/cd8_mouse/'
    sample_id = '72h/GSM8286681_10XSC009-02-RNA'
    condition_name = 'exhausted_d3_'

    adata_72 = create_anndata_object(directory_name, sample_id, condition_name)

    directory_name = '/data_nfs/og86asub/netmap/netmap-evaluation/netmap/data/cd8_mouse/96h'
    sample_id = 'GSM8286680_10XSC009-01-RNA'
    condition_name = 'exhausted_d4_'

    adata_96 = create_anndata_object(directory_name, sample_id, condition_name)

    directory_name = '/data_nfs/og86asub/netmap/netmap-evaluation/netmap/data/cd8_mouse/active'
    sample_id = 'GSM8286683_10XSC011-02-RNA'
    condition_name = 'activated_d1_'

    adata_active = create_anndata_object(directory_name, sample_id, condition_name)

    directory_name = '/data_nfs/og86asub/netmap/netmap-evaluation/netmap/data/cd8_mouse/naive'
    sample_id = 'GSM8286682_10XSC011-01-RNA'
    condition_name = 'naive_d4_'

    adata_naive = create_anndata_object(directory_name, sample_id, condition_name)

    barcode = pd.read_csv('/data_nfs/og86asub/netmap/netmap-evaluation/netmap/data/cd8_mouse/ga_an0682_10x_gex_exvivo_p14_gp33_stim_int_filtered_cell_barcodes.txt')
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

    
    pangalao = pd.read_csv('/data_nfs/og86asub/netmap/netmap-evaluation/netmap/resources/panglaodb_mouse.csv')
    high_ui = pangalao[pangalao['UI']>0.3]['Official gene symbol'].values
    print(high_ui)
    adata.var['gene_name_upper'] = adata.var.index.str.upper()
    adata = adata[:, ~adata.var['gene_name_upper'].isin(high_ui)].copy()
    print(adata.var.gene_name_upper)
    print(f'Anndata after filtering: {adata.shape}')

    adata = adata[:, adata.var.highly_variable == True].copy()
    adata = adata[:, adata.var.pct_dropout_by_counts<97]
    
    print(adata.var.index)
    
    adata.write_h5ad('/data_nfs/og86asub/netmap/netmap-evaluation/netmap/data/blood/reprocessed/cd8tcells.h5ad')
    
    print(f'Anndata after filtering: {adata.shape}')
    gene_names = np.array(adata.var.index)
    data_tensor = adata.X


    if scs.issparse(data_tensor):
        data_tensor = torch.tensor(data_tensor.todense(), dtype=torch.float32)
    else:
        data_tensor = torch.tensor(data_tensor, dtype=torch.float32)

    model_zoo = create_model_zoo(data_tensor,  n_models=10, n_epochs=10000, model_type='NBAutoencoder', latent_dim= 8, dropout_rate=0.1, hidden_dim = [64] )
    grn_adata = inferrence(model_zoo, data_tensor.cuda(), gene_names, xai_method='GradientShap', background_type = 'zeros', backing_file=op.join(output_dir, 'grn_adata.h5'))
    
    model_name = 'cd8_tcells'
    #Save anndata obs to grn obs for reference
    grn_adata.obs = adata.obs
    grn_adata.write_h5ad( op.join(output_dir, f'{model_name}_grn.h5ad'))
    
    grn_adata.var.to_csv(op.join(output_dir, f'{model_name}_var.tsv'), header = '\t')
    # save the original obs
    adata.obs.to_csv(op.join(output_dir, f'{model_name}_obs.tsv'), header = '\t')




