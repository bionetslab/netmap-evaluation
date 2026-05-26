
import anndata as ad
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import combinations
from matplotlib.pyplot import rc_context
import numpy as np
import sys
import os.path as op

import anndata



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

import scipy.sparse as scs
import torch

import time


def train_model(adata, output_dir, model_name):

   
    gene_names = np.array(adata.var.index)
    data_tensor = adata.X # Log normalized, but not standardized data.


    if scs.issparse(data_tensor):
        data_tensor = torch.tensor(data_tensor.todense(), dtype=torch.float32)
    else:
        data_tensor = torch.tensor(data_tensor, dtype=torch.float32)

    model_zoo = create_model_zoo(data_tensor,  n_models=10, n_epochs=10000, model_type='NBAutoencoder', latent_dim= 8, dropout_rate=0.1, hidden_dim = [64] )

    grn_adata = inferrence(model_zoo, data_tensor.cuda(), gene_names, xai_method='GradientShap', background_type = 'zeros', backing_file=op.join(output_dir, f'{model_name}.parquet'), return_in_memory=False)
    
    #Save anndata obs to grn obs for reference
    grn_adata.obs = adata.obs
    grn_adata.write_h5ad( op.join(output_dir, f'{model_name}_grn.h5ad'))
    
    grn_adata.var.to_csv(op.join(output_dir, f'{model_name}_var.tsv'), header = '\t')
    # save the original obs
    adata.obs.to_csv(op.join(output_dir, f'{model_name}_obs.tsv'), header = '\t')



if __name__=='__main__':

    time_tracker = []
    # define your output dir.
    output_dir = "/data_nfs/og86asub/netmap/netmap-evaluation/netmap/case_studies/blood"
    os.makedirs(output_dir, exist_ok=True)
    basedir = '/data_nfs/og86asub/netmap/netmap-evaluation/'

    
    # start = time.monotonic()
    # model_name = 'bd-rhap-rep1-X'
    # model_output_dir = op.join(output_dir, model_name)
    # os.makedirs(model_output_dir, exist_ok=True)
    # adata = sc.read_h5ad(op.join(basedir, 'data/blood/reprocessed/bd-rhap-rep1/bd-rhap-rep1.h5ad'))
    # print(model_name)
    # print(adata)
    # train_model(adata, output_dir=model_output_dir, model_name = model_name)
    # stop = time.monotonic()
    # time_tracker.append([model_name, stop-start])

    # time_tracker_df = pd.DataFrame(time_tracker)
    # time_tracker_df.to_csv(op.join(output_dir, 'time_tracker.tsv'))




    # try:
    #     start = time.monotonic()
    #     model_name = 'bd-rhap-rep2-X'
    #     model_output_dir = op.join(output_dir, model_name)
    #     os.makedirs(model_output_dir, exist_ok=True)
    #     print(model_name)
    #     adata = sc.read_h5ad(op.join(basedir, 'data/blood/reprocessed/bd-rhap-rep2/bd-rhap-rep2.h5ad'))
    #     print(adata)
    #     train_model(adata, output_dir=model_output_dir, model_name = model_name)
    #     stop = time.monotonic()
    #     time_tracker.append([model_name, stop-start])

    #     time_tracker_df = pd.DataFrame(time_tracker)
    #     time_tracker_df.to_csv(op.join(output_dir, 'time_tracker.tsv'))
    # except Exception:
    #     pass
   


    # try:
    #     start = time.monotonic()
    #     model_name = '10x-rep1-kallisto-cellbender-X'
    #     print(model_name)
    #     model_output_dir = op.join(output_dir, model_name)
    #     os.makedirs(model_output_dir, exist_ok=True)
    #     adata = sc.read_h5ad(op.join(basedir, 'data/blood/reprocessed/10x-rep1-kallisto-cellbender/10x-rep1-kallisto-cellbender.h5ad'))
    #     print(adata)
    #     train_model(adata, output_dir=model_output_dir, model_name = model_name)
    #     stop = time.monotonic()
    #     time_tracker.append([model_name, stop-start])

    #     time_tracker_df = pd.DataFrame(time_tracker)
    #     time_tracker_df.to_csv(op.join(output_dir, 'time_tracker.tsv'))
    # except Exception:
    #     pass
    


    try:
        start = time.monotonic()
        model_name = 'blood-10x-rep2-kallisto-cellbender-X'
        print(model_name)
        model_output_dir = op.join(output_dir, model_name)
        os.makedirs(model_output_dir, exist_ok=True)
        adata = sc.read_h5ad(op.join(basedir, 'data/blood/reprocessed/10x-rep2-kallisto-cellbender/10x-rep2-kallisto-cellbender.h5ad'))
        print(adata)
        train_model(adata, output_dir=model_output_dir, model_name = model_name)
        stop = time.monotonic()
        time_tracker.append([model_name, stop-start])

        time_tracker_df = pd.DataFrame(time_tracker)
        time_tracker_df.to_csv(op.join(output_dir, 'time_tracker.tsv'))
    except Exception:
        pass


    time_tracker_df = pd.DataFrame(time_tracker)
    time_tracker_df.to_csv(op.join(output_dir, 'time_tracker.tsv'))


    # try:
    #     for i in range(10):
    #         start = time.monotonic()
    #         model_name = f'10x-rep1-kallisto-cellbender-X-{i}'
    #         print(model_name)
    #         model_output_dir = op.join(output_dir, model_name)
    #         os.makedirs(model_output_dir, exist_ok=True)
    #         adata = sc.read_h5ad(op.join(basedir, 'data/blood/reprocessed/10x-rep1-kallisto-cellbender/10x-rep1-kallisto-cellbender.h5ad'))
    #         print(adata)
    #         train_model(adata, output_dir=model_output_dir, model_name = model_name)
    #         stop = time.monotonic()
    #         time_tracker.append([model_name, stop-start])

    #         time_tracker_df = pd.DataFrame(time_tracker)
    #         time_tracker_df.to_csv(op.join(output_dir, 'time_tracker_repeat.tsv'))
    # except Exception:
    #     pass
    


    # try:
    #     for i in range(10):
    #         start = time.monotonic()
    #         model_name = f'bd-rhap-rep2-X-{i}'
    #         model_output_dir = op.join(output_dir, model_name)
    #         os.makedirs(model_output_dir, exist_ok=True)
    #         print(model_name)
    #         adata = sc.read_h5ad(op.join(basedir, 'data/blood/reprocessed/bd-rhap-rep2/bd-rhap-rep2.h5ad'))
    #         print(adata)
    #         train_model(adata, output_dir=model_output_dir, model_name = model_name)
    #         stop = time.monotonic()
    #         time_tracker.append([model_name, stop-start])

    #         time_tracker_df = pd.DataFrame(time_tracker)
    #         time_tracker_df.to_csv(op.join(output_dir, 'time_tracker_repeat.tsv'))
    # except Exception:
    #     pass