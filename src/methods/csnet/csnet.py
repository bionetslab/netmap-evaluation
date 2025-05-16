import numpy as np
from scipy.sparse import csr_matrix, lil_matrix
from scipy.stats import norm
import scanpy as sc
import os.path as op
import pandas as pd
import anndata
from csnet_config import CsNetConfig
import os

def run_csndm(config):

    ## Load config and setup outputs

    os.makedirs(config.output_directory, exist_ok = True)
    sc.settings.figdir = config.output_directory

    # copy configuration file to the folder
    config.write_yaml(yaml_file=op.join(config.output_directory, 'config.yaml'))
    
    # Load data from path in config
    adata = sc.read_h5ad(config.input_data)
    data = adata.X.T # because rows = genes and columns = cells (for this method)

    #### CSNDM method ###
  
    n1, n2 = data.shape

    # Defining upper and lower neighborhood bounds for gene expression values
    upper = np.zeros((n1, n2))
    lower = np.zeros((n1, n2))

    # Looping through each gene to define its local neighborhood
    for i in range(n1):
        s1 = np.sort(data[i, :])
        s2 = np.argsort(data[i, :])
        n0 = n2 - np.sum(np.sign(s1))
        h = round(config.boxsize / 2 * np.sum(np.sign(s1)))
        k = 0

        while k < n2:
            s = 0
            while k + s + 1 < n2 and s1[k + s + 1] == s1[k]:
                s += 1

            if s >= h:
                upper[i, s2[k:k+s+1]] = data[i, s2[k]]
                lower[i, s2[k:k+s+1]] = data[i, s2[k]]
            else:
                upper[i, s2[k:k+s+1]] = data[i, s2[int(min(n2 - 1, k + s + h))]]
                lower[i, s2[k:k+s+1]] = data[i, s2[int(max(n0 * (n0 > h), k - h))]]

            k += s + 1

    ndm = np.zeros((n1, n2))
    p = -norm.ppf(config.alpha)

    # Construct a cell-specific network for each cell
    for k in range(n2):
        # Check if gene expression values fall within the upper and lower neighborhood
        B = np.logical_and(data <= upper[:, k, np.newaxis],
                           data >= lower[:, k, np.newaxis])
        B = np.logical_and(B, data[:, k][:, np.newaxis] != 0)
        
        a = np.sum(B, axis=1, keepdims=True)
        aaT = np.dot(a, a.T)
        denom = np.sqrt((aaT * ((n2 - a) @ (n2 - a).T)) / (n2 - 1) + np.finfo(float).eps)
        csn = (np.dot(B, B.T) * n2 - aaT) / denom
        csn = csn > p

        # Degree = number of edges per gene (excluding self-connections)
        ndm[:, k] = np.sum(csn, axis=1) - np.diag(csn)
        #print(f"Cell {k+1} is completed")

    if config.normalize:
        ndm_sum = np.sum(ndm, axis=0, keepdims=True)
        ndm = ndm / ndm_sum
        mean_sum = np.mean(np.sum(np.sign(ndm), axis=0))
        ndm *= (mean_sum ** 2) / 2000 #  c=2000 (constant)

    if config.format ==".tsv":
        # Create column names
        row_names = [f"G{i}" for i in range(data.shape[0])]
        df = pd.DataFrame(ndm, columns=[f"Cell_{i}" for i in range(data.shape[1])], index=row_names)

        # Save as TSV file 
        df.to_csv(op.join(op.join(config.output_directory),config.filename+'.tsv'), sep='\t', index=True, header=True)
        
    else:
        # Create a sparse matrix 
        sparse_data = csr_matrix(ndm)
        
        # Create AnnData object
        adata = anndata.AnnData(X=sparse_data)

        # Add cell and gene names
        adata.obs_names =  [f"G{i}" for i in range(data.shape[0])] # Gene names
        adata.var_names = [f"Cell_{i}" for i in range(data.shape[1])]  # Cell names

        # Save as .h5ad file
        adata.write(op.join(op.join(config.output_directory),config.filename+".h5ad"))

    return ndm


def parse_args():
    parser = argparse.ArgumentParser(description="Argument parser with two configuration files.")
    
    parser.add_argument(
        '--config',
        type=str,
        default='configuration.json',
        help='Path to the main configuration file (default: configuration.json)'
    )

    return parser.parse_args()


if __name__ == "__main__":
    import argparse

    args = parse_args()
    print(f"Main Configuration File: {args.config}")
    #print(f"Dataset Configuration File: {args.dataset_config}")

    config = CsNetConfig.read_yaml(args.config)

    run_csndm(config)