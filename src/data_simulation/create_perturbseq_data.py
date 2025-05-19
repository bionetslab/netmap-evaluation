# %%
import scanpy as sc
import pandas as pd
import numpy as np
import random
import matplotlib.pyplot as plt
import warnings
import os.path as op
import os
import anndata
import rpy2.robjects as ro
import rpy2.robjects.pandas2ri as pandas2ri
pandas2ri.activate()
warnings.filterwarnings('ignore')
from sklearn.metrics import *
import itertools as itert
import sys

import yaml
def read_config(file):
    with open(file, "r") as f:
        config = yaml.safe_load(f)
    return config


def write_config(c, file):
    with open(file, "w") as handle:
        yaml.safe_dump(c, handle)


def check_missing(adata, gene_list):
    """Subset an AnnData object to retain only specified genes."""
    missing_genes = [gene for gene in gene_list if gene not in adata.var_names]
    return missing_genes

def preprocess_anndata(adata, ntop=2000):
    """Preprocess the AnnData object."""
    adata.var_names_make_unique()
    adata.raw = adata.copy()
    if adata.raw is None:  # Check if raw is not already set
        raise ValueError("adata.raw is None. Please assign raw counts to adata.raw before proceeding.")
    
    # Annotate mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    # Filter cells based on QC metrics
    adata = adata[adata.obs.n_genes_by_counts < 6000, :]
    adata = adata[adata.obs.pct_counts_mt < 20, :].copy()

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=10)

    # Normalize and log-transform the data
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Identify highly variable genes
    #print(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=ntop, flavor='seurat')
    adata = adata[:, adata.var.highly_variable]

    # Scale the data
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack")
    
    return adata

def preprocess_anndata_2(adata):
    """Preprocess the AnnData object."""
    adata.var_names_make_unique()
    adata.raw = adata.copy()

    adata.obs[['gene_1', 'gene_2']]  = adata.obs['guide_identity'].str.split('_',expand=True)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.scale(adata, max_value=10)

    # Scale the data
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack")

    # visualize
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    
    return adata
    
def visualize_anndata(adata, color_by):
    """Generate UMAP visualization of the AnnData object."""
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=color_by)

def get_subset(perturbed, random_subset):
    perturbed_3genes = perturbed[perturbed.obs['gene_1'].isin(random_subset)].copy()
    adata_preprocessed = preprocess_anndata(perturbed_3genes, ntop = 2000)
    missing = check_missing(adata_preprocessed, random_subset)
    counter = 2
    while len(missing)>0 and counter<=4:
        adata_preprocessed = preprocess_anndata(perturbed_3genes, ntop = 2000*counter)
        counter+=1
        missing = check_missing(adata_preprocessed, random_subset)
    return adata_preprocessed
    
def load_anndata(file_path):
    """Load an AnnData object from a file."""
    adata = sc.read_h5ad(file_path)
    adata.obs[['gene_1', 'gene_2']]  = adata.obs['guide_identity'].str.split('_',expand=True)
    adata = adata[adata.obs['gene_2'] == 'NegCtrl0']
    ctrl = adata[adata.obs['gene_1'].str.contains('NegCtrl')].copy()
    perturbed = adata[~adata.obs['gene_1'].str.contains('NegCtrl')].copy()
    return ctrl, perturbed

def compute_clustering_metric(perturbed, uniquely_perturbed_genes, k = 2):
    """
    Compute clusterin metrics for pairs of transcription factors.
    It subsets the adata object to the relevant condition, preprocesses
    the subset data and clusters. The clustering metrics are returned.
    Runtime dependend on list of genes submitted.
    """
    lili = []
    gg = ""
    for g in itert.combinations(uniquely_perturbed_genes, k):
        if gg != g[0]:
            #print(g[0])
            gg = g[0]
        subi  = get_subset(perturbed, g)
        try:
            silh = silhouette_score(subi.obsm['X_pca'], subi.obs.gene_1)
            cal = calinski_harabasz_score(subi.obsm['X_pca'], subi.obs.gene_1)
            dav_b = davies_bouldin_score(subi.obsm['X_pca'], subi.obs.gene_1)
            res = list(g) + [silh, cal, dav_b]
            lili.append(res)
        except:
            print(f'Failure on {g}')
    
    silhoutte_df = pd.DataFrame(lili)
    cols = [f"gene_{i}" for i in range(k)] + ['silhouette_score', 'calinski_harabasz_score', 'david_bouldin_index']
    silhoutte_df.columns = cols
    return silhoutte_df

def create_pseudobulks(adata, cell_type_key="cell_type", gemgroup_key="gemgroup"):
    if adata.raw is None:
        raise ValueError("adata.raw is None. Please assign raw counts to adata.raw before proceeding.")

    adata.obs["gemgroup"] = adata.obs["gemgroup"].astype("category")
    adata.obs["cell_type"] = adata.obs["gene_1"]
    adata.obs = adata.obs[['gemgroup', 'cell_type']]
    
    counts = adata.raw.X
    groups = adata.obs.groupby([cell_type_key, gemgroup_key])
    
    pseudobulks = [np.array(counts[[adata.obs_names.get_loc(i) for i in g.index], :].sum(axis=0)).flatten()
                   for _, g in groups]

    
    pb = pd.DataFrame(pseudobulks, columns=adata.raw.var_names, 
                        index=[f"{c}_{g}" for c, g in groups.groups.keys()])

    md = pd.DataFrame(pb.index, columns=['samples'])
    md.index = md["samples"]
    md["samples"] = md["samples"].str.split("_").str[0]
    md["samples"] = md["samples"].astype("category")
    
    return pb, md

def network_data(ctrl, data_config, perturbed_adata):
    
    # data for network
    ctrl.var_names_make_unique()
    
    ctrl_vec = [x for x in list(perturbed_adata.obs["gene_2"].unique())+list(perturbed_adata.obs["gene_1"].unique()) if "Ctrl" in x]
    if len(ctrl_vec) > 1:
        raise ValueError("Add the case of two control perturbation vectors")
    
    ctrl_adata = ctrl[
            (ctrl.obs["gene_1"] == ctrl_vec[0]) &  
            (ctrl.obs["gene_2"] == ctrl_vec[0])
            ].copy()

    # Check if raw data exists in both objects
    #print(perturbed_adata.raw)
    #print(ctrl_adata.raw)

    adata_networks = anndata.concat([perturbed_adata, ctrl_adata], join='inner')
    
    adata = sc.read_h5ad(data_config['data_path'])
    adata.var_names_make_unique()
    adata_subset = adata[:, adata_networks.var_names].copy()  # Keep only shared genes (var)
    adata_subset = adata_subset[adata_networks.obs_names, :].copy()  # Keep only shared cells (obs)

    return adata_subset

DE_edgeR_and_pb_Rcode = """
library(edgeR)

counts <- pb
metadata <- md
metadata$group <- as.factor(metadata$samples)

y <- DGEList(counts=counts, group=metadata$samples)

keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

y <- calcNormFactors(y)

design <- model.matrix(~ 0 + metadata$samples)
colnames(design) <- gsub("[^[:alnum:]_]", "", colnames(design))

control_group <- grep("Ctrl", colnames(design), value = TRUE)

experimental_groups <- setdiff(colnames(design), control_group)

if (length(experimental_groups) != 2) {
    stop("There should be exactly two experimental groups.")
}

y <- estimateDisp(y, design)

fit <- glmQLFit(y, design)

contrast1 <- makeContrasts(
  paste0(experimental_groups[1], " - ", control_group), 
  levels = design
)

contrast2 <- makeContrasts(
  paste0(experimental_groups[2], " - ", control_group), 
  levels = design
)

qlf1 <- glmQLFTest(fit, contrast = contrast1)
qlf2 <- glmQLFTest(fit, contrast = contrast2)

topGenes1 <- topTags(qlf1, n = net_size + 1)
topGenes2 <- topTags(qlf2, n = net_size + 1)

topGenes_df1 <- as.data.frame(topGenes1$table)
topGenes_df2 <- as.data.frame(topGenes2$table)
"""

DE_DESeq2_and_NO_pb_Rcode = """
suppressMessages(library(DESeq2))

# Prepare data
counts <- raw_counts
metadata <- meta_data
#print("in R")
#print(dim(counts))
#print(dim(metadata))
#print(head(head(metadata)))
metadata$group <- as.factor(metadata$cell_type)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ group)

# Run DESeq2
dds <- DESeq(dds)

# Extract control and experimental groups
control_group <- grep("Ctrl", levels(metadata$group), value = TRUE)
experimental_groups <- setdiff(levels(metadata$group), control_group)

if (length(experimental_groups) != 2) {
    stop("There should be exactly two experimental groups.")
}

# Perform DEG analysis for each experimental group vs control
res1 <- results(dds, contrast = c("group", experimental_groups[1], control_group))
res2 <- results(dds, contrast = c("group", experimental_groups[2], control_group))

#print(head(res1))

# Convert results to dataframes
topGenes_df1 <- as.data.frame(res1)
topGenes_df2 <- as.data.frame(res2)

# Filter by adjusted p-value (padj < 0.05) and order by absolute log2FoldChange
topGenes_df1 <- topGenes_df1[topGenes_df1$padj < 0.05, ]
topGenes_df2 <- topGenes_df2[topGenes_df2$padj < 0.05, ]

# Order by absolute log2FoldChange and take top N genes
topGenes_df1 <- topGenes_df1[order(abs(topGenes_df1$log2FoldChange), decreasing = TRUE), ]
topGenes_df2 <- topGenes_df2[order(abs(topGenes_df2$log2FoldChange), decreasing = TRUE), ]

# Select the top genes (net_size + 1) for each
topGenes_df1 <- head(topGenes_df1, net_size + 1)
topGenes_df2 <- head(topGenes_df2, net_size + 1)
"""

def DE_GRN(TF, topGenes_df, outdir, case):

    TF = TF
    tglist = list(topGenes_df.index)
    
    network_df = pd.DataFrame({
        "source": [TF] * len(tglist),  # Repeat TF1 for each gene in tglist
        "target": tglist  # List of target genes
    })
    
    filename = op.join(outdir, f"{TF}_network_{case}.csv")
    network_df.to_csv(filename, index=False) 
    return filename
    
def create_subsets(perturbed, genes_of_interest, data_config, netmap_config, ctrl):
    subset_tracker = []
    counter = 0
    while counter<10:
        random_subset = np.random.choice(genes_of_interest, size=data_config['number_tfs'], replace=False)
        dirname = f"{'_'.join(random_subset)}"
        if sorted(random_subset) in subset_tracker:
            continue
        print(dirname)
        counter+=1
        subset_tracker.append(sorted(random_subset))
        adata_preprocessed = get_subset(perturbed, random_subset)
        outdir = op.join(data_config['output_data_dir'], dirname)
        os.makedirs(outdir, exist_ok=True)
        
        if adata_preprocessed.raw is None:  # Check if raw is not already set
            raise ValueError("adata.raw is None. Please assign raw counts to adata.raw before proceeding.")
            
        adata_preprocessed.write_h5ad(op.join(outdir, 'data.h5ad'))

        # SAVE NETWORK DATA
        ctrl.var_names_make_unique()
        ctrl.raw = ctrl.copy()
        adata_networks = network_data(ctrl, data_config, perturbed_adata=adata_preprocessed)
        adata_networks = preprocess_anndata_2(adata_networks)
        adata_networks.write_h5ad(op.join(outdir, 'data_networks.h5ad'))

        # GET DE GRN wilcoxon
        adata = adata_networks.copy()
        adata = preprocess_anndata_2(adata)
        #print(adata)
        #print(adata.obs.columns)
        adata.obs["cell_type"] = adata.obs["gene_1"]
        #print(adata.obs["cell_type"].unique())
        groups__ = adata.obs["cell_type"].unique()
        groups_ = [g for g in groups__ if "Ctrl" not in g]
        control_ = [g for g in groups__ if "Ctrl" in g]
        #print(f"groups_: {groups_}")
        #print(f"control_: {control_}")

        sc.tl.rank_genes_groups(adata, groupby='cell_type', groups=groups_[:1], reference=control_[0], method='wilcoxon')
        l1 = adata.uns['rank_genes_groups']['names'][groups_[:1][0]][:data_config['net_size']]
        topGenes_df1 = pd.DataFrame(l1, columns=['Gene'], index=l1)
        #print(topGenes_df1)

        sc.tl.rank_genes_groups(adata, groupby='cell_type', groups=groups_[1:], reference=control_[0], method='wilcoxon')
        l2 = adata.uns['rank_genes_groups']['names'][groups_[1:][0]][:data_config['net_size']]
        topGenes_df2 = pd.DataFrame(l2, columns=['Gene'], index=l2)
        #print(topGenes_df2)

        tf1 = groups_[:1][0]
        network_1_file = DE_GRN(tf1, topGenes_df1, outdir, case = "wilcoxon_no_pb" )
        netmap_config['evaluation'] = {}
        netmap_config['evaluation']['gene_1'] = tf1
        netmap_config['evaluation']['network_1_0'] = network_1_file
        tf2 = groups_[1:][0]
        network_2_file = DE_GRN(tf2, topGenes_df2, outdir, case = "wilcoxon_no_pb")
        netmap_config['evaluation']['gene_2'] = tf2
        netmap_config['evaluation']['network_2_0'] = network_2_file

        
        # GET DE GRN DE_DESeq2_and_NO_pb
        pb, md = create_pseudobulks(adata_networks)
        #print(pb.shape)
        #print(md.shape)
        pb = pb.T
        #print(f"class pb:{type(pb)}")
        #print(pb.dtypes)
        ro.r.assign('pb', pandas2ri.py2rpy(pb))
        ro.r.assign('md', pandas2ri.py2rpy(md))
        ro.r.assign('net_size', data_config['net_size'])

        print("runing DE_edgeR_and_pb_Rcode")
        ro.r(DE_edgeR_and_pb_Rcode)

        topGenes_df1 = ro.r('topGenes_df1')
        topGenes_df2 = ro.r('topGenes_df2')
        experimental_groups = ro.r('experimental_groups')
        #print(topGenes_df1)
        #print(topGenes_df2)
        print(experimental_groups)
        tf1 = list(experimental_groups)[0].split("metadatasamples")[-1]
        tf2= list(experimental_groups)[1].split("metadatasamples")[-1]
        #print(f"tf1: {tf1} and tf1: {tf2}")
        # END GET DE GRN
        topGenes_df1 = pd.DataFrame(topGenes_df1)
        topGenes_df2 = pd.DataFrame(topGenes_df2)
        print(topGenes_df1.index)

        network_1_file = DE_GRN(tf1, topGenes_df1, outdir, case = "DE_edgeR_and_pb" )
        #netmap_config['evaluation'] = {}
        netmap_config['evaluation']['gene_1'] = tf1
        netmap_config['evaluation']['network_1_1'] = network_1_file
        
        network_2_file = DE_GRN(tf2, topGenes_df2, outdir, case = "DE_edgeR_and_pb")
        netmap_config['evaluation']['gene_2'] = tf2
        netmap_config['evaluation']['network_2_1'] = network_2_file

        # GET DE GRN DE_DESeq2_and_NO_pb
        #print(adata_networks)
        raw_counts = adata_networks.raw.X
        # Convert sparse to dense
        raw_counts = pd.DataFrame(raw_counts.todense(), index=adata_networks.obs_names, columns=adata_networks.var_names)  
        # Extract metadata
        meta_data = adata_networks.obs  
        meta_data = meta_data[['cell_type']]  

        ro.r.assign('meta_data', pandas2ri.py2rpy(meta_data))
        raw_counts = raw_counts.T
        ro.r.assign('raw_counts', pandas2ri.py2rpy(raw_counts))
        ro.r(DE_DESeq2_and_NO_pb_Rcode)
        topGenes_df1 = ro.r('topGenes_df1')
        topGenes_df2 = ro.r('topGenes_df2')
        experimental_groups = ro.r('experimental_groups')

        topGenes_df1 = pd.DataFrame(topGenes_df1)
        topGenes_df2 = pd.DataFrame(topGenes_df2)

        network_1_file = DE_GRN(tf1, topGenes_df1, outdir, case = "DE_DESeq2_and_NO_pb" )
        netmap_config['evaluation']['gene_1'] = tf1
        netmap_config['evaluation']['network_1_2'] = network_1_file
        network_2_file = DE_GRN(tf2, topGenes_df2, outdir, case = "DE_DESeq2_and_NO_pb")
        netmap_config['evaluation']['gene_2'] = tf2
        netmap_config['evaluation']['network_2_2'] = network_2_file


        # Write the updated netmap config to file
        os.makedirs(data_config['netmap_configuration_dir'], exist_ok=True)
        config_file =  op.join(data_config['netmap_configuration_dir'], f"{dirname}.yaml")
        netmap_config['data']['input_data'] = op.join(outdir, 'data.h5ad')
        netmap_config['results']['output_directory'] = op.join(data_config['results_dir'], dirname )
        write_config(netmap_config, config_file)
    
    return adata_preprocessed, outdir

# %%
def parse_args():
    parser = argparse.ArgumentParser(description="Argument parser with two configuration files.")
    
    parser.add_argument(
        '--config',
        type=str,
        default='configuration.json',
        help='Path to the main configuration file (default: configuration.json)'
    )
    parser.add_argument(
        '--dataset_config',
        type=str,
        default='dataset_configuration.json',
        help='Path to the dataset configuration file (default: dataset_configuration.json)'
    )
    
    return parser.parse_args()


if __name__ == "__main__":
    import argparse

    args = parse_args()
    print(f"Main Configuration File: {args.config}")

    
    data_config = read_config(
        file=args.config
    )
    
    ctrl, perturbed = load_anndata(data_config['data_path'])
    tfs = pd.read_csv(data_config['tf_file'], sep = '\t')

    netmap_config = read_config(data_config['netmap_base_config'])


    uniquely_perturbed_genes = perturbed.obs['gene_1'].unique()
    # Stable perturbations according to UMAP representation in paper
    genes_of_interest = [
            'CEBPE', 'KLF1', 'FOXA3', 'TBX2',
            'MAPK1', 'CEBPA']
    # These do not fall within stable GIs according to authors.
    unstable_perturbations = ['SNAI1', 'BCL2L11', 'PRDM1', 'ETS2', 'FOXA1' , 'FOSB' , 'JUN', 'IRF1','TBX3','AHR','FOXL2', 'ZBTB25', 'POU3F2','PLK4',  'LHX1']


    # Subset list of TFs from privided gene annotation
    tfs[tfs['gene type'] == 'transcription factor'].iloc[:, 0]
    # Intersect with the list of perturbed genes to get the list of relevant TFs
    uniquely_perturbed_tfs  = list(set(tfs[tfs['gene type'] == 'transcription factor'].iloc[:, 0]).intersection(set(uniquely_perturbed_genes)))
    
    create_subsets(perturbed, genes_of_interest, data_config, netmap_config, ctrl )
    
    
    
    silhoutte_df = compute_clustering_metric(perturbed, genes_of_interest, k = 2)
    silhoutte_df = silhoutte_df.sort_values('david_bouldin_index', ascending=False)
    os.makedirs(op.join(data_config['results_dir']), exist_ok = True)
    silhoutte_df.to_csv(op.join(data_config['results_dir'], 'clustering_metrics_2tfs.tsv'), sep= '\t')


