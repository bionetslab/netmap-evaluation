reticulate::py_install("anndata")
library(optparse)
library(yaml)
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(Matrix)
library(ggplot2)
library(ggraph)
library(scatterpie)
library(ggrepel)
library(ggpubr)
library(igraph)
library(SingleCellExperiment)
options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages("BiocManager")
if (!requireNamespace("GenomeInfoDbData", quietly = TRUE)) {
  BiocManager::install("GenomeInfoDbData")
}
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")
remotes::install_github("cellgeni/sceasy")
source("~/netmap-evaluation/NeighbourNet/tests/script.new.R") #nolint
print("Source script loaded.")


run_neighbournet <- function(config, dataset_config) {
  # Load the dataset
  adatapath <- config$input_data
  outdir <- config$output_directory
  SeuratDisk::Convert(adatapath, dest = "h5seurat", overwrite = TRUE)
  h5seurat_path <- sub("\\.h5ad$", ".h5seurat", adatapath)
  obj <- SeuratDisk::LoadH5Seurat(h5seurat_path, verbose = FALSE)

  # Load priors?
  load("~/netmap-evaluation/NeighbourNet/data/gene.list.rda", envir = .GlobalEnv) #nolint
  load("~/netmap-evaluation/NeighbourNet/data/sig.graph.rda", envir = .GlobalEnv) #nolint
  load("~/netmap-evaluation/NeighbourNet/data/gr.graph.rda", envir = .GlobalEnv) #nolint
  load("~/netmap-evaluation/NeighbourNet/data/receptor.ppr.rda", envir = .GlobalEnv) #nolint
  # Preprocess
  rt.ppr <- get.ppr()                        # receptor‑target prior matrix
  genes  <- select.gene(obj, min.cells = 10) # QC → TF / target lists

  obj <- obj |>
    prepare.seurat(genes = genes$genes) |>   # scale + PCA
    prepare.graph() |>                       # 30‑NN graph
    select.cell(all=TRUE) |>                         # subsample
    prepare.reg(predictors = genes$tfs,      # local variance scaffolding
                responses  = genes$targets)

  # Run NeighbourNet
  responses <- genes$targets
  obj <- run.nn.reg(obj, responses = responses, return.p.val = TRUE) # nolint
  #results <- Misc(obj, "mod")

  #obj <- build.meta.network(obj)


  effect <- obj@misc$mod$effect
  pval <- obj@misc$mod$p.val

  source_genes <- dimnames(effect)[[1]]
  target_genes <- dimnames(effect)[[2]]
  cells <- dimnames(effect)[[3]]

  # Create all gene pair names
  gene_pairs <- expand.grid(source = source_genes, target = target_genes, stringsAsFactors = FALSE) # nolint
  colnames_pairs <- paste(gene_pairs$source, gene_pairs$target, sep = "_")

  # Reshape: for each cell, flatten the [target, source] matrix into a vector
  effect_mat <- sapply(cells, function(cell) {
    as.vector(effect[, , cell])
  })
  pval_mat <- sapply(cells, function(cell) {
    as.vector(pval[, , cell])
  })

  # Transpose so rows = cells, columns = gene pairs
  effect_mat <- t(effect_mat)
  colnames(effect_mat) <- colnames_pairs
  rownames(effect_mat) <- cells
  pval_mat <- t(pval_mat)
  colnames(pval_mat) <- colnames_pairs
  rownames(pval_mat) <- cells

  # Add to sce object
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(
      effect = Matrix::Matrix(effect_mat, sparse = TRUE),
      pval = Matrix::Matrix(pval_mat, sparse = TRUE)
    )
  )

  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  outpath <- paste0(outdir, "/nnet_results.h5ad")
  # Convert to AnnData format
  sceasy::convertFormat(sce, from = "sce", to = "anndata",
                        outFile = outpath)
  print(paste("Saved anndata file to ", outpath))
}

# Argument parsing
option_list <- list(
  make_option(c("-c", "--config"), type = "character",
              help = "Path to the YAML configuration file", metavar = "FILE"),
  make_option(c("-d", "--dataset_config"), type = "character",
              help = "Path to the dataset configuration file", metavar = "FILE")
)

'
Rscript ~/netmap-evaluation/src/methods/neighbourNet/neighbourNet.R -c ~/netmap-evaluation/configurations/neighbourNet/nnet_config.yaml -d ~/netmap-evaluation/configurations/data_simulation/config_easy.yaml
'
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$config)) {
  stop("Error: Please provide a YAML config file using the --config option.")
} else if (is.null(opt$dataset_config)) {
  stop("Error: Please provide a dataset config file using the --dataset_config option.") # nolint
}
# config <- yaml::read_yaml("~/netmap-evaluation/configurations/neighbourNet/nnet_config.yaml") #nolint
# dataset_config <- yaml::read_yaml("~/netmap-evaluation/configurations/data_simulation/config_easy.yaml") #nolint
config <- yaml::read_yaml(opt$config)
dataset_config <- yaml::read_yaml(opt$dataset_config)

run_neighbournet(config, dataset_config)