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


install_seuratdisk <- function() {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  remotes::install_github("mojaveazure/seurat-disk")
}
source("~/netmap-evaluation/NeighbourNet/tests/script.new.R") #nolint

run_neighbournet <- function(config, dataset_config) {
  # Load the dataset 
  adatapath <- config$input_data
  SeuratDisk::Convert(adatapath, dest = "h5seurat", overwrite = TRUE)
  h5seurat_path <- sub("\\.h5ad$", ".h5seurat", adatapath)
  obj <- SeuratDisk::LoadH5Seurat(h5seurat_path, verbose = FALSE)

  # Load priors?
  load("~/netmap-evaluation/NeighbourNet/data/gene.list.rda") #nolint
  load("~/netmap-evaluation/NeighbourNet/data/sig.graph.rda") #nolint
  load("~/netmap-evaluation/NeighbourNet/data/gr.graph.rda") #nolint
  load("~/netmap-evaluation/NeighbourNet/data/receptor.ppr.rda") #nolint

  # Preprocess
  rt.ppr <- get.ppr()                        # receptor‑target prior matrix
  genes  <- select.gene(obj, min.cells = 10) # QC → TF / target lists

  obj <- obj |>
    prepare.seurat(genes = genes$genes) |>   # scale + PCA
    prepare.graph() |>                       # 30‑NN graph
    select.cell() |>                         # subsample
    prepare.reg(predictors = genes$tfs,      # local variance scaffolding
                responses  = genes$targets)

  # Run NeighbourNet
  responses <- genes$targets[1:10]       # for debugging/testing, delete later
  obj <- run.nn.reg(obj, responses = responses, return.p.val = TRUE)
  results <- Misc(obj, "mod")
  
  obj2 <- build.meta.network(obj)
}

# Argument parsing
option_list <- list(
  make_option(c("-c", "--config"), type = "character",
              help = "Path to the YAML configuration file", metavar = "FILE"),
  make_option(c("-d", "--dataset_config"), type = "character",
              help = "Path to the dataset configuration file", metavar = "FILE")
)

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

run_neighbourNet(config, dataset_config)