load_packages <- function() {
    pkgs <- c("Seurat","SeuratData", "SeuratDisk","dplyr","Matrix","ggplot2","ggraph",
            "scatterpie","ggrepel","ggpubr","igraph","optparse","yaml")
    if(any(miss <- !pkgs %in% installed.packages()[,1]))
        install.packages(pkgs[miss], repos = "https://cloud.r-project.org")
    invisible(lapply(pkgs, library, character.only = TRUE))

    # Source NeighbourNet (in develop)
    source("tests/script.new.R")
}

run_neighbourNet <- function(config, dataset_config) {

    # Load the dataset
    adatapath <- config$input_data
    Convert(adatapath, dest = "h5seurat", overwrite = TRUE)
    obj <- LoadH5Seurat(adatapath, verbose = FALSE)


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
    obj <- run.nn.reg(obj, responses = responses, return.p.val = T)

}


#### Main function to run the neighbourNet method
#' @param config Path to the YAML configuration file
#' @param dataset_config Path to the dataset configuration file

load_packages()

option_list <- list(
  make_option(
    c("-c", "--config"),
    type = "character",
    default = NULL,
    help = "Path to the YAML configuration file",
    metavar = "FILE"
  ),
  make_option(
    c("-d", "--dataset_config"),
    type = "character",
    default = NULL,
    help = "Path to the dataset configuration file",
    metavar = "FILE"
  ),
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$config)) {
  stop("Error: Please provide a YAML config file using the --config option.")
}else if (is.null(opt$dataset_config)) {
  stop("Error: Please provide a dataset config file using the --dataset_config option.")
}

config <- yaml.load_file(opt$config)
dataset_config <- yaml.load_file(opt$dataset_config)

run_neighbourNet(config, dataset_config)
