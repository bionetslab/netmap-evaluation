suppressMessages(require(grn2gex))
suppressMessages(require(scMultiSim))
suppressMessages(require(data.table))
suppressMessages(require(optparse))
suppressMessages(library(yaml))
suppressMessages(require(SingleCellExperiment))
suppressMessages(library(sceasy))
suppressMessages(library(reticulate))

## Save as h5ad file
use_python('/opt/conda/bin/python')

use_condaenv('base')
print(py_config())
# 
# Define command line options
option_list <- list(
  make_option(
    c("-c", "--config"),
    type = "character",
    default = NULL,
    help = "Path to the YAML configuration file",
    metavar = "FILE"
  )
)
data_sim
# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if the config file is provided
if (is.null(opt$config)) {
  stop("Error: Please provide a YAML config file using the --config option.")
}

# Load the YAML configuration file
print(opt$config)
config <- yaml.load_file(opt$config)



###### DATA SIMULATION #######################################################################
# Define required variables within 'data_simulation'
required_vars_simulation <- c("n_cells", "n_celltypes", "edgelist", 'nodelist', "seed", "modification_type", "n_genes_per_module", 'n_celltypes', 'n_modules')

# Check for missing variables in the 'data_simulation' section
missing_vars <- required_vars_simulation[!required_vars_simulation %in% names(config)]
if (length(missing_vars) > 0) {
  stop(paste("Error: Missing required variables", paste(missing_vars, collapse = ", ")))
}

# # Check if 'input_grn' exists and is a valid file path
# if (!file.exists(config$edgelist)) {
#   stop(paste("Error: Specified file", config$edgelist, "does not exist."))
# }
# # Check if 'input_grn' exists and is a valid file path
# if (!file.exists(config$nodelist)) {
#   stop(paste("Error: Specified file", config$nodelist, "does not exist."))
# }


output_dir<-file.path('/usr/src/app/output')


edgelist<-fread(file.path('/usr/src/app/input', config$edgelist), sep='\t')
nodelist<-fread(file.path('/usr/src/app/input', config$nodelist), sep = '\t')


# config <- yaml.load_file("/home/bionets-og86asub/Documents/netmap/thenetmap/dockerize-grn2gex/config.yaml")
# edgelist<-fread(file.path("/home/bionets-og86asub/Documents/netmap/thenetmap/NetMap_LRP/data/simulation/collectri_subnetworks_2/net_135_44105/", config$edgelist), sep='\t')
# nodelist<-fread(file.path("/home/bionets-og86asub/Documents/netmap/thenetmap/NetMap_LRP/data/simulation/collectri_subnetworks_2/net_135_44105/", config$nodelist), sep = '\t')
# output_dir<-file.path("/home/bionets-og86asub/Documents/netmap/thenetmap/NetMap_LRP/data/simulation", config$dataset_id)

print('UPDATED')

dataset<-create_gex_data(net = edgelist, 
                         node_labels = nodelist, 
                         net_name = config$dataset_id,
                         nr_grns = config$n_celltypes,
                         nr_modules = config$n_modules,
                         nr_genes_per_module = config$n_genes_per_module,
                         seed = config$seed,
                         tf_effect = config$tf_effect,
                         weight_delta = config$weight_delta,
                         base_effect = config$base_effect,
                         mean = config$mean,
                         sd = config$sd,
                         disregulation_type = config$modification_type)


gex.dir<-file.path(output_dir,config$dataset_id )
gex.dir<-save_generated_data(dataset$net, dataset$counts, dataset$meta, gex_dir =gex.dir, disregulated_info = dataset$disregulated_regulators)

print(gex.dir)


sce <- SingleCellExperiment(list(counts=dataset$counts), colData=dataset$meta, rowData = DataFrame(genes = rownames(dataset$counts)))

out_path <- file.path(gex.dir, 'data.h5ad')

sce <- sceasy::convertFormat(sce, from="sce", to="anndata",
                      outFile=out_path)