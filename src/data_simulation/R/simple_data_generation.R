suppressMessages(require(grn2gex))
suppressMessages(require(scMultiSim))
suppressMessages(require(data.table))
suppressMessages(require(optparse))
suppressMessages(library(yaml))
suppressMessages(require(SingleCellExperiment))
suppressMessages(library(sceasy))
suppressMessages(library(reticulate))


py_require()

# Define command line options
option_list <- list(
  make_option(
    c("-c", "--config"),
    type = "character",
    default = NULL,
    help = "Path to the YAML configuration file",
    metavar = "FILE"
  ),
    make_option(
    c("-d", "--data_output_dir"),
    type = "character",
    default = NULL,
    help = "Path to the YAML configuration file",
    metavar = "FILE"
  )
)


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


print(opt$data_output_dir)
###### DATA SIMULATION #######################################################################
# Define required variables within 'data_simulation'
required_vars_simulation <- c("n_cells", "n_celltypes", "edgelist", 'nodelist', "seed", 'noise')

# Check for missing variables in the 'data_simulation' section
missing_vars <- required_vars_simulation[!required_vars_simulation %in% names(config)]
if (length(missing_vars) > 0) {
  stop(paste("Error: Missing required variables in 'data_simulation' section:", paste(missing_vars, collapse = ", ")))
}

# # Check if 'input_grn' exists and is a valid file path
# if (!file.exists(config$data_simulation$edgelist)) {
#   stop(paste("Error: Specified file", config$data_simulation$edgelist, "does not exist."))
# }
# # Check if 'input_grn' exists and is a valid file path
# if (!file.exists(config$data_simulation$nodelist)) {
#   stop(paste("Error: Specified file", config$nodelist, "does not exist."))
# }

read_dataframes<-function(filelist){
    list_edgelist<-list()
    counter<-1
    for(e in filelist){
        edgelist<-fread(file.path( e), sep='\t')
        edgelist$module<-counter
        list_edgelist[[counter]]<-edgelist
        counter<-counter+1
    }
    edgelist<-rbindlist(list_edgelist)

    return(edgelist)
}


output_dir<-file.path(opt$data_output_dir)



edgelist<-read_dataframes(config$edgelist )
nodelist<-read_dataframes(config$nodelist)

commonlist<-read_dataframes(config$common_edges )


print('UPDATED')
print(config$noise)
dataset<-create_gex_data_easy(net = edgelist, 
                         common_net = commonlist,
                         net_name = config$dataset_id,
                         seed = config$seed,
                         base_effect = config$base_effect,
                         mean = config$mean,
                         sd = config$sd,
                         noise=config$noise)

gex.dir<-file.path(output_dir)
gex.dir<-save_generated_data(dataset$net, dataset$counts, dataset$meta, gex_dir =gex.dir, disregulated_info = NULL)

print(gex.dir)


sce <- SingleCellExperiment(list(counts=dataset$counts), colData=dataset$meta, rowData = DataFrame(genes = rownames(dataset$counts)))

out_path <- file.path(gex.dir, 'data.h5ad')

sce <- sceasy::convertFormat(sce, from="sce", to="anndata",
                      outFile=out_path)