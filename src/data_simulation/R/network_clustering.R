
require(grn2gex)
require(data.table)
require(optparse)
library(yaml)
require(igraph)
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


# Check if the 'data_simulation' section exists
if (!"network_clustering" %in% names(config)) {
  stop("Error: Missing 'data_simulation' section in config file.")
}

###### DATA SIMULATION #######################################################################
# Define required variables within 'network_clustering'
required_vars_simulation <- c("input_network_file", 'output_network_dir')

# Check for missing variables in the 'network_clustering' section
missing_vars <- required_vars_simulation[!required_vars_simulation %in% names(config$network_clustering)]
if (length(missing_vars) > 0) {
  stop(paste("Error: Missing required variables in 'network_clustering' section:", paste(missing_vars, collapse = ", ")))
}




net.dir<-file.path('/usr/src/app/output', config$network_clustering$output_network_dir)
collectri.file <- file.path('/usr/src/app/input', config$network_clustering$input_network_file)


dir.create(net.dir, recursive = TRUE)

print(collectri.file)
collectri <- loadOrDownloadCollectTRI(collectri.file)
colnames(collectri)<- c('source', 'target', 'direction')
graph_list<-clusterNetwork(collectri)

print('Simplifying graph')
graph_list<-sapply(graph_list, function(g) simplify(g, remove.loops = TRUE))

graph_list<-remove_irrelevant_networks(graph_list)
orig_graph<-createGraph(collectri)
for (g in graph_list){
  g<-restore_original_directionality(g, orig_graph)
  save_small_subgraph(g$subgraph, net.dir, orig_graph)
}

