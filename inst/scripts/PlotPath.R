#!/usr/bin/env Rscript

##############################################################
# PlotPath.R
# Adrian C
#
# Command-line wrapper for plotPath.R function
#
##############################################################

############## Library loading ###############################
options(stringsAsFactors = F)
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(optparse))
library(kboolnet)

################## Argument loading/parsing ################
# If not interactive, get config file
if (!interactive()) {
  args = commandArgs(trailingOnly = TRUE)
  if (length(args) != 1) {
    stop("Please provide path to config file as the only argument to the script.")
  } else if (!file.exists(normalizePath(args))) {
    stop("File ", args, " does not exist.")
  }

  source(normalizePath(args))
}

# Check that config is loaded, print it
if (!exists("config")) {
  stop("Config file must be loaded in before running the script.")
}
print(config)

# Set default args if they are not already set
default <- list(out="./plot.pdf", nodes=c(), nodomains=FALSE, ratio=0.8)
opt <- setDefaults(config, default)

# Stop if no file provided
if (is.na(opt$file)) {
  stop("Please provide the file to be plotted.")
}

# Process the nodes argument into a list
nodes <- opt$nodes

################# Path loading/processing/plotting ##################
# Load path and set first col as row names
path <- read.csv(opt$file, header=TRUE)
rownames(path) <- path[,1]
path <- path[,2:ncol(path), drop=F]

# Keep only the nodes that are wanted
if (length(nodes) > 0) {
  # First, check that every nodes argument actually exists in the path
  if (any(!(nodes %in% rownames(path))) ){
    stop("Node(s) ", paste0(nodes[!(nodes %in% rownames(path))], collapse=", "), " are not in the path/attractor file")
  }

  # Get indices of nodes to keep and only keep them
  path <- path[nodes, ,drop=F]
}

# Remove domain names if requested
if (opt$nodomains) {
  rownames(path) <- removeDomains(rownames(path))
}

plotPath(path, opt$out, ratio=opt$ratio)
cat("Done! Output written to", opt$path, "\n")
