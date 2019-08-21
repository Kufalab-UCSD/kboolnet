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

############### Argument parsing and setup ###################
# Get commandline args
option_list = list(
  make_option("--config", action="store", default=NA, type="character",
              help="Path of config file. You can specify parameters here instead of passing them as command-line
              arguments"),
  make_option("--kboolnetPath", action="store", default=NA, type="character",
              help="Path to root directory of kboolnet repository"),
  make_option(c("--file", "-f"), action="store", default=NA, type="character",
              help="Name of csv file containing path/attractor to be plotted"),
  make_option(c("--nodes", "-n"), action="store", default=NA, type="character",
              help="Comma-separated list of nodes to be displayed in plot. [default: all nodes]"),
  make_option(c("--out", "-o"), action="store", default=NA, type="character",
              help="Name of PDF file to which plot should be written. [default: plot.pdf]"),
  make_option(c("--nodomains", "-d"), action="store_true", default=NA, type="logical",
              help="Remove domains from rxncon node names. [default: don't remove]")
)
opt <- parse_args(OptionParser(option_list=option_list))
opt <- opt[!is.na(opt)] # Discard NA values

# Load config file if provided
if ("config" %in% names(opt)) {
  source(opt$config)

  if (!exists("config")) {
    stop("No config object found in config file")
  }

  config <- config[!is.na(config)] # Discard NA values

  # Keep only config values that were not passed as command line options
  config <- config[!(names(config) %in% names(opt))]

  # Merge command line args and config values
  opt <- c(opt, config)
}

# Set default args if they are not already set
default <- list(out="./plot.pdf", nodes="", nodomains=FALSE)
default <- default[!(names(default) %in% names(opt))]
opt     <- c(opt, default)

# Stop if no file provided
if (is.na(opt$file)) {
  stop("Please provide the file to be plotted.")
}

# Process the nodes argument into a list
nodes <- trimws(strsplit(opt$nodes, ",")[[1]])

# Normalize paths
kboolnetPath  <- paste0(normalizePath(opt$kboolnetPath), "/")

# Load plot function
suppressMessages(source(paste0(kboolnetPath, "R/functions/plotPath.R")))

################# Path loading/processing/plotting ##################
# Load path and set first col as row names
path <- read.csv(opt$file, header=TRUE)
rownames(path) <- path[,1]
path <- path[,2:ncol(path), drop=F]

# Remove domain names if requested
if (opt$nodomains) {
  newNames <- gsub("_\\[.*?\\]", "", rownames(path))
  
  # If there are ambigious names due to domain simplification, replace them with the old names
  ambigNames <- duplicated(newNames) | duplicated(newNames, fromLast = T)
  newNames[ambigNames] <- rownames(path)[ambigNames]
  
  rownames(path) <- newNames
}

# Keep only the nodes that are wanted
if (length(nodes > 0)) {
  # First, check that every nodes argument actually exists in the path
  if (any(!(nodes %in% rownames(path))) ){
    stop("Node(s) ", paste0(nodes[!(nodes %in% rownames(path))], collapse=", "), " are not in the path/attractor file")
  }
  
  # Get indices of nodes to keep and only keep them
  path <- path[nodes, ,drop=F]
}

plotPath(path, opt$out)
