#!/usr/bin/env Rscript
#############################################
# plotMatrixComp.R
# Adrian C
#
# Takes two binary matrices, aligns them as best
# as possible, and then creates a plot showing
# the differences between the two.
############################################

#################### Library loading ################

suppressMessages(library(ggplot2))
library(tidyr)
library(optparse)
library(kboolnet)

############# Argument parsing and handling ############

option_list = list(
  make_option(c("--pathA", "-a"), action="store", default=NA, type="character",
              help="First of the two paths to be compared."),
  make_option(c("--pathB", "-b"), action="store", default=NA, type="character",
              help="Second of the two paths to be compared."),
  make_option(c("--config", "-c"), action="store", default=NA, type="character",
              help="Path of config file. You can specify parameters here instead of passing them as command-line
              arguments"),
  make_option(c("-p", "--path"), action="store_true", default=NA,
              help="Indicate that input files are paths, not attractors. This slightly changes how alignment is performed."),
  make_option(c("-o", "--output"), action="store", default=NA,
              help="Base name for output files. [default: ./combined]"),
  make_option(c("-n", "--nodes"), action="store", default=NA,
              help="Comma-separated ordered list of nodes to plot. [default: all]"),
  make_option(c("--nodomains", "-d"), action="store_true", default=NA, type="logical",
              help="Remove domains from rxncon node names. [default: don't remove]")
)
opt <- parse_args(OptionParser(option_list=option_list))
opt <- opt[!is.na(opt)] # Discard NA values

# Load config file if provided
if ("config" %in% names(opt)) {
  opt <- loadConfig(opt, config)
}

# Set default args if they are not already set
default <- list(path=FALSE, output="./combined", nodomains=FALSE, nodes="")
opt <- setDefaults(opt, default)

# Stop if no file provided
if (is.na(opt$pathA) | is.na(opt$pathB)) {
  stop("Please two files to be plotted.")
}

output <- opt$output
nodes <- trimws(strsplit(opt$nodes, ",")[[1]])

############## The Actual Code™ ###############

# Read in the matrices
mat1 <- read.csv(opt$pathA, header = TRUE)
mat2 <- read.csv(opt$pathB, header = TRUE)

# Set the first column as the row names
rownames(mat1) <- mat1[,1]
mat1 <- mat1[,2:ncol(mat1), drop=F]
rownames(mat2) <- mat2[,1]
mat2 <- mat2[,2:ncol(mat2), drop=F]

mats <- alignMats(mat1, mat2)
mats[[1]][is.na(mats[[1]])] <- 0
mats[[2]][is.na(mats[[1]])] <- 0

# Order rows according to nodes argument
if (length(nodes) != 0) {
  missing <- nodes[!nodes %in% rownames(mats[[1]])]
  if (length(missing) > 0) {
    stop(paste0("Node(s) ", paste(missing, collapse = ", "), " not found in attractor. Make sure to reference full names, including domains."))
  }
  mats[[1]] <- mats[[1]][nodes,,drop=F]
  mats[[2]] <- mats[[2]][nodes,,drop=F]
}

# Remove domain information if requested
if (opt$nodomains) {
  rownames(mats[[1]]) <- removeDomains(rownames(mats[[1]]))
  rownames(mats[[2]]) <- removeDomains(rownames(mats[[2]]))
}

# Get rows with not matching values
diffRows <- rowSums(mats[[1]] - mats[[2]]) != 0

# Covert to data frame
pathRowNames <- rownames(mats[[1]])
mats[[1]] <- as.data.frame(mats[[1]])
mats[[2]] <- as.data.frame(mats[[2]])
rownames(mats[[1]]) <- pathRowNames
rownames(mats[[2]]) <- pathRowNames


############### Plotting ################
# Plot entire path
plotPathOverlap(mats[[1]], mats[[2]], paste0(output, "_all.pdf"))

# Plot only different rows
if (sum(diffRows) == 0) {
  cat("No difference between attractors for selected nodes! Skipping difference plot. \n")
} else {
  plotPathOverlap(mats[[1]][diffRows,,drop=F], mats[[2]][diffRows,,drop=F], paste0(output, "_diff.pdf"))
}
