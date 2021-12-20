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
default <- list(path=FALSE, output="./combined", nodomains=FALSE, nodes=c())
opt <- setDefaults(config, default)

# Stop if no file provided
if (is.na(opt$pathA) | is.na(opt$pathB)) {
  stop("Please two files to be plotted.")
}

output <- opt$output
nodes <- opt$nodes

############## The Actual Code™ ###############
# Read in the matrices
mat1 <- read.csv(opt$pathA, header = TRUE)
mat2 <- read.csv(opt$pathB, header = TRUE)

# Set the first column as the row names
rownames(mat1) <- mat1[,1]
mat1 <- mat1[,2:ncol(mat1), drop=F]
rownames(mat2) <- mat2[,1]
mat2 <- mat2[,2:ncol(mat2), drop=F]

comparePaths(mat1, mat2, output = opt$output, nodes = nodes, nodomains = opt$nodomains)

cat("Done! Output written to", opt$output, "\n")
