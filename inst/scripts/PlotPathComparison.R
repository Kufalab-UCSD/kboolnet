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
  opt <- loadConfig(opt)
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

comparePaths(mat1, mat2, output = opt$output, nodes = nodes, nodomains = opt$nodomains)
