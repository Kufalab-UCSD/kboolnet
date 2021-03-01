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
default <- list(path=FALSE, output="./combined", nodomains=FALSE, nodes="")
default <- default[!(names(default) %in% names(opt))]
opt     <- c(opt, default)

# Stop if no file provided
if (is.na(opt$pathA) | is.na(opt$pathB)) {
  stop("Please two files to be plotted.")
}

output <- opt$output
print(opt$nodes)
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

if (any(rownames(mat1) != rownames(mat2))) {
  stop("Path matrices must be in same order and have same number of rows.")
}
numRows <- nrow(mat1)

mats <- list(as.matrix(mat1, nrow=nrow(mat1), ncol=ncol(mat1)), as.matrix(mat2, nrow=nrow(mat2), ncol=ncol(mat2)))

# Figure out which matrix is long and which one is short
longMat <- which.max(lapply(mats, length))
if (longMat == 1) shortMat = 2
if (longMat == 2) shortMat = 1

# Align and score the matrices
scores <- numeric()
shifts <- integer()
orders <- list()
orders[1] <- list(1:ncol(mat1))
for (i in 1:ncol(mats[[longMat]])) {
  # Create new ordering
  if (i != 1) {
    orders[[i]] <- c((ncol(mats[[longMat]]) - (i - 2)):ncol(mats[[longMat]]), 1:(ncol(mats[[longMat]]) - (i - 1)))
  }

  # Apply new ordering to long matrix and convert to vector
  longVec <- as.vector(mats[[longMat]][,orders[[i]]])

  # Calculate scores for all horizontal shifts for short matrix
  shiftScores <- numeric()
  for (j in 0:(ncol(mats[[longMat]]) - ncol(mats[[shortMat]]))) {
    tmp <- mats[[shortMat]]
    # Add columns to the left of short matrix if necessary
    if (j > 0) {
      tmp <- cbind(matrix(nrow=numRows, ncol=j), tmp)
    }
    # Add columns to the right of short matrix if necessary
    if (j < (ncol(mats[[longMat]]) - ncol(mats[[shortMat]]))) {
      tmp <- cbind(tmp, matrix(nrow=numRows, ncol=((ncol(mats[[longMat]]) - ncol(mats[[shortMat]])) - j)))
    }

    # Calculate score
    shortVec <- as.vector(tmp)
    shiftScores[j+1] <- sum(longVec == shortVec, na.rm=TRUE) / length(longVec)
  }

  # Find highest scoring horizontal "shift" for short matrix and keep that score
  shifts[i] <- which.max(shiftScores) - 1
  scores[i] <- max(shiftScores)

  # If perfect score, we can stop there
  if (scores[i] == 1) break

  # If this is a path, only do this with original ordering
  if (opt$path) break
}

# Apply highest-scoring order to long matrix
bestOrder <- orders[[which.max(scores)]]
mats[[longMat]] <- mats[[longMat]][,bestOrder, drop=F]

# Apply highest-scoring shift to short matrix
bestShift <- shifts[which.max(scores)]
tmp <- mats[[shortMat]]
if (bestShift > 0) { # Add to right
  tmp <- cbind(matrix(nrow=numRows, ncol=bestShift), tmp)
}
if (bestShift < (ncol(mats[[longMat]]) - ncol(mats[[shortMat]]))) { # Add to left
  tmp <- cbind(tmp, matrix(nrow=numRows, ncol=((ncol(mats[[longMat]]) - ncol(mats[[shortMat]])) - bestShift)))
}
tmp[is.na(tmp)] <- 0
mats[[shortMat]] <- tmp

# Order rows according to nodes argument
if (length(nodes) != 0) {
  mats[[1]] <- mats[[1]][nodes,,drop=F]
  mats[[2]] <- mats[[2]][nodes,,drop=F]
}

# Remove domain information if requested
if (opt$nodomains) {
  newNames <- gsub("_\\[.*?\\]", "", rownames(mats[[1]]))

  # If there are ambigious names due to domain simplification, replace them with the old names
  ambigNames <- duplicated(newNames) | duplicated(newNames, fromLast = T)
  newNames[ambigNames] <- rownames(mats[[1]])[ambigNames]

  rownames(mats[[1]]) <- newNames
  rownames(mats[[2]]) <- newNames
}

# Get rows with not matching values
diffRows <- rowSums(mats[[1]] - mats[[2]]) != 0

############### Plotting ################
# Plot entire path
plotPathOverlap(mats[[1]], mats[[2]], paste0(output, "_all.pdf"))

# Plot only different rows
if (sum(diffRows) == 0) {
  cat("No difference between attractors for selected nodes! Skipping difference plot. \n")
} else {
  plotPathOverlap(mats[[1]][diffRows,,drop=F], mats[[2]][diffRows,,drop=F], paste0(output, "_diff.pdf"))
}
