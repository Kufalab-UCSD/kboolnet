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

#################### Function definition ################

plotPath <- function(path, filePath = "", ratio = 0.8) {
  path <- as.data.frame(path)
  path_length   <- ncol(path)
  colnames(path) <- sprintf("t%03d", 0:(path_length-1)) # Replace t..0 with t000
  path$symbols  <- rownames(path)
  
  # Rearranges data frame
  pathGather          <- gather(path, "t", "value", 1:path_length)
  pathGather$symbols  <- factor(pathGather$symbols, levels = rownames(path)[nrow(path):1])
  pathGather$t        <- factor(pathGather$t)
  pathGather$value    <- factor(pathGather$value, levels = c("0", "1", "2", "3", "4"))
  
  # Plot first simulation
  compColors <- c("0" = "white", "1" = "red", "2" = "steelblue", "3" = "purple", "4" = "lightgrey")
  p         <- ggplot(pathGather, aes(t, symbols)) + geom_tile(aes(fill = value), colour = "white") +
                      scale_fill_manual(values = compColors, labels=c("Neither", "1", "2", "1+2", "NA"), drop=FALSE)
  base_size <- 8
  ratio <- 0.8
  p         <- p + theme_grey(base_size = base_size) + labs(x = "", y = "") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
                coord_fixed(ratio=ratio) + guides(fill = guide_legend(title="Active in path:"))
  
  # Save plot to file
  if (filePath != "") {
    filePath <- gsub(".pdf$", "", filePath, ignore.case = TRUE) # Remove extension if present
    ggsave(paste0(filePath, ".pdf"), plot = last_plot(), height=(5 + length(levels(p$data$symbols)) * 0.25), width=(10 + length(levels(p$data$t)) * 0.45), scale = 1, units = "cm")
  }
  
  return(p)
}

############# Argument parsing and handling ############

option_list = list(
  make_option(c("-p", "--path"), action="store_true", default=FALSE,
              help="Indicate that input files are paths, not attractors. This slightly changes how alignment is performed."),
  make_option(c("-o", "--output"), action="store", default="./combined",
              help="Base name for output files. [default %]")
)
usage <- "usage: %prog [options] FILE1 FILE2\n\n FILE1 and FILE2 are the .csv files to be compared"
opt <- parse_args(OptionParser(usage=usage, option_list=option_list), positional_arguments = 2)

output <- opt$options$output

############## The Actual Code™ ###############

# Read in the matrices
mat1 <- read.csv(opt$args[1], header = TRUE)
mat2 <- read.csv(opt$args[2], header = TRUE)

# Set the first column as the row names
rownames(mat1) <- mat1[,1]
mat1 <- mat1[,2:ncol(mat1)]
rownames(mat2) <- mat2[,1]
mat2 <- mat2[,2:ncol(mat2)]

if (any(rownames(mat1) != rownames(mat2))) {
  stop("Path matrices must be in same order and have same number of rows.")
}
numRows <- nrow(mat1)

mats <- list(as.matrix(mat1), as.matrix(mat2))

# # If there is a mismatch in number of columns, add empty columns to the shorter matrix
# if (ncol(mat1) < ncol(mat2)) {
#   mat1 <- cbind(mat1, matrix(nrow=nrow(mat1), ncol=ncol(mat2)-ncol(mat1)))
# } else if (ncol(mat2) < ncol(mat1)) {
#   mat2 <- cbind(mat2, matrix(nrow=nrow(mat2), ncol=ncol(mat1)-ncol(mat2)))
# }
# 
# # Set NAs to equal 0
# mat1[is.na(mat1)] <- 0
# mat2[is.na(mat2)] <- 0

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
  if (opt$options$path) break
}

# Apply highest-scoring order to long matrix
bestOrder <- orders[[which.max(scores)]]
mats[[longMat]] <- mats[[longMat]][,bestOrder]

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

# Multiply mat2 by 2 and then add the matrices. As a result, values will be as follows:
# mat1 and mat2 on: 3, only mat2 on: 2, only mat1 on: 1, neither on: 0
mat_combined_all <- mats[[1]] + (mats[[2]] * 2)

# Get rows with not matching values
mat_combined_diff <- mat_combined_all[apply(mat_combined_all, 1, function(x) {
  any(sapply(x, function(y) y == 1 | y == 2))
}),]


############### Plotting ################
# Plot entire path
plotPath(mat_combined_all, paste0(output, "_all.pdf"))
  
# Plot only different rows
plotPath(mat_combined_diff, paste0(output, "_diff.pdf"))
