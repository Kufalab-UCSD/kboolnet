#!/usr/bin/env Rscript
#############################################
# plotMatrixComp.R
# Adrian C
#
# Takes two binary matrices, aligns them as best
# as possible, and then creates a plot showing
# the differences between the two.
############################################

suppressMessages(library(ggplot2))
library(tidyr)

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
  
  # Plot first simulation
  compColors <- c("0" = "white", "1" = "red", "2" = "steelblue", "3" = "purple", "4" = "lightgrey")
  p         <- ggplot(pathGather, aes(t, symbols)) + geom_tile(aes(fill = factor(value)), colour = "white") + scale_fill_manual(values = compColors, labels =c("Neither", "1", "2", "1+2", "NA") )
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

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Please provide two .csv files as arguments to this function.")
}

############## The Actual Code™ ###############

# Read in the matrices
mat1 <- read.csv(args[1], header = TRUE)
mat2 <- read.csv(args[2], header = TRUE)

# Set the first column as the row names
rownames(mat1) <- mat1[,1]
mat1 <- mat1[,2:ncol(mat1)]
rownames(mat2) <- mat2[,1]
mat2 <- mat2[,2:ncol(mat2)]

if (any(rownames(mat1) != rownames(mat2))) {
  stop("Path matrices must be in same order.")
}

mat1 <- as.matrix(mat1)
mat2 <- as.matrix(mat2)

# If there is a mismatch in number of columns, add empty columns to the shorter matrix
if (ncol(mat1) < ncol(mat2)) {
  mat1 <- cbind(mat1, matrix(nrow=nrow(mat1), ncol=ncol(mat2)-ncol(mat1)))
} else if (ncol(mat2) < ncol(mat1)) {
  mat2 <- cbind(mat2, matrix(nrow=nrow(mat2), ncol=ncol(mat1)-ncol(mat2)))
}

# Set NAs to equal 0
mat1[is.na(mat1)] <- 0
mat2[is.na(mat2)] <- 0

# Align and score the matrices
scores <- numeric()
orders <- list()
orders[1] <- list(1:ncol(mat1))
for (i in 1:ncol(mat1)) {
  # Create new ordering
  if (i != 1) {
    orders[[i]] <- c((ncol(mat1) - (i - 2)):ncol(mat1), 1:(ncol(mat1) - (i - 1)))
  }
  
  # Apply new ordering to matrix 2 and calculate the similarity (matching/total)
  vec1 <- as.vector(mat1)
  vec2 <- as.vector(mat2[,orders[[i]]])
  scores[i] <- sum(vec1 == vec2, na.rm=T)/length(vec1)
  
  # If perfect score, we can stop there
  if (scores[i] == 1) break
}

# Apply highest-scoring order to mat2
mat2 <- mat2[,orders[[which.max(scores)]]]

# Multiply mat2 by 2 and then add the matrices. As a result, values will be as follows:
# mat1 and mat2 on: 3, only mat2 on: 2, only mat1 on: 1, neither on: 0
mat_combined_all <- mat1 + (mat2 * 2)

# Get rows with not matching values
mat_combined_diff <- mat_combined_all[apply(mat_combined_all, 1, function(x) {
  any(sapply(x, function(y) y == 1 | y ==2))
}),]


############### Plotting ################
# Plot entire path
plotPath(mat_combined_all, "combined_all.pdf")
  
# Plot only different rows
plotPath(mat_combined_diff, "combined_diff.pdf")