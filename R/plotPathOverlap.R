alignMats <- function(mat1, mat2) {
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
  orders[1] <- list(1:ncol(mats[[longMat]]))
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

    # # If this is a path, only do this with original ordering
    # if (opt$path) break
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
  mats[[shortMat]] <- tmp

  return(mats)
}

plotPathOverlap <- function(path1, path2, filePath = "", ratio = 0.8) {
  path_length   <- ncol(path1)
  colnames(path1) <- sprintf("t%03d", 0:(path_length-1)) # Replace t..0 with t000
  colnames(path2) <- sprintf("t%03d", 0:(path_length-1)) # Replace t..0 with t000
  path1$symbols  <- rownames(path1)
  path2$symbols  <- rownames(path2)

  # Rearranges data frame 1
  pathGather1         <- gather(path1, "t", "value", 1:path_length)
  pathGather1$symbols <- factor(pathGather1$symbols, levels = rownames(path1)[nrow(path1):1])
  pathGather1$t       <- factor(pathGather1$t)
  pathGather1$value   <- factor(pathGather1$value, levels = c("0", "1"))

  # Rearranges data frame 2
  pathGather2         <- gather(path2, "t", "value", 1:path_length)
  pathGather2$symbols <- factor(pathGather2$symbols, levels = rownames(path2)[nrow(path2):1])
  pathGather2$t       <- factor(pathGather2$t)
  pathGather2$value   <- factor(pathGather2$value, levels = c("0", "1", "2"))
  pathGather2$value[pathGather2$value == "1"] <- "2"

  # Plot first simulation
  compColors <- c("0" = "transparent", "1" = "red", "2" = "steelblue")
  p          <- ggplot() +
                  geom_tile(data = pathGather1, mapping = aes(t, symbols, fill = value), colour = "white", alpha = 0.7) +
                  geom_tile(data = pathGather2, mapping = aes(t, symbols, fill = value), colour = "white", alpha = 0.5) +
                  scale_fill_manual(values = compColors, labels=c("Neither", "1", "2", "1+2", "NA"), drop=FALSE)
  base_size <- 8
  ratio <- 0.8
  p         <- p + theme_grey(base_size = base_size) + theme(panel.background = element_blank()) + labs(x = "", y = "") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
                coord_fixed(ratio=ratio) + guides(fill = guide_legend(title="Active in path:"))

  # Save plot to file
  if (filePath != "") {
    filePath <- gsub(".pdf$", "", filePath, ignore.case = TRUE) # Remove extension if present
    ggsave(paste0(filePath, ".pdf"), plot = last_plot(), height=(5 + length(levels(pathGather1$symbols)) * 0.25), width=(10 + length(levels(pathGather1$t)) * 0.45), scale = 1, units = "cm")
  }

  return(p)
}

comparePaths <- function(path1, path2, output = "./combined", nodes = c(), nodomains = FALSE) {
  mats <- alignMats(path1, path2)
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
  if (nodomains) {
    rownames(mats[[1]]) <- removeDomains(rownames(mats[[1]]))
    rownames(mats[[2]]) <- removeDomains(rownames(mats[[2]]))
  }

  # Get rows with not matching values
  diffRows <- rowSums(abs(mats[[1]] - mats[[2]])) != 0

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
}
