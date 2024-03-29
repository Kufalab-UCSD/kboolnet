#############################################################
# plotPath.R
# Adrian C
#
# Function to plot attractor/path matrix generated by
# BoolNet's getPathToAttractor()
#
# NOTE: The attractor/path matrix needs to be first transposed
# using t() such that columns = timepoints and rows = nodes
#
# Args:
#   - path: Transposed BoolNet attractor/path matrix
#   - (filePath): Path to save PDF of graph, default: don't save
#   - (ratio): Width/length ratio of plot grid, default: 0.8
#
# Dependencies: ggplot2, dplyr, tidyr
############################################################

plotPath <- function(path, filePath = "", ratio = 0.8, attractor_idxs = NA) {
  if(!is.data.frame(path)) {
    stop("Path must be a data frame")
  }

  # If attractor indexes given, make them a different color
  if (any(!is.na(attractor_idxs))) {
    for (idx in attractor_idxs) {
      path[path[,idx] == 1, idx] <- 2
    }
  }

  # Add columnvector with symbols and proper column names
  path          <- as.data.frame(path)
  path_length   <- ncol(path)
  colnames(path) <- sprintf("t%03d", 0:(path_length-1)) # Replace t..0 with t000
  path$symbols  <- rownames(path)

  # Rearranges data frame
  pathGather          <- tidyr::gather(path, "t", "value", all_of(1:path_length))
  pathGather$symbols  <- factor(pathGather$symbols, levels = rownames(path)[nrow(path):1])
  pathGather$t        <- factor(pathGather$t)

  # Plot first simulation
  p         <- ggplot2::ggplot(pathGather, aes(t, symbols)) +
    ggplot2::geom_tile(aes(fill = value), colour = "white") +
    ggplot2::scale_fill_gradient2(limits = c(0,2), low = "white", mid = "steelblue", high = "purple", midpoint = 1)
  base_size <- 8
  p         <- p + ggplot2::theme_grey(base_size = base_size) + ggplot2::labs(x = "", y = "") +
    ggplot2::scale_x_discrete(expand = c(0, 0)) + ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::coord_fixed(ratio=ratio) + ggplot2::theme(legend.position = "none")

  # If file path provided, save plot to file
  if (filePath != "") {
    filePath <- gsub(".pdf$", "", filePath, ignore.case = TRUE) # Remove extension if present
    ggplot2::ggsave(paste0(filePath, ".pdf"), plot = last_plot(), height=(5 + length(levels(p$data$symbols)) * 0.25), width=(10 + length(levels(p$data$t)) * 0.45), scale = 1, units = "cm", limitsize=FALSE)
  }

  return(p)
}

