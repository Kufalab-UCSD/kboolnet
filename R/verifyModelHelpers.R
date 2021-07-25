plotCompMat <- function(mat) {
  # Plot comparison matrices (w/ ligand)
  df <- as.data.frame(mat)
  colnames(df) <- 1:ncol(df)
  rows <- nrow(df)
  df <- cbind(1:rows, df)
  colnames(df)[1] <- "y"
  df <- gather(df, "x", "value", 2:(rows+1), na.rm = T) # Rearrange data for plotting
  df <- transform(df, x = as.numeric(x))
  p <- ggplot(df, aes(x, y)) + geom_tile(aes(fill = value), color = "lightgrey") +
    scale_x_continuous(position = "top", expand = c(0,0), breaks=1:rows) + scale_y_reverse(expand=c(0,0), breaks=rows:1) +
    scale_fill_gradient2(limits=c(min(df$value)-0.0001,1), low = "red", mid="white", high="steelblue", midpoint=(1 + min(df$value) - 0.0001)/2) +
    theme_classic() + theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
                            axis.text = element_text(size = 10), panel.grid.minor = element_line(color = "lightgrey")) + coord_equal()
  return(p)
}

simulateWithLigands <- function(network, initStates, ligands) {
  # Add ligands in fully neutral state
  for (lig in 1:length(ligands)){
    initStates$state[grepl(paste0(ligands[lig], "_.*--0$"), initStates$name)] <- 1
    initStates$state[grepl(paste0(ligands[lig], "_.*-\\{0\\}$"), initStates$name)] <- 1
  }

  # Simulate and return
  return(getPathAndAttractor(network, initStates$state, initStates$name))
}

simulateWithoutLigands <- function(network, initStates, ligands) {
  # Unbind ligands from complexes and remove any ligand
  for (lig in 1:length(ligands)){
    initStates$state <- unbindLigand(ligands[lig], initStates$name, initStates$state)
    initStates$state[initStates$ID %in% mapLigToUnboundNodes(ligands, initStates)] <- 0
  }

  # Simulate and return
  return(getPathAndAttractor(network, initStates$state, initStates$name))
}


escapeRegex <- function(regex) {
  regex <- gsub('([.|()\\^{}+$?]|\\[|\\])', '\\\\\\1', regex)
  regex <- gsub('\\*', '\\.\\*\\?', regex)
  regex <- paste0("^", regex, "$")
  return(regex)
}

mapLigToUnboundNodes <- function(ligands, initStates) {
  ligNodes <- character()
  for (i in 1:length(ligands)) {
    # Make sure the ligand exists in a neutral state
    if (!(any(grepl(paste0("^", ligands[i], "(_.*--0|_.*-\\{0\\})$"), initStates$name)))) {
      stop("No neutral state found for ligand ", ligands[i], ". Please verify that ", ligands[i], " is a valid component in the rxncon system.")
    }

    # Make a regex matching unbound, and modified forms of ligand
    ligRegex <- paste0("(", c(paste0(ligands[i], "_.*--0"), paste0("^", ligands[i], "$"),
                              paste0(ligands[i], "_.*-\\{.*\\}")), ")", collapse="|")
    ligMatch <- grepl(ligRegex, initStates$name)
    ligNodes <- c(ligNodes, initStates$ID[ligMatch])

    cat("Ligand", ligands[i], "matched to node(s)", paste0(symbolMapping$name[ligMatch], collapse=", "), "\n")
  }

  return(ligNodes)
}

# Merge results into a single list, one entry per run
mergeResults <- function(results) {
  merged_results <- list()
  for (i in 1:length(results)) {
    run_name <- names(results)[i]
    merged_results[[run_name]] <- results[[run_name]][[1]]$noLig$path[,0]
    for (j in 1:length(results[[run_name]])) {
      x <- results[[run_name]][[j]]$noLig$path
      colnames(x) <- rep(paste0(j, "_", "noLig_path"), ncol(x))
      attr_idxs <- (ncol(x) + 1 - length(results[[run_name]][[j]]$noLig$attractor)):ncol(x)
      colnames(x)[attr_idxs] <- paste0(j, "_", "noLig_attractor")
      y <- results[[run_name]][[j]]$lig$path
      colnames(y) <- rep(paste0(j, "_", "lig_path"), ncol(y))
      colnames(y)[(ncol(y) + 1 - length(results[[run_name]][[j]]$lig$attractor)):ncol(y)] <- paste0(j, "_", "lig_attractor")
      merged_results[[run_name]] <- cbind(merged_results[[run_name]], x, y)
    }
  }

  return(merged_results)
}

identicalScores <- function(scores) {
  res <- list()
  visited <- numeric()
  for (i in 1:nrow(scores)) {
    identical <- numeric()
    for (j in i:ncol(scores)) {
      if (j %in% visited) {
        next
      } else if (scores[i, j] == 1) {
        identical <- c(identical, j)
        visited <- c(visited, j)
      }
    }
    if (length(identical) > 0) {
      res[[length(res) + 1]] <- identical
    }
  }
  return(res)
}

saveMergedResults <- function(merged_results, outPath) {
  for (i in 1:length(merged_results)) {
    write.csv(merged_results[[i]], file = paste0(outPath, "run_", names(merged_results)[i], ".csv"))
    attractor_idxs <- which(grepl("_attractor$", colnames(merged_results[[i]])))
    plotPath(merged_results[[i]], filePath = paste0(outPath, "run_", names(merged_results)[i], ".pdf"), attractor_idxs = attractor_idxs)
  }
}
