pathToGraph <- function(path, attractor = FALSE) {
  # Covert path to a vector of strings
  pathVec <- apply(path, 2, paste0, collapse="")

  # If there is no loop, make last state loop on itself
  if (!(pathVec[length(pathVec)] %in% pathVec[-length(pathVec)])) {
    pathVec <- c(pathVec, pathVec[length(pathVec)])
  }

  # Convert vector to an edge list
  df <- data.frame(source=numeric(), target=numeric())
  for (i in 1:(length(pathVec) - 1)) {
    df[nrow(df) + 1,] <- c(pathVec[i], pathVec[i + 1])
  }

  # Turn into a graph
  g <- graph_from_data_frame(df)

  return(g)
}

makeStateSpaceFromResults <- function(results) {
  stateSpaces <- list()
  # First add all the nodes to the graph
  for (i in 1:length(results)) {
    stateSpaces[[i]] <- graph.empty()

    # First make the base graphs, labeling attractor nodes as such
    attractor_nodes <- c()
    init_state <- NA
    for (j in 1:length(results[[i]])) {
      if (j == 1) {
        init_state <- paste0(results[[i]][[j]]$noLig$path[,1], collapse = "")
      }
      stateSpaces[[i]] <- igraph::union(stateSpaces[[i]], pathToGraph(cbind(results[[i]][[j]]$noLig$path, results[[i]][[j]]$noLig$attractor)))
      stateSpaces[[i]] <- igraph::union(stateSpaces[[i]], pathToGraph(cbind(results[[i]][[j]]$lig$path, results[[i]][[j]]$lig$attractor)))

      attractor_nodes <- union(attractor_nodes, apply(results[[i]][[j]]$noLig$attractor, 2, paste0, collapse=""))
      attractor_nodes <- union(attractor_nodes, apply(results[[i]][[j]]$lig$attractor, 2, paste0, collapse=""))
    }
    V(stateSpaces[[i]])$type <- "trajectory"
    V(stateSpaces[[i]])[attractor_nodes]$type <- "attractor"
    V(stateSpaces[[i]])[init_state]$type <- "initial_state"

    # Connect no ligand attractors to w/ ligand paths
    for (j in 1:length(results[[i]])) {
      noLigEnd <- paste0(results[[i]][[j]]$noLig$path[,ncol(results[[i]][[j]]$noLig$path)], collapse = "")
      ligStart <- paste0(results[[i]][[j]]$lig$path[,1], collapse = "")
      stateSpaces[[i]] <- stateSpaces[[i]] + edge(c(noLigEnd, ligStart), action = "ligand_addition")
    }

    # Then connect w/ ligand attractors to no ligand paths
    if (length(results[[i]]) > 1) {
      for (j in 1:(length(results[[i]]) - 1)) {
        ligEnd <- paste0(results[[i]][[j]]$lig$path[,ncol(results[[i]][[j]]$lig$path)], collapse = "")
        noLigStart <- paste0(results[[i]][[j + 1]]$noLig$path[,1], collapse = "")
        stateSpaces[[i]] <- stateSpaces[[i]] + edge(c(ligEnd, noLigStart), action = "ligand_removal")
      }
    }

    stateSpaces[[i]] <- simplify(stateSpaces[[i]], remove.multiple = FALSE)
  }
  names(stateSpaces) <- names(results)

  return(stateSpaces)
}
