#!/usr/bin/env Rscript

#################################################
# VerifyModel.R
# Adrian C
#
# Script to download model from Google Drive
# and run a "sanity check" round of simulations.
#################################################

################# Library loading ##################
options(stringsAsFactors = F)
suppressMessages(library(BoolNet))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(openxlsx))
suppressMessages(library(googledrive))
suppressMessages(library(optparse))
suppressMessages(library(tidyr))
suppressMessages(library(numbers))
suppressMessages(library(igraph))
library(kboolnet)

############## Argument parsing ####################
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
default <- list(modules="", out="./out/", minQuality=0, ligands=NA, file=NA, driveFile=NA, rounds=20, inhib="")
opt <- setDefaults(config, default)

# Create out dir if it does not exist
if (!dir.exists(opt$out)) {
  dir.create(opt$out)
}

# Normalize paths
outPath <- paste0(normalizePath(opt$out), "/")

################ Load and process rxncon file ###################

# If GDrive file provided
if (!(is.na(opt$driveFile))) {
  cat("Downloading rxncon file from Google Drive...", "\n")

  # Download file
  masterFile  <- paste0(outPath, "master.xlsx")
  driveDownload(driveFile = opt$driveFile, out = masterFile, type = "spreadsheet")

# If local file provided
} else {
  masterFile <- normalizePath(opt$file)

  # File verification
  if (!grepl("\\.xlsx$", masterFile)) { # Make sure file is Excel file
    stop("rxncon file must be an Excel file (.xslx extension)")
  } else if (!file.exists(masterFile)) { # Make sure file exists
    stop("rxncon file does not exist")
  }
}

# Extract modules from master file, write to modules file
modulesFile <- paste0(outPath, "modules.xlsx")
cat("Extracting modules... ")
callExtractModules(masterFile, modulesFile, opt$modules)
cat("Done.\n")

# Pass files to rxncon for processing
cat("Generating BoolNet files... ")
netFilePrefix <- gsub("\\.xlsx$", "", modulesFile)
callRxncon2Boolnet(modulesFile, netFilePrefix)
cat("Done.\n")

################# Load BoolNet files ######################
# Load network
network    <- loadNetwork(paste0(netFilePrefix, ".boolnet"), symbolic=TRUE)

# Load symbols and remove spaces
symbolMapping       <- read.csv(paste0(netFilePrefix, "_symbols.csv"), header=F, col.names=c("ID", "name"))
symbolMapping$name  <- gsub("[[:space:]]", "", symbolMapping$name)
symbolMapping       <- symbolMapping[,1:2]

# Load initial states from file
initStates <- read.csv(paste0(netFilePrefix, '_initial_vals.csv'), col.names=c("ID","state","name"), header=F)
initStates$name <- gsub("# ", "", initStates$name) # Clean up names
initStates$name <- gsub(" ", "", initStates$name)

############## Simulation time #####################
results <- list()
scores <- list()
for (i in 1:length(opt$runs)) {
  run <- opt$runs[[i]]
  runNet <- network
  run_res <- list(list())
  scoreLig <- matrix(nrow = opt$rounds, ncol = opt$rounds)
  scoreNoLig <- matrix(nrow = opt$rounds, ncol = opt$rounds)
  scoreLig[1,1] <- 1
  scoreNoLig[1,1] <- 1

  # Load initial state if given
  run_initStates <- initStates
  if ("initial" %in% names(run)) {
    # TODO: proper input checking
    run_initStates <- read.csv(run$initial, col.names=c("ID","state","name"), header=F)
    run_initStates$name <- gsub("# ", "", run_initStates$name) # Clean up names
    run_initStates$name <- gsub(" ", "", run_initStates$name)
  }

  # Apply on/off nodes
  if (length(run$on) > 0) {
    runNet <- fixedNetwork(runNet, run$on, initStates$name, initStates$ID, 1, regex=TRUE)
  }
  if (length(run$off) > 0) {
    runNet <- fixedNetwork(runNet, run$off, initStates$name, initStates$ID, 0, regex=TRUE)
  }

  # # Equilibrate first
  # tmp <- simulateWithLigands(runNet, run_initStates, run$toggle)
  # run_initStates$state <- tmp$attractor[,1]

  # First round of simulation
  run_res[[1]]$noLig <- simulateWithoutLigands(runNet, run_initStates, run$toggle)
  run_initStates$state <- run_res[[1]]$noLig$attractor[,1]
  run_res[[1]]$lig <- simulateWithLigands(runNet, run_initStates, run$toggle)
  run_initStates$state <- run_res[[1]]$lig$attractor[,1]

  # Following rounds of simulation
  for (j in 2:opt$rounds) {
    run_res[[j]] <- list()

    ####### NO LIG ######
    # Simulate using every possible initial state from the previous w/ ligand attractor
    noLigRes <- list()
    for (k in 1:ncol(run_res[[j - 1]]$lig$attractor)) {
      run_initStates$state <- run_res[[j - 1]]$lig$attractor[,k]
      noLigRes[[k]] <- simulateWithoutLigands(runNet, run_initStates, run$toggle)
    }

    # Compare these to each other
    consistencyScores <- matrix(nrow = length(noLigRes), ncol = length(noLigRes))
    for (k in 1:length(noLigRes)) {
      for (z in k:length(noLigRes)) {
        consistencyScores[k,z] <- 1 - attractorDistance(noLigRes[[k]]$attractor, noLigRes[[z]]$attractor)
      }
    }

    # Check if identical. If not, dump results and error out
    identicalResults <- identicalScores(consistencyScores)
    if (length(identicalResults) > 1) {
      # Save full results, merge them, and save that too
      merged_results <- mergeResults(results)
      saveMergedResults(merged_results, outPath)

      # Make state space from the results
      results[[names(opt$runs)[i]]] <- run_res[1:j-1]
      stateSpaces <- makeStateSpaceFromResults(results)

      # Add the inconsistent attractors to it
      attractor_nodes <- c()
      for (k in 1:length(noLigRes)) {
        stateSpaces[[i]] <- igraph::union(stateSpaces[[i]], pathToGraph(cbind(noLigRes[[k]]$path, noLigRes[[k]]$attractor)))
        attractor_nodes <- union(attractor_nodes, apply(noLigRes[[k]]$attractor, 2, paste0, collapse=""))

        ligEnd <- paste0(run_res[[j - 1]]$lig$attractor[,k], collapse = "")
        attractor_nodes <- union(attractor_nodes, ligEnd)
        noLigStart <- paste0(noLigRes[[k]]$path[,1], collapse = "")
        stateSpaces[[i]] <- stateSpaces[[i]] + edge(c(ligEnd, noLigStart), action = "ligand_removal")
      }
      V(stateSpaces[[i]])[is.na(V(stateSpaces[[i]])$type)]$type <- "trajectory"
      V(stateSpaces[[i]])[attractor_nodes]$type <- "attractor"

      # Save state spaces
      for (i in 1:length(stateSpaces)) {
        stateSpaces[[i]] <- simplify(stateSpaces[[i]], remove.multiple = FALSE)
        write_graph(stateSpaces[[i]], file = paste0(outPath, "run_", names(stateSpaces)[i], "_state_space.graphml"), format = "graphml")
      }

      # Save initial states that led to inconsistent attractors, as well as full paths
      for (k in 1:length(identicalResults)) {
        idx <- identicalResults[[k]][1]
        run_initStates$state <- noLigRes[[idx]]$path[,1]
        write.csv(run_initStates, file = paste0(outPath, "inconsistent_initState_", k, ".csv"))
      }

      # Compare the first two of the inconsistent attractors
      comparePaths(noLigRes[[identicalResults[[1]][1]]]$attractor, noLigRes[[identicalResults[[2]][1]]]$attractor,
                   output = paste0(outPath, "inconsistent_attractor"))


      save(results, file = paste0(outPath, "results.RData"))
      stop("Inconsistent no-ligand attractors were found on round ", j, " of run ", names(opt$runs)[i],
           ". Dumping results as well as the initial states which resulted in inconsistent attractors.")
    }

    # Save result
    run_res[[j]]$noLig <- noLigRes[[1]]

    ###### WITH LIG ######
    # Simulate using every possible initial state from the previous no ligand attractor
    ligRes <- list()
    for (k in 1:ncol(run_res[[j]]$noLig$attractor)) {
      run_initStates$state <- run_res[[j]]$noLig$attractor[,k]
      ligRes[[k]] <- simulateWithLigands(runNet, run_initStates, run$toggle)
    }

    # Compare these to each other
    consistencyScores <- matrix(nrow = length(noLigRes), ncol = length(noLigRes))
    for (k in 1:length(noLigRes)) {
      for (z in k:length(noLigRes)) {
        consistencyScores[k,z] <- 1 - attractorDistance(noLigRes[[k]]$attractor, noLigRes[[z]]$attractor)
      }
    }

    # Check if identical. If not, dump results and error out
    identicalResults <- identicalScores(consistencyScores)
    if (length(identicalResults) > 1) {
      # Save full results, merge them, and save that too
      merged_results <- mergeResults(results)
      saveMergedResults(merged_results, outPath)

      # Make state space from the results
      results[[names(opt$runs)[i]]] <- run_res[1:j-1]
      stateSpaces <- makeStateSpaceFromResults(results)

      # Add the inconsistent attractors to it
      attractor_nodes <- c()
      for (k in 1:length(noLigRes)) {
        stateSpaces[[i]] <- igraph::union(stateSpaces[[i]], pathToGraph(cbind(noLigRes[[k]]$path, noLigRes[[k]]$attractor)))
        attractor_nodes <- union(attractor_nodes, apply(noLigRes[[k]]$attractor, 2, paste0, collapse=""))

        ligEnd <- paste0(run_res[[j - 1]]$lig$attractor[,k], collapse = "")
        attractor_nodes <- union(attractor_nodes, ligEnd)
        noLigStart <- paste0(noLigRes[[k]]$path[,1], collapse = "")
        stateSpaces[[i]] <- stateSpaces[[i]] + edge(c(ligEnd, noLigStart), action = "ligand_removal")
      }
      V(stateSpaces[[i]])[is.na(V(stateSpaces[[i]])$type)]$type <- "trajectory"
      V(stateSpaces[[i]])[attractor_nodes]$type <- "attractor"

      # Save state spaces
      for (i in 1:length(stateSpaces)) {
        stateSpaces[[i]] <- simplify(stateSpaces[[i]], remove.multiple = FALSE)
        write_graph(stateSpaces[[i]], file = paste0(outPath, "run_", names(stateSpaces)[i], "_state_space.graphml"), format = "graphml")
      }

      # Save initial states that led to inconsistent attractors, as well as full paths
      for (k in 1:length(identicalResults)) {
        idx <- identicalResults[[k]][1]
        run_initStates$state <- noLigRes[[idx]]$path[,1]
        write.csv(run_initStates, file = paste0(outPath, "inconsistent_initState_", k, ".csv"))
      }

      # Compare the first two of the inconsistent attractors
      comparePaths(noLigRes[[identicalResults[[1]][1]]]$attractor, noLigRes[[identicalResults[[2]][1]]]$attractor,
                   output = paste0(outPath, "inconsistent_attractor"))


      save(results, file = paste0(outPath, "results.RData"))
      stop("Inconsistent with-ligand attractors were found on round ", j, " of run ", names(opt$runs)[i],
           ". Dumping results as well as the initial states which resulted in inconsistent attractors.")
    }

    run_res[[j]]$lig <- ligRes[[1]]

    # Calculate score
    for (k in 1:j) {
      scoreLig[k,j] <- 1 - attractorDistance(run_res[[k]]$lig$attractor, run_res[[j]]$lig$attractor)
      scoreNoLig[k,j] <- 1 - attractorDistance(run_res[[k]]$noLig$attractor, run_res[[j]]$noLig$attractor)
    }

    if (j != 1 & any(scoreLig[1:(j-1), j] == 1) & any(scoreNoLig[1:(j-1), j] == 1)) {
      cat("Meta-attractor found for run", names(opt$runs)[i], "! Stopping simulation.", "\n")

      if (j > 2) {
        warning("Meta-attractor found for run ", names(opt$runs)[i], " after more than 2 rounds. This indicates possible traps/irreversible modifications.")
      }
      break
    } else if (j == opt$rounds) {
      cat("Max number of simulation rounds reached! Stopping simulation run", names(opt$runs)[i], "\n")
    }
  }

  results[[names(opt$runs)[i]]] <- run_res
  scores[[names(opt$runs)[i]]] <- list(noLig = scoreNoLig, lig = scoreLig)
}

################### Create state space representation ################
stateSpaces <- makeStateSpaceFromResults(results)

# Merge into a single global state space
globalStateSpace <- graph.empty()
for (g in stateSpaces) {
  globalStateSpace <- igraph::union(globalStateSpace, g)
}

# Copy over edge attributes into the global state space
E(globalStateSpace)$action <- NA
for (attr in edge_attr_names(globalStateSpace)[grepl("action_", edge_attr_names(globalStateSpace))]) {
  for (e in E(globalStateSpace)) {
    edge <- E(globalStateSpace)[[e]]
    new_attr <- get.edge.attribute(globalStateSpace, attr, edge)
    if (!is.na(new_attr)) {
      E(globalStateSpace)[[e]]$action <- new_attr
    }
  }
}

########### Create output RData file ##############
# Merge results into a single list, one entry per run
merged_results <- mergeResults(results)

# Create ordering for plots based on first simulation results
orderSimPath <- cbind(merged_results[[1]], 1:nrow(merged_results[[1]]))
for(i in (ncol(orderSimPath)-1):1){
  orderSimPath <- orderSimPath[order(orderSimPath[,i], decreasing = T), ]
}

# Move ligands to first position of order
ligands <- c()
for (i in 1:length(config$runs)) {
  ligands <- c(ligands, config$runs[[i]]$toggle)
}
ligands <- unique(ligands)
for (i in 1:length(ligands)) {
  # Find row numbers of unbound forms of ligand
  inds <- grep(ligands[i], rownames(orderSimPath))

  orderSimPath <- rbind(orderSimPath[inds,], orderSimPath[-inds,])
  # origInds <- inds
  #
  # # Move each of those rows to the beginning of the order (unless they're already at the front)
  # for (j in length(origInds):1) {
  #   if (origInds[j] > j) {
  #     if (origIinds[j] != 1) { # If not at top
  #       orderSimPath <- orderSimPath[c(inds[j], 1:(inds[j]-1), (inds[j]+1):nrow(orderSimPath)),]
  #     } else if (inds[j] == nrow(orderSimPath)) { # If all the way at the bottom
  #       orderSimPath <- orderSimPath[c(inds[j], 1:(inds[j]-1)),]
  #     }
  #     inds <- grep(ligands[i], rownames(orderSimPath))
  #   }
  # }
}

# Save order
plotOrder <- orderSimPath[,ncol(orderSimPath)]

# Apply order to all runs
for (i in 1:length(merged_results)) {
  merged_results[[i]] <- merged_results[[i]][plotOrder,]
}

for (i in 1:length(results)) {
  for (j in 1:length(results[[i]])) {
    results[[i]][[j]]$noLig$path <- results[[i]][[j]]$noLig$path[plotOrder,]
    results[[i]][[j]]$noLig$attractor <- results[[i]][[j]]$noLig$attractor[plotOrder,]
    results[[i]][[j]]$lig$path <- results[[i]][[j]]$lig$path[plotOrder,]
    results[[i]][[j]]$lig$attractor <- results[[i]][[j]]$lig$attractor[plotOrder,]
  }
}

################# Write things to file ###############
save(results, file = paste0(outPath, "results.RData"))

saveMergedResults(merged_results, outPath)

for (i in 1:length(stateSpaces)) {
  write_graph(stateSpaces[[i]], file = paste0(outPath, "run_", names(stateSpaces)[i], "_state_space.graphml"), format = "graphml")
}
write_graph(globalStateSpace, file = paste0(outPath, "global_state_space.graphml"), format = "graphml")

