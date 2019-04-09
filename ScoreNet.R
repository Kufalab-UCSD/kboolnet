options(stringsAsFactors = F)
invisible(library(dplyr))
invisible(library(BoolNet))
invisible(library(tidyr))
invisible(library(CellNOptR))
invisible(library(readr))
invisible(library(plyr))
invisible(source("./plotPath.R"))

#######################################################
# ScoreNet.R
# Adrian C
#
# Script to score similarity between BoolNet 
# simulation of a network and experimental data about 
# the network.
#
#######################################################

########## Data and network loading #################

# accept network as commandline argument
# args <- commandArgs(trailingOnly = TRUE)
# if(length(args) == 0) {
#   stop("Please provide the name of the model.")
# }
# 
# netFilePrefix <- args[1]

netFilePrefix <- "~/Programming/internship/rxncon/CCR2_migration/v11/CCR2_migration_11"
dataFilePath  <- "/home/adrian/Programming/internship/rxncon/CCR2_migration/v11/Experiment_CCR2_130319 - GrantExample_CXCR4.csv"

# Load network
network <- loadNetwork(paste0(netFilePrefix, ".boolnet"), symbolic=TRUE)

# Load symbols and remove spaces
symbolMapping       <- read.csv(paste0(netFilePrefix, "_symbols.csv"), header=F, col.names=c("ID", "name"))
symbolMapping$name  <- gsub("[[:space:]]", "", symbolMapping$name)
symbolMapping       <- symbolMapping[,1:2]

# Load initial states from file
initStates <- read.csv(paste0(netFilePrefix, '_initial_vals.csv'), col.names=c("ID","state","name"), header=F)
initStates$name <- gsub("# ", "", initStates$name) # Clean up names
initStates$name <- gsub(" ", "", initStates$name)

# # Load data file
# dataFile <- read.csv(dataFilePath, header=F)
# dataFile <- dataFile[,!sapply(dataFile, function(x) all(is.na(x) | x == ""))] # Remove empty columns
# dataFile <- dataFile[!apply(dataFile == "", 1, all),] # Remove empty rows
# 
# # Extract subtables
# headerRows <- grep("TR:", dataFile[,1]) # Get rows which are headers of subtables
# subtables <- vector("list", length(headerRows))
# if (length(subtables) > 1) {
#   for (i in 1:(length(subtables)-1)) {
#     subtables[i] <- list(dataFile[headerRows[i]:(headerRows[i+1]-1),]) # Add subtables to list
#   }
# }
# subtables[length(subtables)] <- list(dataFile[headerRows[length(subtables)]:nrow(dataFile),]) # Add last subtable to list
# 
# # Process column names for subtables
# for (i in 1:length(subtables)) {
#   # COME BACK TO THIS LATER
#   # # Column name adjustments (add quotes around rxncon names to prevent confusion w/ inhibitors)
#   # subtables[[i]][1,] <- sapply(subtables[[i]][1,], function(x) {
#   #   if(!grepl("CellLine", x)) {
#   #     gsub("TR:", 'TR:"', x)
#   #     
#   #     gsub("")
#   #   }
#   # })
#   # subtables[[i]][[1,!grepl("CellLine", subtables[[i]])]] <- gsub("TR:", 'TR:"', subtables[[i]][1,!grepl("CellLine", subtables[[i]])])
#   
#   # Set first row as column names
#   colnames(subtables[[i]]) <- subtables[[i]][1,]
#   subtables[[i]] <- subtables[[i]][-1,]
#   subtables[[i]] <- subtables[[i]][,!colnames(subtables[[i]]) == ""] # Remove empty rows
# }
# 
# # Join all tables into single table
# joinTable <- Reduce(function(dtf1, dtf2) full_join(dtf1, dtf2, all = TRUE), subtables)
# 
# joinTablePath <- gsub(".csv", "_joined.csv", dataFilePath)
# write.csv(joinTable, joinTablePath)
# 
# expDataMIDAS <- readMIDAS("/home/adrian/Programming/internship/rxncon/CCR2_migration/sheet1_joined.csv", verbose=TRUE) 
# expData <- makeCNOlist(expDataMIDAS, subfield=FALSE, verbose=TRUE)

expData <- readMIDAS(dataFilePath, verbose=TRUE) %>% makeCNOlist(subfield=TRUE, verbose=TRUE)
plotCNOlist(expData)

#############  First simulation: Find stable unperturbed state  #############

# Turn off all stimulator nodes
initStates$state[initStates$name %in% expData$namesStimuli] <- 0 

# Simulate to attractor and take last state as new steady state#
path <- t(getPathToAttractor(network, initStates$state))

# TODO: I'm pretty sure this is still broken
# # Sort in order of activation, guarantees consistent ordering of graphs
# pathOrder <- c(1:nrow(path))
# for(i in ncol(path):1){
#   newOrder    <- order(path[,i], decreasing = T)
#   path        <- path[newOrder, ]
#   pathOrder   <- pathOrder[newOrder]
# }
# initStates    <- initStates[pathOrder,]
# symbolMapping <- symbolMapping[pathOrder,]

# Save last state of path as new steady initial state
steadyInitStates       <- initStates
steadyInitStates$state <- path[,ncol(path)]

########### Simulate different perturbed states, scoring happens here too ###################

# Create dfs to store results
scores                <- data.frame(inhib=numeric(), stim=numeric())
attractors            <- data.frame(matrix(ncol = 6, nrow = ncol(expData$valueCues)))
colnames(attractors)  <- c("inhib", "stim", "inhibPath", "inhibAttract", "stimPath", "stimAttract")
simData               <- expData # Copy expData CNOList to store simulation data

for(i in 1:ncol(expData$valueCues)) { # For each combination of cues
  # Init empty rows in results dfs
  scores[i,] <- c(NA, NA)
  
  # Prepare inhibited network
  inhibOn         <- expData$namesInhibitors[which(expData$valueInhibitors[i,] == 1)] # See which inhibitors are present
  attractors$inhib[i] <- paste(inhibOn,",")
  inhibNodes      <- symbolMapping$ID[symbolMapping$name %in% inhibOn] # Get their BoolNet node IDs
  inhibNetwork    <- fixGenes(network, inhibNodes, 0) # Create new network with inhib nodes permanently off
  
  # Simulate to stable inhibited state
  inhibInitStates       <- steadyInitStates
  inhibPath             <- t(getPathToAttractor(inhibNetwork, steadyInitStates$state))
  # inhibPath             <- inhibPath[pathOrder,]
  rownames(inhibPath)   <- symbolMapping$name
  inhibInitStates$state <- inhibPath[,ncol(inhibPath)] # Save last attractor state as new init state
  
  attractors$inhibPath[i] <- list(inhibPath)
  
  # Simulate stable inhibited attractor
  inhibAttractor <- t(getPathToAttractor(inhibNetwork, inhibInitStates$state))
  inhibAttractor <- as.data.frame(inhibAttractor)
  # inhibAttractor <- inhibAttractor[pathOrder,]
  rownames(inhibAttractor)    <- symbolMapping$name
  
  attractors$inhibAttract[i]  <- list(inhibAttractor)
  
  # Compare inhibited attractor to experimental data
  inhibAttractor      <- cbind(inhibAttractor, rowMeans(inhibAttractor)) # Average attractor values
  colnames(inhibAttractor)[ncol(inhibAttractor)] <- "avg"
  inhibAttractor$name <- symbolMapping$name
  inhibData           <- data.frame(value=expData$valueSignals[[1]][i,]) # Get data from first timepoint
  inhibData$name      <- expData$namesSignals
  inhibData_Atrract   <- join(inhibData, inhibAttractor, by="name") # Match experimental data
  simData$valueSignals[[1]][i,] <- inhibData_Atrract$avg # Save attractor average to CNOList
  scores$inhib[i]     <- mean((inhibData_Atrract$value - inhibData_Atrract$avg)^2, na.rm=TRUE) # Calculate MSE
  
  # Prepare stimulated network
  stimOn          <- expData$namesStimuli[which(expData$valueStimuli[i] == 1)] # See which stimulators are present
  attractors$stim[i] <- paste(stimOn,",")
  stimInitStates  <- inhibInitStates
  stimInitStates$state[stimInitStates$name %in% expData$namesStimuli] <- 1 # Turn on stimulus nodes
  stimPath        <- t(getPathToAttractor(inhibNetwork, stimInitStates$state))
  # stimPath <- stimPath[pathOrder,]
  
  attractors$stimPath[i] <- list(stimPath)
  
  # Simulate stable stimulated attractor
  stimAttractor <- t(getPathToAttractor(inhibNetwork, stimPath[,ncol(stimPath)]))
  stimAttractor <- as.data.frame(stimAttractor)
  # stimAttractor <- stimAttractor[pathOrder,]
  rownames(stimAttractor) <- symbolMapping$name
  
  attractors$stimAttract[i] <- list(stimAttractor)
  
  # Compare stim attractor to experimental data
  stimAttractor       <- cbind(stimAttractor, rowMeans(stimAttractor)) # Average attractor values
  colnames(stimAttractor)[ncol(stimAttractor)] <- "avg"
  stimAttractor$name  <- symbolMapping$name
  stimData            <- data.frame(value=expData$valueSignals[[2]][i,]) # Get data from second timepoint
  stimData$name       <- expData$namesSignals
  stimData_Atrract    <- join(stimData, stimAttractor, by="name") # Match experimental data
  simData$valueSignals[[2]][i,] <- stimData_Atrract$avg
  scores$stim[i]      <- mean((stimData_Atrract$value - stimData_Atrract$avg)^2, na.rm=TRUE) # Calculate MSE
}

# Generate some graphs
toyNodes   <- c("CCR2_[lig]_i+_CCL2_[rec]", "CCL2_[rec]--CCR2_[lig]", "PI3K_p+_PIP2_[(P)]", "PIP2_[(P)]-{p}",
              "cdc42GEF_gef_cdc42_[(GnP)]", "cdc42_[(GnP)]-{gtp}", "PRex1_gef_Rac_[(GnP)]", "Rac_[(GnP)]-{gtp}", "[Migration]")
for(i in 1:length(toyNodes)) { # Get index numbers of toy nodes from symbol mapping
  toyNodes[i] <- grep(toyNodes[i], symbolMapping$name, fixed = TRUE)[1]
}
toyNodes <- as.numeric(toyNodes)
plotPath(attractors$stimAttract[[1]][toyNodes,])

