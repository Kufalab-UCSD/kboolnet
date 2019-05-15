options(stringsAsFactors = F)
invisible(library(BoolNet))
invisible(source("./plotPath.R"))

##############################################################
# SimplifyAttractorEx.R 
# Adrian C
#
# Short example script that shows how to extract only certain
# nodes from a BoolNet attractor matrix for easier visualization
# with plotPath()
#
##############################################################

############### Network loading and simulation ####################
# Load network
netFilePrefix <- "~/Programming/internship/rxncon/CCR2_migration/v11/CCR2_migration_11"
network       <- loadNetwork(paste0(netFilePrefix, ".boolnet"), symbolic=TRUE)

# Load symbols and remove spaces
symbolMapping       <- read.csv(paste0(netFilePrefix, "_symbols.csv"), header=F, col.names=c("ID", "name"))
symbolMapping$name  <- gsub("[[:space:]]", "", symbolMapping$name)
symbolMapping       <- symbolMapping[,1:2]

# Load initial states from file
initStates      <- read.csv(paste0(netFilePrefix, '_initial_vals.csv'), col.names=c("ID","state","name"), header=F)
initStates$name <- gsub("# ", "", initStates$name) # Clean up names
initStates$name <- gsub(" ", "", initStates$name)

# Simulate to attractor
attractor           <- t(getPathToAttractor(network, initStates$state)) # NOTE that the output of getPath is transposed here
rownames(attractor) <- initStates$name

################ "Simplifying" attractor #####################

nodes   <- c("CCR2_[lig]_i+_CCL2_[rec]", "CCL2_[rec]--CCR2_[lig]", "PI3K_p+_PIP2_[(P)]", "PIP2_[(P)]-{p}",  # Nodes to be graphed
              "cdc42GEF_gef_cdc42_[(GnP)]", "cdc42_[(GnP)]-{gtp}", "PRex1_gef_Rac_[(GnP)]", "Rac_[(GnP)]-{gtp}", "[Migration]")
nodesInd <- numeric(length(nodes)) # Vector to store row numbers of selected nodes

# Get row numbers of selected nodes
for(i in 1:length(nodes)) { 
  nodesInd[i] <- grep(nodes[i], rownames(attractor), fixed = TRUE)[1]
}

simpleAttractor <- attractor[nodesInd,] # Get only selected rows of attractor
plotPath(simpleAttractor)
