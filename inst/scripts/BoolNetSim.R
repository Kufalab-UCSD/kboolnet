#!/usr/bin/env Rscript

options(stringsAsFactors = F)
suppressMessages(library(BoolNet))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(optparse))
library(kboolnet)

#############################################################
# Disclaimer
#
# This is a basic script that simulates a bipartite Boolean
# network created from a rxncon model. It requires three files:
#   model.boolnet
#   model_initial_vals.csv
#   model_symbols.csv
#
# The script runs two simulations and plots the trajectories.
#############################################################

################## Argument loading/parsing ################
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
default <- list(file=NA, driveFile=NA, on=c(), off=c(), initState=NA, out="./BoolNetSim_out")
opt <- setDefaults(config, default)

# Normalize paths
# Create out dir if it does not exist
if (!dir.exists(opt$out)) {
  dir.create(opt$out)
}
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
if (is.na(opt$initState)) {
  initStates <- read.csv(paste0(netFilePrefix, '_initial_vals.csv'), col.names=c("ID","state","name"), header=F)
} else {
  initStates <- read.csv(opt$initState, col.names=c("ID","state","name"), header=F)
}
initStates$name <- gsub("# ", "", initStates$name) # Clean up names
initStates$name <- gsub(" ", "", initStates$name)

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


################# Simulate and write to file ##################
# Get the path and the attractor
res <- getPathAndAttractor(network, initStates$state, initStates$name)

# Write to file
write.table(res$path, file=paste0(outPath, 'path.csv'), sep=",", col.names=F)
write.table(res$attractor, file=paste0(outPath, 'attractor.csv'), sep=",", col.names=F)
initStates$state <- res$attractor[,1]
write.table(initStates, file=paste0(outPath, 'new_initState.csv'), sep=",", col.names=F, row.names=F)

# Get orders for both the path and attractor
getOrder <- function(path) {
  orderSimPath <- cbind(path, 1:nrow(path))
  for(i in (ncol(orderSimPath)-1):1){
    orderSimPath <- orderSimPath[order(orderSimPath[,i], decreasing = T), ]
  }
  return(orderSimPath[,ncol(orderSimPath)])
}
pathOrder <- getOrder(res$path)
attractorOrder <- getOrder(res$attractor)

# Save plots
plotPath(res$path[pathOrder,,drop=F], paste0(outPath, 'path.pdf'))
plotPath(res$attractor[attractorOrder,,drop=F], paste0(outPath, 'attractor.pdf'))

cat("Done! Output written to directory", outPath, "\n")
