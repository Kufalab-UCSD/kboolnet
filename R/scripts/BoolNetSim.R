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

# accept network as commandline argument
option_list = list(
  make_option(c("--model", "-m"), action="store", default=NA, type='character',
              help="Prefix of model to be simulated"),
  make_option(c("--noplot", "-n"), action="store_true", default=FALSE,
              help="Do not write plots to pdf. May be necessary for larger models")
)
opt = parse_args(OptionParser(option_list=option_list))

# Provide network name
if(!is.na(opt$model)) {
  filePrefix <- opt$model
} else {
  stop("Please provide the name of the model.")
}

# Load network
network    <- loadNetwork(paste0(filePrefix, ".boolnet"), symbolic=TRUE)

# Load symbols and remove spaces
symbolMapping       <- read.csv(paste0(filePrefix, "_symbols.csv"), header=F, col.names=c("ID", "name"))
symbolMapping$name  <- gsub("[[:space:]]", "", symbolMapping$name)
symbolMapping       <- symbolMapping[,1:2]

# Load initial states from file
initStates <- read.csv(paste0(filePrefix, '_initial_vals.csv'), col.names=c("ID","state","name"), header=F)
initStates$name <- gsub("# ", "", initStates$name) # Clean up names
initStates$name <- gsub(" ", "", initStates$name)

# Get the path and the attractor
res <- getPathAndAttractor(network, initStates$state, initStates$name)

# Write to file
write.table(res$path, file=paste0(filePrefix, '_path.csv'), sep=",", col.names=F)
write.table(res$attractor, file=paste0(filePrefix, '_attractor.csv'), sep=",", col.names=F)

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
if (!opt$noplot) {
  plotPath(res$path[pathOrder,,drop=F], paste0(filePrefix, '_path.pdf'))
  plotPath(res$attractor[attractorOrder,,drop=F], paste0(filePrefix, '_attractor.pdf'))
}
