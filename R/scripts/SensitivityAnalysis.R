#!/usr/bin/env Rscript

#################################################
# SensitivityAnalysis.R
# Adrian C
#
# Script to run a "sensitivity analysis" on a
# rxncon network to see how modulating nodes/components
# affects the entire network.
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

################ Function definitions #################
# Performs sensitivity analysis for all/selected components from a given initial state
KOanalysis <- function(network, names, symbols, states, outputs, components = c()) {
  # If components was not defined, include all components in the network
  if (length(components) == 0) {
    allStates <- symbolMapping$name[startsWith(symbolMapping$ID, "S") &! startsWith(symbolMapping$name, "[")]
    components <- unique(gsub("_.*", "", allStates))
  }
  
  # Create df to store results 
  diffs <- as.data.frame(matrix(ncol = length(names), nrow = length(components)))
  rownames(diffs) <- components
  colnames(diffs) <- names

  # Simulate with the original network
  origAttr <- getPathAndAttractor(network, states, names)$attractor
  origAttrMeans <- rowMeans(origAttr)
  
  # Simulate with each component knocked out
  for (i in 1:length(components)) {
    KOnodes <- symbols[grepl(paste0(components[i], "(_|#|$)"), names)]
    KOnetwork <- fixGenes(network, KOnodes, 0)
    
    KOAttr <- getPathAndAttractor(KOnetwork, states, names)$attractor
    KOAttrMeans <- rowMeans(KOAttr)
    
    # Save the difference to diffs
    diffs[i,] <- KOAttrMeans - origAttrMeans
  }
  
  diffs <- diffs[, colnames(diffs) %in% outputs] # Keep only the requested output nodes
  diffs$MSE <- rowMeans(diffs^2)
  
  return(diffs)
}

# Similar to KOanalysis, except doing it on a per-reaction basis
reactionAnalysis <- function(network, names, symbols, states, outputs, reactions = c()) {
  # If components was not defined, include all reactions in the network
  if (length(reactions) == 0) {
    reactions <- symbolMapping$name[startsWith(symbolMapping$ID, "R")]
  }
  
  # Create df to store results 
  diffs <- as.data.frame(matrix(ncol = length(names), nrow = length(reactions)))
  rownames(diffs) <- reactions
  colnames(diffs) <- names

  # Simulate with the original network
  origAttr <- getPathAndAttractor(network, states, names)$attractor
  origAttrMeans <- rowMeans(origAttr)
  
  # Simulate with each reaction inhibited
  for (i in 1:length(reactions)) {
    reactionNode <- symbols[names == reactions[i]]
    inhibNetwork <- fixGenes(network, reactionNode, 0)
    
    inhibAttr <- getPathAndAttractor(inhibNetwork, states, names)$attractor
    inhibAttrMeans <- rowMeans(inhibAttr)
    
    # Save the difference to diffs
    diffs[i,] <- inhibAttrMeans - origAttrMeans
  }
  
  diffs <- diffs[, colnames(diffs) %in% outputs] # Keep only the requested output nodes
  diffs$MSE <- rowMeans(diffs^2)
  
  return(diffs)
}
################# Argument parsing #################
# Get commandline args
option_list = list(
  make_option(c("--config", "-c"), action="store", default=NA, type="character",
              help="Path of config file. You can specify parameters here instead of passing them as command-line
              arguments"),
  make_option("--kboolnetPath", action="store", default=NA, type="character",
              help="Path to root directory of kboolnet repository"),
  make_option("--rxnconPath", action="store", default=NA, type="character",
              help="Path to directory containing rxncon scripts"),
  make_option("--file", action="store", default=NA, type="character",
              help="Path of master rxncon file (local)"),
  make_option("--driveFile", action="store", default=NA, type="character",
              help="File name or path of master rxncon file (on Google Drive)"),
  make_option("--modules", action="store", default=NA, type="character",
              help="Comma-separated modules to be loaded from master rxncon file [default: load all modules]"),
  make_option("--minQuality", action="store", default=NA, type="integer",
              help="Minimum quality for rule to be loaded [default: 0]"),
  make_option(c("--out", "-o"), action="store", default=NA, type="character",
              help="Folder to which output files will be written [default: ./out/]"),
  make_option(c("--ligands", "-l"), action="store", default=NA, type="character",
              help="Comma-separated rxncon name(s) of ligand component(s) to be toggled in simulation."),
  make_option("--inhib", action="store", default=NA, type="character",
              help="Comma-separated rxncon name(s) of node(s) to be inhibited in simulation"),
  make_option("--outputs", action="store", default=NA, type="character",
              help="Comma-separated rxncon name(s) of nodes to be considered as outputs when performing sensitivity analysis.")
)
opt <- parse_args(OptionParser(option_list=option_list))
opt <- opt[!is.na(opt)] # Discard NA values

# Load config file if provided
if ("config" %in% names(opt)) {
  source(opt$config)

  if (!exists("config")) {
    stop("No config object found in config file")
  }

  config <- config[!is.na(config)] # Discard NA values

  # Keep only config values that were not passed as command line options
  config <- config[!(names(config) %in% names(opt))]

  # Merge command line args and config values
  opt <- c(opt, config)
}

# Set default args if they are not already set
default <- list(modules="", out="./out/", minQuality=0, ligands=NA, file=NA, driveFile=NA, outputs=NA, inhib="")
default <- default[!(names(default) %in% names(opt))]
opt     <- c(opt, default)

# Create out dir if it does not exist
if (!dir.exists(opt$out)) {
  dir.create(opt$out)
}

# Normalize paths
outPath       <- paste0(normalizePath(opt$out), "/")
kboolnetPath  <- paste0(normalizePath(opt$kboolnetPath), "/")
rxnconPath    <- paste0(normalizePath(opt$rxnconPath), "/")

# Load functions
suppressMessages(source(paste0(kboolnetPath, "R/functions/plotPath.R")))
suppressMessages(source(paste0(kboolnetPath, "R/functions/unbindLigand.R")))
suppressMessages(source(paste0(kboolnetPath, "R/functions/compMatrix.R")))
suppressMessages(source(paste0(kboolnetPath, "R/functions/driveDownload.R")))
suppressMessages(source(paste0(kboolnetPath, "R/functions/inhibitedNetwork.R")))
suppressMessages(source(paste0(kboolnetPath, "R/functions/getPathAndAttractor.R")))

# Parse modules option to a list
modules <- trimws(strsplit(opt$modules, ",")[[1]])

# Same for ligands option
if (is.na(opt$ligands)) {
  stop("Please provide ligand(s) to be toggled in simulation rounds")
}
ligands <- trimws(strsplit(opt$ligands, ",")[[1]])

# Same for outputs
if (is.na(opt$outputs)) {
  stop("Please provide node(s) to be used as output in sensitivity analysis")
}
outputs <- trimws(strsplit(opt$outputs, ",")[[1]])

# Same for inhibitors
inhib <- trimws(strsplit(opt$inhib, ",")[[1]])

# Escape inhibitor regexes
inhib <- gsub('([.|()\\^{}+$?]|\\[|\\])', '\\\\\\1', inhib)
inhib <- gsub('\\*', '\\.\\*\\?', inhib)

# Make sure input file path was provided
if (is.na(opt$file) & is.na(opt$driveFile)){ # If neither file was provided
  stop("Please provide a path to a local rxncon file with --file or a Google Drive file with --driveFile")
} else if ((!is.na(opt$file)) & (!is.na(opt$driveFile))) { # If both files were provided
  stop("Please provide only one path with EITHER --file or --driveFile")
}

minQuality <- opt$minQuality
if (opt$minQuality < 0) {
  stop("minQuality must be >= 0")
}

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
cat("Extracting modules...", "\n")
suppressWarnings(stderr <- system2(command = "python3", args = c(paste0(kboolnetPath, "Python/extract_modules.py"), "--file", masterFile,
                                                                 "--modules", paste0('"', paste0(modules, collapse=","), '"'), "--quality", minQuality,
                                                                 "--output", modulesFile), stderr = TRUE, stdout = ""))
if (any(grepl("Error", stderr, ignore.case = TRUE))) {
  cat(paste(stderr, "\n"))
  stop("Error during module extraction. Please run extract_modules.py on its own with the -v DEBUG flag.")
}


# Pass files to rxncon for processing
netFilePrefix <- gsub("\\.xlsx$", "", modulesFile)
suppressWarnings(stderr <- system2("python3", args = c(paste0(rxnconPath, "rxncon2boolnet.py"), modulesFile, "--output",
                                                       netFilePrefix), stderr = TRUE, stdout = ""))
if (any(grepl("Error", stderr, ignore.case = TRUE))) {
  cat(paste(stderr, "\n"))
  stop("Error during BoolNet file generation. Please run rxncon2boolnet.py on its own with the -v DEBUG flag.")
}

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

# Create a BoolNet network in which nodes have been inhibited
if (length(inhib) > 0) {
  network <- inhibitedNetwork(network, inhib, symbolMapping$name, symbolMapping$ID)
}

# Create init states without the ligand
noLigStates <- initStates
for (lig in ligands) {
  noLigStates$state[grepl(paste0("^", lig, "(_|$).*"), noLigStates$name)] <- 0
}

################ Doing the analyses ####################
# KO analysis
cat("Performing component KO sensitivity analysis... ")
noLigKO <- KOanalysis(network, symbolMapping$name, symbolMapping$symbol, noLigStates$state, outputs)
ligKO <- KOanalysis(network, symbolMapping$name, symbolMapping$symbol, initStates$state, outputs)
write.csv(noLigKO, file = paste0(outPath, "no_ligand_KO_analysis.csv"))
write.csv(ligKO, file = paste0(outPath, "ligand_KO_analysis.csv"))
cat("Done.", "\n")

# Reaction analysis
cat("Performing reaction inhibition sensitivity analysis... ")
noLigReaction <- reactionAnalysis(network,symbolMapping$name, symbolMapping$symbol, noLigStates$state, outputs)
ligReaction <- reactionAnalysis(network,symbolMapping$name, symbolMapping$symbol, initStates$state, outputs)
write.csv(noLigReaction, file = paste0(outPath, "no_ligand_reaction_analysis.csv"))
write.csv(noLigReaction, file = paste0(outPath, "ligand_reaction_analysis.csv"))
cat("Done.", "\n")