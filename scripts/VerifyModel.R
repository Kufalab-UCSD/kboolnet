#!/usr/bin/env Rscript

# This is the root directory of the project. Should be set by an install script of sorts in the future.
rootDir <- "<YOUR ROOT DIR HERE>"

options(stringsAsFactors = F)
suppressMessages(library(BoolNet))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(openxlsx))
suppressMessages(library(googledrive))
suppressMessages(library(optparse))
suppressMessages(library(tidyr))
suppressMessages(source(paste0(rootDir, "functions/extractModules.R")))
suppressMessages(source(paste0(rootDir, "functions/plotPath.R")))

#################################################
# VerifyModel.R
# Adrian C
#
# Script to download model from Google Drive
# and run a "sanity check" round of simulations.
#################################################

################# Config ###########################
# TODO: find a better way of implementing this
# Path to rxncon scripts
rxnconPath <- "~/.local/bin/"

################# Argument parsing #################
option_list = list(
  make_option("--file", action="store", default=NA, type="character",
              help="Path of master rxncon file (local)"),
  make_option("--driveFile", action="store", default=NA, type="character",
              help="File name or path of master rxncon file (on Google Drive)"),
  make_option("--modules", action="store", default="", type="character",
              help="Comma-separated modules to be loaded from master rxncon file [default: load all modules]"),
  make_option("--minQuality", action="store", default=0, type="integer",
              help="Minimum quality for rule to be loaded [default: %default]"),
  make_option("--out", action="store", default="./", type="character",
              help="Folder to which output files will be written [default: %default]"),
  make_option("--ligands", action="store", default=NA, type="character",
              help="Comma-separated rxncon name(s) of ligand node(s) to be toggled in verification simulation")
)
opt = parse_args(OptionParser(option_list=option_list))

# DEBUG ONLY
# opt = list(file=NA, out="asdf/", minQuality=0, file=NA, modules="",
           # ligands="CCL2_[rec]--0", driveFile="rxncon test")

# Ensure path ends with /
outPath <- suppressWarnings(normalizePath(opt$out))
outPath <- gsub("/$", "/", outPath)
outPath <- paste0(outPath, "/")

# Create out dir if it does not exist
if (!dir.exists(outPath)) {
  dir.create(outPath)
}

# Parse modules option to a list
modules <- strsplit(opt$modules, ",")[[1]]
modules <- gsub("^ *", "", modules) # Remove leading spaces
modules <- gsub(" *$", "", modules) # Remove trailing spaces

# Same for ligands option
ligands <- strsplit(opt$ligands, ",")[[1]]
ligands <- gsub("^ *", "", ligands) # Remove leading spaces
ligands <- gsub(" *$", "", ligands) # Remove trailing spaces

# Make sure input file path was provided
if (is.na(opt$file) && is.na(opt$driveFile)){ # If neither file was provided
  stop("Please provide a path to a local rxncon file with --file or a Google Drive file with --driveFile")
} else if ((!is.na(opt$file)) && (!is.na(opt$driveFile))) { # If both files were provided
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
  
  # Find the file in Drive
  gDriveID <- drive_find(pattern = opt$driveFile, type = "spreadsheet")[1,]
  if(is.na(gDriveID$name)) { # If file does not exist
    stop("rxncon file does not exist in Google Drive")
  }
  
  # Download file
  masterFile  <- paste0(outPath, gDriveID$name, ".xlsx")
  drive_download(gDriveID, path = masterFile, type = "xlsx", overwrite = T) # Download the file
  cat("Downloaded.", "\n")
  
# If local file provided
} else {
  masterFile <- suppressWarnings(normalizePath(opt$file))
  
  # File verification
  if (!grepl("\\.xlsx$", masterFile)) { # Make sure file is Excel file
    stop("rxncon file must be an Excel file (.xslx extension)")
  } else if (!file.exists(masterFile)) { # Make sure file exists
    stop("rxncon file does not exist")
  }
}

# Extract modules from master file, write to modules file
modulesFile <- paste0(outPath, "modules.xlsx")
extractModules(inPath = masterFile, outPath = modulesFile, minQuality = minQuality, modules = modules)
cat("Modules written to", modulesFile, "\n")

# Pass files to rxncon for processing
command <- paste0("cd ", outPath, " && ",
                  "python3 ", rxnconPath, "rxncon2boolnet.py ", modulesFile)
system(command)

################# Load BoolNet files ######################
# Load network
netFilePrefix <- gsub("\\.xlsx$", "", modulesFile)
network    <- loadNetwork(paste0(netFilePrefix, ".boolnet"), symbolic=TRUE)

# Load symbols and remove spaces
symbolMapping       <- read.csv(paste0(netFilePrefix, "_symbols.csv"), header=F, col.names=c("ID", "name"))
symbolMapping$name  <- gsub("[[:space:]]", "", symbolMapping$name)
symbolMapping       <- symbolMapping[,1:2]

# Load initial states from file
initStates <- read.csv(paste0(netFilePrefix, '_initial_vals.csv'), col.names=c("ID","state","name"), header=F)
initStates$name <- gsub("# ", "", initStates$name) # Clean up names
initStates$name <- gsub(" ", "", initStates$name)

####################### Simulate #########################
### This is hacky and temporary, need to rewrite

# First round
initStates$state[initStates$name %in% ligands] <- 0 # Remove ligands
path1 <- getPathToAttractor(network, initStates$state) %>% t()
initStates$state <- path1[,ncol(path1)] # Take last value of sim as new init state
path1Attract <- getPathToAttractor(network, initStates$state) %>% t()
rownames(path1) <- symbolMapping$name
rownames(path1Attract) <- symbolMapping$name
plotPath(path1)
plotPath(path1Attract, paste0(outPath, "1.pdf"))

# Second round
initStates$state[initStates$name %in% ligands] <- 1 # Add ligands
path2 <- getPathToAttractor(network, initStates$state) %>% t()
initStates$state <- path2[,ncol(path2)] # Take last value of sim as new init state
rownames(path2) <- symbolMapping$name
plotPath(path2, paste0(outPath, "2.pdf"))

# Third round
initStates$state[initStates$name %in% ligands] <- 0 # Remove ligands
path3 <- getPathToAttractor(network, path1[,ncol(path1)]) %>% t()
rownames(path3) <- symbolMapping$name
plotPath(path3, paste0(outPath, "3.pdf"))
