#!/usr/bin/env Rscript

#################################################
# PlotModules.R
# Adrian C
#
# Script to download a model and plot each of its
# modules along with
#################################################

################# Library loading ##################
options(stringsAsFactors = F)
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(openxlsx))
suppressMessages(library(googledrive))
suppressMessages(library(optparse))
suppressMessages(library(tidyr))
suppressMessages(library(numbers))
suppressMessages(library(xml2))
library(kboolnet)

################ Function definitions #################


################# Argument parsing #################
# Get commandline args
option_list = list(
  make_option(c("--config", "-c"), action="store", default=NA, type="character",
              help="Path of config file. You can specify parameters here instead of passing them as command-line
              arguments"),
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
  make_option(c("--nodomains", "-n"), action="store_true", default=NA,
              help="Remove domains from graph outputs [default: don't remove]")
)
opt <- parse_args(OptionParser(option_list=option_list))
opt <- opt[!is.na(opt)] # Discard NA values

# Load config file if provided
if ("config" %in% names(opt)) {
  opt <- loadConfig(opt, config)
}

# Set default args if they are not already set
default <- list(modules="", out="./out/", minQuality=0, file=NA, driveFile=NA, nodomains=FALSE)
opt <- setDefaults(opt, default)

# Create out dir if it does not exist
if (!dir.exists(opt$out)) {
  dir.create(opt$out)
}

# Normalize paths
outPath       <- paste0(normalizePath(opt$out), "/")

# Parse modules option to a list
modules <- trimws(strsplit(opt$modules, ",")[[1]])

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

# If modules were not provided, get list of all modules from excel file
master_rxn <- read.xlsx(masterFile, sheet = 1)
modulesCol_rxn <- which(master_rxn[1,] == "!Module")
modules <- c(modules, unlist(strsplit(master_rxn[2:nrow(master_rxn), modulesCol_rxn], ",")))
master_cont <- read.xlsx(masterFile, sheet = 2)
modulesCol_cont <- which(master_cont[1,] == "!Module")
modules <- c(modules, unlist(strsplit(master_cont[2:nrow(master_cont), modulesCol_cont], ",")))
modules <- trimws(modules)
modules <- modules[!is.na(modules)]
modules <- unique(modules)

############## Module extraction and plotting ####################
for (module in modules) {
  # Extract the module
  modulesFile <- paste0(outPath, "module_", module, ".xlsx")
  cat(paste0("Extracting module ", module, "..."), "\n")
  callExtractModules(masterFile, modulesFile, module)

  cat(paste0("Plotting module ", module, "..."), "\n")
  regFile <- paste0(outPath, "module_", module, "_reg.xgmml")
  callRxncon2Reg(modulesFile, regFile)


  if (opt$nodomains) {
    removeDomainsXGMML(regFile, regFile)
  }
}

############# Generate node - module mapping #####################
# Get all module file
modulesFile <- paste0(outPath, "all_modules.xlsx")
cat("Extracting all modules ...", "\n")
callExtractModules(masterFile, modulesFile, "")

# Plot regulatory graph
regFile <- paste0(outPath, "all_modules_reg.xgmml")
cat(paste0("Plotting all modules..."), "\n")
callRxncon2Reg(modulesFile, regFile)

if (opt$`nodomains`) {
  removeDomainsXGMML(paste0(outPath, "all_modules_reg.xgmml"),)
}

# Load xgmml file for all module file
all_modules <- paste0(outPath, "all_modules_reg.xgmml") %>% read_xml()
all_modules_nodes <- all_modules %>% xml_find_all(".//*[@name='rxnconID']/@value") %>% xml_text()
all_nodes <- vector(mode = "list", length = length(all_modules_nodes))
names(all_nodes) <- all_modules_nodes

# Load xgmml file for each module and append module name to list if node found in it
for (module in modules) {
  module_file <- paste0(outPath, "module_", module, "_reg.xgmml") %>% read_xml()
  module_nodes <- module_file %>% xml_find_all(".//*[@name='rxnconID']/@value") %>% xml_text()
  for (node in module_nodes) {
    all_nodes[[node]] <- c(all_nodes[[node]], module)
  }
}

# Make a node to module mapping
mapping <- data.frame(node = names(all_nodes), module = sapply(all_nodes, paste0, collapse = ","))
write.csv(mapping, file = paste0(outPath, "node_modules.csv"), row.names = F)

cat("Done!", "\n")