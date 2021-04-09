#!/usr/bin/env Rscript

#################################################
# VerifyModel.R
# Adrian C
#
# Script to download model from Google Drive and
# simualte it using NFsim
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
library(kboolnet)

################ Function definitions #################
insert <- function(vector, elems, index) {
  if (index == 1) {
    return(c(elems, vector))
  } else if (index == length(vector) + 1) {
    return(c(vector, elems))
  } else {
    return(c(vector[1:(index - 1)], elems, vector[index:length(vector)]))
  }
}

################# Argument parsing #################
# Get commandline args
option_list = list(
  make_option(c("--config", "-c"), action="store", default=NA, type="character",
              help="Path of config file. You can specify parameters here instead of passing them as command-line
              arguments"),
  make_option("--BNGPath", action="store", default=NA, type="character",
              help="Path to BioNetGen install directory"),
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
              help="Comma-separated rxncon name(s) of ligand component(s) to be toggled in simulation"),
  make_option("--inhib", action="store", default=NA, type="character",
              help="Comma-separated rxncon names of reactions to be inhibited in simulation"),
  make_option("--KO", action="store", default=NA, type="character",
              help="Comma-separated rxncon names of components to be knocked out in simulation"),
  make_option("--time", action="store", default=NA, type="numeric",
              help="Time to run the simulation")
)
opt <- parse_args(OptionParser(option_list=option_list))
opt <- opt[!is.na(opt)] # Discard NA values

# Load config file if provided
if ("config" %in% names(opt)) {
  opt <- loadConfig(opt)
}

# Set default args if they are not already set
default <- list(modules="", out="./out/", minQuality=0, ligands=NA, file=NA, driveFile=NA, time=20, inhib="", KO="")
opt <- setDefaults(opt, default)

# Create out dir if it does not exist
if (!dir.exists(opt$out)) {
  dir.create(opt$out)
}

# Normalize paths
outPath       <- paste0(normalizePath(opt$out), "/")

# Parse modules option to a list
modules <- trimws(strsplit(opt$modules, ",")[[1]])

# Same for ligands option
if (is.na(opt$ligands)) {
  stop("Please provide ligand(s) to be toggled in simulation rounds")
}
ligands <- trimws(strsplit(opt$ligands, ",")[[1]])

# Same for inhibitors and KOs
inhib <- trimws(strsplit(opt$inhib, ",")[[1]])
KOs <- trimws(strsplit(opt$KO, ",")[[1]])

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

simTime <- opt$time
if (opt$time <= 0) {
  stop("time must be > 0")
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
callExtractModules(masterFile, modulesFile, modules, minQuality)

# Pass files to rxncon for processing
cat("Creating bngl files... ")
netFilePrefix <- gsub("\\.xlsx$", "", modulesFile)
callRxncon2BNG(modulesFile, netFilePrefix)

############### Generate NFsim files ######################
# Read in the bngl file
bngl <- readLines(paste0(netFilePrefix, ".bngl"))

# Add functions section and the new parameters to bngl
bngl <- insert(bngl, c("begin functions", "end functions"), which(bngl == "begin reaction rules"))
bngl <- insert(bngl, c("lig_prod 0", "lig_deg 0", "MaxLigand 1000"), which(bngl == "end parameters"))

# Generate everything we need to toggle the ligands on and off
seed_species_idxs <- which(grepl("seed species", bngl))
seed_species <- bngl[seed_species_idxs[1]:seed_species_idxs[2]]
for (lig in ligands) {
  # First we find the neutral representation of the ligand
  neutral_lig <- regmatches(seed_species, regexpr(paste0("^", lig, "\\(.*\\)\\s"), seed_species))

  # Check that the ligand was found
  if (length(neutral_lig) == 0) {
    stop(paste("Ligand ", lig, "not present in bngl file!"))
  }

  # Create observable
  bngl <- insert(bngl, paste("Molecules", lig, paste0(lig, "()"), sep="\t"), which(bngl == "end observables"))

  # Create functions to calculate rates
  prod_func <- paste0(lig, "_prod_func() = lig_prod * abs(MaxLigand - ", lig, ")")
  deg_func <- paste0(lig, "_deg_func() = lig_deg * ", lig)
  bngl <- insert(bngl, c(prod_func, deg_func), which(bngl == "end functions"))

  # Create synthesis/degradation reactions
  syn <- paste0("0 -> ", neutral_lig, " ", lig, "_prod_func()")
  deg <- paste0(lig, " -> 0 ", lig, "_deg_func() DeleteMolecules")
  bngl <- insert(bngl, c(syn, deg), which(bngl == "end reaction rules"))

  # Set initial concentration of the ligand to 0
  bngl[grepl(paste0("^Num", lig, "\\s"), bngl)] <- paste0("Num", lig, " 0")
}

# Adjust reaction rates
param_idxs <- which(bngl == "begin parameters" | bngl == "end parameters")
reaction_param_idxs <- which(grepl("^k([_][0-9]+)+", bngl))
reaction_param_idxs <- reaction_param_idxs[reaction_param_idxs > param_idxs[1] & reaction_param_idxs < param_idxs[2]]

rxnMapFile <- paste0(netFilePrefix, "_mapping.csv")
callReactionMapping(modulesFile, rxnMapFile)
mapping <- read.csv(rxnMapFile)

for (i in 1:nrow(mapping)) {
  rxn_idx <- intersect(which(grepl(mapping$rxnconName[i], bngl, fixed=TRUE)), reaction_param_idxs)
  if (length(rxn_idx) > 1) {
    stop(paste0("Reaction ", mapping$rxnconName[i], " is non-unique!"))
  } else if (length(rxn_idx) == 0) {
    stop(paste0("Reaction ", mapping$rxnconName[i], " is not present in bngl file."))
  }
  bngl[rxn_idx] <- gsub("\\s([0-9]+|[0-9]+\\.|\\.[0-9]+|[0-9]+\\.[0-9]+|\\.[0-9]+)\\s", mapping$rate[i], bngl[rxn_idx])
}

# Inhibit the reactions with inhib
for (inhibitor in inhib) {
  inhib_idxs <- intersect(which(grepl(inhibitor, bngl)), reaction_param_idxs)
  if (length(inhib_idxs) < 1) {
    stop(paste0("Reaction ", inhibitor, " is not present in the bngl file."))
  }
  bngl[inhib_idxs] <- gsub("\\s([0-9]+|[0-9]+\\.|\\.[0-9]+|[0-9]+\\.[0-9]+|\\.[0-9]+)\\s", "0.0", bngl[inhib_idxs])
}

# Remove KO components from the system
init_param_idxs <- which(grepl("^Num[0-9A-Za-z]+\\s+[0-9]+", bngl))
init_param_idxs <- init_param_idxs[init_param_idxs > param_idxs[1] & init_param_idxs < param_idxs[2]]
for (KO in KOs) {
  KO_idxs <- intersect(which(grepl(paste0("^Num", KO, "\\s"), bngl)), init_param_idxs)
  if (length(KO_idxs) < 1) {
    stop(paste0("Knocked-out component ", KO, " is not present in the bngl file."))
  }
  bngl[KO_idxs] <- gsub("\\s[0-9]+", " 0", bngl[KO_idxs])
}

# Change the last line to output XML instead
bngl[length(bngl)] <- "writeXML();"

# Add some spaces after the section ends
# ends <- which(grepl("^end ", bngl))
# for (i in 1:length(ends)) {
#   bngl <- insert(bngl, "", ends[i] + i)
# }

# Write the new bngl to file and generate the XML
writeLines(bngl, paste0(netFilePrefix, ".bngl"))
callBNG(paste0(netFilePrefix, ".bngl"), outPath)

# Generate the rnf file and run the simulation
rnf <- c(paste0("-xml ", netFilePrefix, ".xml"),
         paste0("-o ", netFilePrefix, ".gdat"),
         "-v",
         "begin",
         paste("sim", simTime/4, "100"),
         "set lig_prod 1000",
         "update",
         paste("sim", simTime/2, "100"),
         "set lig_deg 1000",
         "set lig_prod 0",
         "update",
         paste("sim", simTime/4, "100"),
         "end")
writeLines(rnf, paste0(netFilePrefix, ".rnf"))
cat("Done.", "\n")
cat("Simulating the model for ", simTime, " seconds...")
callNFSim(paste0(netFilePrefix, ".rnf"))

################# Read and plot results ###################
# We need to remove the # at the beginning of the gdat file to read it in properly
cat("Plotting results... ")
gdat <- readLines(paste0(netFilePrefix, ".gdat"))
gdat[1] <- gsub("#", "", gdat[1])
tmp <- tempfile()
writeLines(gdat, tmp)
data <- read.table(tmp, header = TRUE)

# Now we can plot
data_pivot <- pivot_longer(data, -time)
p <- ggplot(data_pivot, aes(x=time, y=value, color=name)) + geom_line() +
  geom_vline(xintercept = simTime/4) + geom_vline(xintercept = simTime * (3/4)) +
  labs(color="Observable")
ggsave(paste0(outPath, "results.pdf"), width = 7, height = 5, p)

# Write to file as well
write.csv(data, paste0(netFilePrefix, ".csv"))

cat("Done.", "\n")
