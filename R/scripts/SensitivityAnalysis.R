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
suppressMessages(library(egg))
suppressMessages(library(dplyr))
suppressMessages(library(openxlsx))
suppressMessages(library(googledrive))
suppressMessages(library(optparse))
suppressMessages(library(tidyr))

################ Function definitions #################
# Performs sensitivity analysis for all/selected components from a given initial state
KOanalysis <- function(network, names, symbols, states, components = c()) {
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

  # diffs <- diffs[, colnames(diffs) %in% outputs] # Keep only the requested output nodes
  # diffs$MSE <- rowMeans(diffs^2)

  return(diffs)
}

# Similar to KOanalysis, except doing it on a per-reaction basis
reactionAnalysis <- function(network, names, symbols, states, reactions = c()) {
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

  # diffs <- diffs[, colnames(diffs) %in% outputs] # Keep only the requested output nodes
  # diffs$MSE <- rowMeans(diffs^2)

  return(diffs)
}

plotDiffs <- function(diffs, file = NA) {
  # Processing before plotting
  diffs <- as.data.frame(diffs)
  diffs <- diffs[rowSums(diffs == 0) != ncol(diffs),] # Keep only rows which are different
  diffsGather <- diffs
  diffsGather$MSE <- rowMeans(diffs^2)
  diffsGather$KO <- rownames(diffsGather) # Add rownames as a column
  diffsGather <- gather(diffsGather, key = "signal", value = "diff", -KO, -MSE) # Gather into key-values
  diffsGather$KO <- factor(diffsGather$KO, levels=unique(diffsGather$KO[order(diffsGather$MSE, decreasing = TRUE)])) # Order KO by decreasing MSE

  # Plot the differences
  diffPlot <- ggplot(data = diffsGather, aes(x = 0.5, y = 0.5, height = 1 , width = 1)) +
    facet_grid(cols = vars(signal), rows = vars(KO), switch = "y") +
    geom_tile(aes(fill = diff), colour = "black", size = 1) +
    scale_fill_gradient2(limits = c(-1, 1), low = "red", mid = "white", high = "steelblue") +
    scale_y_continuous(limits = c(0,1), expand = c(0,0)) + scale_x_continuous(limits = c(0,1), expand = c(0,0)) +
    # xlab("Output") + ylab("KO") +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), strip.background = element_blank(),
          strip.text.x = element_text(angle = 90, hjust = 0), strip.text.y = element_text(angle = 180, hjust = 0))

  if (!is.na(file)) {
    width = length(unique(diffsGather$signal)) * .5 + 3
    height = length(unique(diffsGather$KO)) * .5 + 2
    ggsave(filename = file, plot = diffPlot, width = width, height = height)
  }

  return(diffPlot)
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
  opt <- loadConfig(opt, config)
}

# Set default args if they are not already set
default <- list(modules="", out="./out/", minQuality=0, ligands=NA, file=NA, driveFile=NA, outputs=NA, inhib="")
opt <- setDefaults(opt, default)

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
suppressMessages(source(paste0(kboolnetPath, "R/functions/fixedNetwork.R")))
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
  network <- fixedNetwork(network, inhib, symbolMapping$name, symbolMapping$ID)
}

# Create init states without the ligand
noLigStates <- initStates
for (lig in ligands) {
  if (sum(grepl(paste0("^", lig, "(_|$).*"), noLigStates$name)) == 0) {
    stop(paste0("Ligand ", lig, " not found!"))
  }
  noLigStates$state[grepl(paste0("^", lig, "(_|$).*"), noLigStates$name)] <- 0
}

################ Doing the analyses ####################
# KO analysis w/out ligand
cat("Performing component KO sensitivity analysis without ligand... ")
noLigKO <- KOanalysis(network, symbolMapping$name, symbolMapping$ID, noLigStates$state)
write.csv(noLigKO, file = paste0(outPath, "no_ligand_KO_analysis.csv"))
noLigKO <- noLigKO[,colnames(noLigKO) %in% outputs] # Keep only the outputs for plotting
if(sum(noLigKO) == 0) {
  cat("No difference in outputs detected between KO and non-KO attractors, skipping plotting.", "\n")
} else {
  plotDiffs(noLigKO, paste0(outPath, "no_ligand_KO_analysis.pdf"))
  cat("Done.", "\n")
}

# KO analysis w/ ligand
cat("Performing component KO sensitivity analysis with ligand... ")
ligKO <- KOanalysis(network, symbolMapping$name, symbolMapping$ID, initStates$state)
write.csv(ligKO, file = paste0(outPath, "ligand_KO_analysis.csv"))
ligKO <- ligKO[,colnames(ligKO) %in% outputs] # Keep only the outputs
if(sum(ligKO) == 0) {
  cat("No difference in ouputs detected between KO and non-KO attractors, skipping plotting.", "\n")
} else {
  plotDiffs(ligKO, paste0(outPath, "ligand_KO_analysis.pdf"))
  cat("Done.", "\n")
}

# Reaction analysis w/out ligand
cat("Performing reaction inhibition sensitivity analysis without ligand... ")
noLigReaction <- reactionAnalysis(network, symbolMapping$name, symbolMapping$ID, noLigStates$state)
write.csv(noLigReaction, file = paste0(outPath, "no_ligand_reaction_analysis.csv"))
noLigReaction <- noLigReaction[,colnames(noLigReaction) %in% outputs] # Keep only the outputs for plotting
if(sum(noLigReaction) == 0) {
  cat("No difference in ouputs detected between inhibited and uninhibted attractors, skipping plotting.", "\n")
} else {
  plotDiffs(noLigReaction, paste0(outPath, "no_ligand_reaction_analysis.pdf"))
  cat("Done.", "\n")
}

# Reaction analysis w/ ligand
cat("Performing reaction inhibition sensitivity analysis with ligand... ")
ligReaction <- reactionAnalysis(network, symbolMapping$name, symbolMapping$ID, initStates$state)
write.csv(noLigReaction, file = paste0(outPath, "ligand_reaction_analysis.csv"))
ligReaction <- ligReaction[,colnames(ligReaction) %in% outputs] # Keep only the outputs for plotting
if(sum(ligReaction) == 0) {
  cat("No difference in ouputs detected between inhibited and uninhibted attractors, skipping plotting.", "\n")
} else {
  plotDiffs(ligReaction, paste0(outPath, "ligand_reaction_analysis.pdf"))
  cat("Done.", "\n")
}
