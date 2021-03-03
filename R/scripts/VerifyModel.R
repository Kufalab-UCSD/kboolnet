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
library(kboolnet)

################ Function definitions #################

plotCompMat <- function(mat) {
  # Plot comparison matrices (w/ ligand)
  df <- as.data.frame(mat)
  colnames(df) <- 1:ncol(df)
  rows <- nrow(df)
  df <- cbind(1:rows, df)
  colnames(df)[1] <- "y"
  df <- gather(df, "x", "value", 2:(rows+1), na.rm = T) # Rearrange data for plotting
  df <- transform(df, x = as.numeric(x))
  p <- ggplot(df, aes(x, y)) + geom_tile(aes(fill = value), color = "lightgrey") +
    scale_x_continuous(position = "top", expand = c(0,0), breaks=1:rows) + scale_y_reverse(expand=c(0,0), breaks=rows:1) +
    scale_fill_gradient2(limits=c(0,1), low = "red", mid="white", high="steelblue", midpoint = 0.5) +
    theme_classic() + theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
                            axis.text = element_text(size = 10), panel.grid.minor = element_line(color = "lightgrey")) + coord_equal()
  return(p)
}

simulateWithLigands <- function(network, states, names, ligands) {
  # Add ligands in fully neutral state
  for (lig in 1:length(ligands)){
    states[grepl(paste0(ligands[lig], "_.*--0$"), names)] <- 1
    states[grepl(paste0(ligands[lig], "_.*-\\{0\\}$"), names)] <- 1
  }

  # Simulate and return
  return(getPathAndAttractor(network, states, names))
}

simulateWithoutLigands <- function(network, states, names, ligands) {
  # Unbind ligands from complexes and remove any ligand
  for (lig in 1:length(ligands)){
    states <- unbindLigand(ligands[lig], names, states)
    states[names %in% ligNodesNames] <- 0
  }

  # Simulate and return
  return(getPathAndAttractor(network, states, names))
}

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
  make_option(c("--ligands", "-l"), action="store", default=NA, type="character",
              help="Comma-separated rxncon name(s) of ligand component(s) to be toggled in verification simulation"),
  make_option("--rounds", action="store", default=NA, type="integer",
              help="Maximum number of ligand/no-ligand simulation pairs that should be run (default: 20)"),
  make_option("--inhib", action="store", default=NA, type="character",
              help="Comma-separated rxncon name(s) of node(s)/component(s) to be inhibited in verification simulation")
)
opt <- parse_args(OptionParser(option_list=option_list))
opt <- opt[!is.na(opt)] # Discard NA values

# Load config file if provided
if ("config" %in% names(opt)) {
  opt <- loadConfig(opt, config)
}

# Set default args if they are not already set
default <- list(modules="", out="./out/", minQuality=0, ligands=NA, file=NA, driveFile=NA, rounds=20, inhib="")
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

maxrounds <- opt$rounds
if (opt$rounds < 2) {
  stop("rounds must be >= 2")
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
cat("Extracting modules... ")
path <- paste0(system.file(package="kboolnet"), "/python/extract_modules.py")
suppressWarnings(stderr <- system2(command = "python3", args = c(path, "--file", masterFile,
                                                                 "--modules", paste0('"', paste0(modules, collapse=","), '"'), "--quality", minQuality,
                                                                 "--output", modulesFile), stderr = TRUE, stdout = ""))
if (any(grepl("Error", stderr, ignore.case = TRUE))) {
  cat(paste(stderr, "\n"))
  stop("Error during module extraction. Please run extract_modules.py on its own with the -v DEBUG flag.")
}
cat("Done.\n")


# Pass files to rxncon for processing
cat("Generating BoolNet files... ")
netFilePrefix <- gsub("\\.xlsx$", "", modulesFile)
path <- paste0(system.file(package="kboolnet"), "/python/rxncon2boolnet.py")
suppressWarnings(stderr <- system2("python3", args = c(path, modulesFile, "--output",
                                                       netFilePrefix), stderr = TRUE, stdout = ""))
print(stderr)
if (any(grepl("Error", stderr, ignore.case = TRUE))) {
  cat(paste(stderr, "\n"))
  stop("Error during BoolNet file generation. Please run rxncon2boolnet.py on its own with the -v DEBUG flag.")
}
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

# Create a BoolNet network in which nodes have been inhibited
if (length(inhib) > 0) {
  network <- fixedNetwork(network, inhib, symbolMapping$name, symbolMapping$ID)
}

# Find the nodes that the ligands correspond to
ligNodes <- character()
ligNodesNames <- character()
for (i in 1:length(ligands)) {
  # Make sure the ligand exists in a neutral state
  if (!(any(grepl(paste0("^", ligands[i], "(_.*--0|_.*-\\{0\\})$"), initStates$name)))) {
    stop("No neutral state found for ligand ", ligands[i], ". Please verify that ", ligands[i], " is a valid component in the rxncon system.")
  }

  # Make a regex matching unbound, and modified forms of ligand
  ligRegex <- paste0("(", c(paste0(ligands[i], "_.*--0"), paste0("^", ligands[i], "$"),
                            paste0(ligands[i], "_.*-\\{.*\\}")), ")", collapse="|")
  ligMatch <- grepl(ligRegex, symbolMapping$name)
  ligNodes <- c(ligNodes, symbolMapping$ID[ligMatch])
  ligNodesNames <- c(ligNodesNames, symbolMapping$name[ligMatch])

  cat("Ligand", ligands[i], "matched to node(s)", paste0(symbolMapping$name[ligMatch], collapse=", "), "\n")
}

####################### Simulate #########################
# Lists to store simulation results
noLigAttr <- list()
noLigPath <- list()
ligAttr   <- list()
ligPath   <- list()

# Matrices to store comparison results
scoreLig   <- matrix(nrow = maxrounds, ncol = maxrounds)
scoreNoLig <- matrix(nrow = maxrounds, ncol = maxrounds)

# Simulation rounds
cat("Starting simulations...", "\n")
rounds <- 1

# First round to set things up
cat("Simulation round 1 ... ")
noLigResults <- simulateWithoutLigands(network, initStates$state, initStates$name, ligands)
noLigAttr[[1]] <- noLigResults$attractor
noLigPath[[1]] <- noLigResults$path
ligResults <- simulateWithLigands(network, noLigAttr[[1]][,ncol(noLigAttr[[1]])], initStates$name, ligands)
ligAttr[[1]] <- ligResults$attractor
ligPath[[1]] <- ligResults$path
scoreLig[1,1] <- 1
scoreNoLig[1,1] <- 1
cat("Done.", "\n")

for (i in 2:maxrounds) {
  rounds <- rounds + 1
  cat("Simulation round", rounds, "... ")

  roundLigAttr   <- list()
  roundLigPath    <- list()
  roundNoLigAttr <- list()
  roundNoLigPath  <- list()

  # Simulate without ligand using each state of previous w/ ligand attractor as an initial state
  for (j in 1:ncol(ligAttr[[i-1]])) {
    res <- simulateWithoutLigands(network, ligAttr[[i-1]][,j], initStates$name, ligands)
    roundNoLigAttr[[j]] <- res$attractor
    roundNoLigPath[[j]] <- res$path
  }

  # Compare no ligand attractors to each other
  roundNoLigCompScores <- matrix(nrow = length(roundNoLigAttr), ncol = length(roundNoLigAttr))
  for (j in 1:length(roundNoLigAttr)) {
    for (k in j:length(roundNoLigAttr)) {
      roundNoLigCompScores[j,k] <- 1 - attractorDistance(roundNoLigAttr[[j]], roundNoLigAttr[[k]])
    }
  }

  if (any(roundNoLigCompScores != 1, na.rm=T)) { # If dissimilar attractors found
    warning("No-ligand attractors in round ", rounds, " are not consistent, possibly indicating \"trap\" states! Using most dissimilar attractor to continue.")
    # Compare no ligand attractors for this round to previous round's no ligand attractor, save the most dissimilar attractor
    roundNoLigScores <- numeric(length(roundNoLigAttr))
    for (j in 1:length(roundNoLigAttr)) {
      roundNoLigScores[j] <- 1 - attractorDistance(roundNoLigAttr[[j]], noLigAttr[[i-1]])
    }
    noLigAttr[[i]] <- roundNoLigAttr[[which.min(roundNoLigScores)]]
    noLigPath[[i]] <- roundNoLigPath[[which.min(roundNoLigScores)]]
  } else { # Else just save the first one
    noLigAttr[[i]] <- roundNoLigAttr[[1]]
    noLigPath[[i]] <- roundNoLigPath[[1]]
  }

  # Simulate with ligand using each state of previous no ligand attractor as an initial state
  for (j in 1:ncol(noLigAttr[[i]])) {
    res <- simulateWithLigands(network, noLigAttr[[i]][,j], initStates$name, ligands)
    roundLigAttr[[j]] <- res$attractor
    roundLigPath[[j]] <- res$path
  }

  # Compare with ligand attractors to each other
  roundLigCompScores <- matrix(nrow = length(roundLigAttr), ncol = length(roundLigAttr))
  for (j in 1:length(roundLigAttr)) {
    for (k in j:length(roundLigAttr)) {
      roundLigCompScores[j,k] <- 1 - attractorDistance(roundLigAttr[[j]], roundLigAttr[[k]])
    }
  }

  if (any(roundLigCompScores != 1, na.rm=T)) { # If dissimilar attractors found
    warning("With-ligand attractors in round ", rounds, " are not consistent, possibly indicating \"trap\" states! Using most dissimilar attractor to continue.")
    # Compare no ligand attractors for this round to previous round's no ligand attractor, save the most dissimilar attractor
    roundLigScores <- numeric(length(roundLigAttr))
    for (j in 1:length(roundLigAttr)) {
      roundLigScores[j] <- 1 - attractorDistance(roundLigAttr[[i]], ligAttr[[i-1]])
    }
    ligAttr[[i]] <-roundLigAttr[[which.min(roundLigScores)]]
    ligPath[[i]] <-roundLigPath[[which.min(roundLigScores)]]
  } else { # Else just save the first one
    ligAttr[[i]] <- roundLigAttr[[1]]
    ligPath[[i]] <- roundLigPath[[1]]
  }

  cat("Done.", "\n")

  # Compare this round of simulation to previous rounds
  for (j in 1:rounds) {
    scoreLig[j,rounds] <- 1 - attractorDistance(ligAttr[[rounds]], ligAttr[[j]])
    scoreNoLig[j,rounds] <- 1 - attractorDistance(noLigAttr[[rounds]], noLigAttr[[j]])
  }

  # If both ligand and no-ligand simulations are identical to a previous round, then we know a
  # "meta-attractor" has been found and we can stop simulation
  if (rounds != 1 & any(scoreLig[1:rounds-1,rounds] == 1) & any(scoreNoLig[1:rounds-1,rounds] == 1)) {
    cat("Meta-attractor found! Stopping simulation.", "\n")

    if (rounds > 2) {
      warning("Meta-attractor found after more than 2 rounds. This indicates possible traps/irreversible modifications.")
    }
    break
  } else if (rounds == maxrounds) {
    cat("Max number of simulation rounds reached! Stopping simulation.", "\n")
  }
}

# Remove extra rows/columnds from comparison matrices
scoreLig   <- scoreLig[1:rounds,1:rounds]
scoreNoLig <- scoreNoLig[1:rounds,1:rounds]

##################### Plot and save data ##############################
cat("Simulations complete after", rounds, "rounds.", "\n")
# Create output subdirs
unlink(paste0(outPath, "lig"), recursive = TRUE)
unlink(paste0(outPath, "nolig"), recursive = TRUE)
dir.create(paste0(outPath, "lig/"))
dir.create(paste0(outPath, "lig/attractor/"))
dir.create(paste0(outPath, "lig/path/"))
dir.create(paste0(outPath, "nolig/"))
dir.create(paste0(outPath, "nolig/attractor/"))
dir.create(paste0(outPath, "nolig/path/"))

# Write comparison data
write.csv(scoreLig, file = paste0(outPath, "lig/simulation_comparison.csv"))
write.csv(scoreNoLig, file = paste0(outPath, "nolig/simulation_comparison.csv"))

# Plot the attractor comparison matrices
scoreLigPlot   <- plotCompMat(scoreLig)
scoreNoLigPlot <- plotCompMat(scoreNoLig)
width <- 1 + 0.75 * ncol(scoreLig)
height <- 1 + 0.75 * nrow(scoreLig)
suppressMessages(ggsave(paste0(outPath, "lig/simulation_comparison.pdf"), device="pdf", plot=scoreLigPlot, width=width, height=height))
suppressMessages(ggsave(paste0(outPath, "nolig/simulation_comparison.pdf"), device="pdf", plot=scoreNoLigPlot, width=width, height=height))

# Create ordering for plots based on first simulation results
orderSimPath <- cbind(ligAttr[[1]], 1:nrow(ligAttr[[1]]))
for(i in (ncol(orderSimPath)-1):1){
  orderSimPath <- orderSimPath[order(orderSimPath[,i], decreasing = T), ]
}

# Move ligands to first position of order
for (i in 1:length(ligNodesNames)) {
  # Find row numbers of unbound forms of ligand
  inds <- grep(ligNodesNames[i], rownames(orderSimPath), fixed=T)
  origInds <- inds

  # Move each of those rows to the beginning of the order (unless they're already at the front)
  for (j in 1:length(origInds)) {
    if (origInds[j] > j) {
      if (inds[j] == 1) { # If not at top
        orderSimPath <- orderSimPath[c(inds[j], 1:(inds[j]-1), (inds[j]+1):nrow(orderSimPath)),]
      } else if (inds[j] == nrow(orderSimPath)) { # If all the way at the bottom
        orderSimPath <- orderSimPath[c(inds[j], 1:(inds[j]-1)),]
      }
      inds <- grep(ligNodesNames[i], rownames(orderSimPath), fixed=T)
    }
  }
}

# Save order
plotOrder <- orderSimPath[,ncol(orderSimPath)]

for (i in 1:rounds) {
  # Reorder the paths
  ligPath[[i]] <- ligPath[[i]][plotOrder,]
  ligAttr[[i]] <- ligAttr[[i]][plotOrder,]
  noLigPath[[i]] <- noLigPath[[i]][plotOrder,]
  noLigAttr[[i]] <- noLigAttr[[i]][plotOrder,]

  # Plot ligand paths/attractors
  plotPath(path = ligPath[[i]], filePath = paste0(outPath, "lig/path/" , i, ".pdf"))
  plotPath(path = ligAttr[[i]], filePath = paste0(outPath, "lig/attractor/", i, ".pdf"))
  write.csv(ligPath[[i]], file = paste0(outPath, "lig/path/" , i, ".csv"))
  write.csv(ligAttr[[i]], file = paste0(outPath, "lig/attractor/" , i, ".csv"))

  # Plot no ligand paths/attractors
  plotPath(path = noLigPath[[i]], filePath = paste0(outPath, "nolig/path/" , i, ".pdf"))
  plotPath(path = noLigAttr[[i]], filePath = paste0(outPath, "nolig/attractor/", i, ".pdf"))
  write.csv(noLigPath[[i]], file = paste0(outPath, "nolig/path/" , i, ".csv"))
  write.csv(noLigAttr[[i]], file = paste0(outPath, "nolig/attractor/" , i, ".csv"))
}

cat("Done.", "\n")