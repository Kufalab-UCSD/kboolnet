#!/usr/bin/env Rscript

#################################################
# ScoreNet.R
# Adrian C
#
# Script to download model and experimental database 
# from Google Drive and generate comparison score. 
#################################################

################# Library loading ##################
options(stringsAsFactors = F)
options(warn = 1)
suppressMessages(library(BoolNet))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(openxlsx))
suppressMessages(library(googledrive))
suppressMessages(library(optparse))
suppressMessages(library(tidyr))
suppressMessages(library(egg))

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
  make_option("--rxnconFile", action="store", default=NA, type="character",
              help="Path to master rxncon file (local)"),
  make_option("--rxnconDriveFile", action="store", default=NA, type="character",
              help="File name or path of master rxncon file (on Google Drive)"),
  make_option("--MIDASFile", action="store", default=NA, type="character",
              help="Path to MIDAS file (local)"),
  make_option("--MIDASDriveFile", action="store", default=NA, type="character",
              help="File name or path of MIDAS file (on Google Drive)"),
  make_option(c("--modules", "-m"), action="store", default=NA, type="character",
              help="Comma-separated modules to be loaded from master rxncon file [default: load all modules]"),
  make_option("--minQuality", action="store", default=NA, type="integer",
              help="Minimum quality for rule to be loaded [default: 0]"),
  make_option(c("--out", "-o"), action="store", default=NA, type="character",
              help="Folder to which output files will be written [default: ./out/]"),
  make_option(c("--bin", "-b"), action="store", default=NA, type="integer",
              help="Bin to use as final timepoint when scoring. 0 = All data after t = 0, 1 = (0, 3], 2 = (3, 30], 3 = (30, 300], 4 = (300, Inf) [default: 0]"),
  make_option(c("--normalize", "-n"), action="store_true", default=NA, type="logical",
              help="Normalize output of simulations such that the maximum value for a measurement across allsimulation is 1."),
  make_option(c("--pretreat", "-p"), action="store_true", default=NA, type="logical",
              help="Option to cause inhibitors to be applied BEFORE stimuli are applied"),
  make_option(c("--celltype"), action="store", default=NA, type="character",
              help="Comma-separated list of cell types to filter original data by"),
  make_option(c("--height"), action="store", default=NA, type="character",
              help="Height in inches of the plots"),
  make_option(c("--width"), action="store", default=NA, type="character",
              help="Width in inches of the plots")
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
default <- list(modules="", out="./out/", minQuality=0, rxnconFile=NA, rxnconDriveFile=NA, width=NA, height= NA,
                MIDASFile=NA, MIDASDriveFile=NA, bin=0, normalize=FALSE, pretreat=FALSE, celltype="all")
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
suppressMessages(source(paste0(kboolnetPath, "R/functions/driveDownload.R")))
suppressMessages(source(paste0(kboolnetPath, "R/functions/fixedNetwork.R")))
suppressMessages(source(paste0(kboolnetPath, "R/functions/readMIDASExcel.R")))
suppressMessages(source(paste0(kboolnetPath, "R/functions/plotMIDAS.R")))
suppressMessages(source(paste0(kboolnetPath, "R/functions/plotMIDASComp.R")))
suppressMessages(source(paste0(kboolnetPath, "R/functions/getPathAndAttractor.R")))

# Parse modules option to a list
modules <- trimws(strsplit(opt$modules, ",")[[1]])
cell_types <- trimws(strsplit(opt$celltype, ",")[[1]])

# Make sure input file paths were provided
if (is.na(opt$rxnconFile) & is.na(opt$rxnconDriveFile)){ # If neither file was provided
  stop("Please provide a path to a local rxncon file with --rxnconFile or a Google Drive file with --rxnconDriveFile")
} else if ((!is.na(opt$rxnconFile)) & (!is.na(opt$rxnconDriveFile))) { # If both files were provided
  stop("Please provide only one path with EITHER --rxnconFile or --rxnconDriveFile")
} else if (is.na(opt$MIDASFile) & is.na(opt$MIDASDriveFile)){ # If neither file was provided
  stop("Please provide a path to a local MIDAS file with --rxnconFile or a Google Drive file with --rxnconDriveFile")
} else if ((!is.na(opt$MIDASFile)) & (!is.na(opt$MIDASDriveFile))) { # If both files were provided
  stop("Please provide only one path with EITHER --MIDASFile or --MIDASDriveFile")
} 

minQuality <- opt$minQuality
if (opt$minQuality < 0) {
  stop("minQuality must be >= 0")
}

################ Load and process rxncon file ###################

# If GDrive file provided
if (!(is.na(opt$rxnconDriveFile))) {
  cat("Downloading rxncon file from Google Drive... ", "\n")
  
  # Download file
  masterFile  <- paste0(outPath, "master.xlsx")
  driveDownload(driveFile = opt$rxnconDriveFile, out = masterFile, type = "spreadsheet")
} else { # If local file provided
  masterFile <- normalizePath(opt$rxnconFile)
  
  # File verification
  if (!grepl("\\.xlsx$", masterFile)) { # Make sure file is Excel file
    stop("rxncon file must be an Excel file (.xslx extension)")
  } else if (!file.exists(masterFile)) { # Make sure file exists
    stop("rxncon file ", masterFile, " does not exist")
  }
}

# Extract modules from master file, write to modules file
modulesFile <- paste0(outPath, "modules.xlsx")
cat("Extracting modules... ")
suppressWarnings(stderr <- system2(command = "python3", args = c(paste0(kboolnetPath, "Python/extract_modules.py"), "--file", masterFile,
                                                                 "--modules", paste0('"', paste0(modules, collapse=","), '"'), "--quality", minQuality,
                                                                 "--output", modulesFile), stderr = TRUE, stdout = ""))
if (any(grepl("Error", stderr, ignore.case = TRUE))) {
  cat(paste(stderr, "\n"))
  stop("Error during module extraction. Please run extract_modules.py on its own with the -v DEBUG flag.")
}
cat("Done.", "\n")

# Pass files to rxncon for processing
cat("Generating BoolNet files... ")
netFilePrefix <- gsub("\\.xlsx$", "", modulesFile)
suppressWarnings(stderr <- system2("python3", args = c(paste0(rxnconPath, "rxncon2boolnet.py"), modulesFile, "--output",
                                                       netFilePrefix), stderr = TRUE, stdout = ""))
if (any(grepl("Error", stderr, ignore.case = TRUE))) {
  cat(paste(stderr, "\n"))
  stop("Error during BoolNet file generation. Please run rxncon2boolnet.py on its own with the -v DEBUG flag.")
}
cat("Done.", "\n")

############### Load and process MIDAS file ################
# If GDrive file provided
if (!(is.na(opt$MIDASDriveFile))) {
  cat("Downloading MIDAS file from Google Drive...", "\n")
  
  # Download file
  MIDASFile  <- paste0(outPath, "MIDAS.xlsx")
  driveDownload(driveFile = opt$MIDASDriveFile, out = MIDASFile, type = "spreadsheet")
} else { # If local file provided
  MIDASFile <- normalizePath(opt$MIDASFile)
  
  # File verification
  if (!grepl("\\.xlsx$", MIDASFile)) { # Make sure file is Excel file
    stop("MIDAS file must be an Excel file (.xslx extension)")
  } else if (!file.exists(MIDASFile)) { # Make sure file exists
    stop("MIDAS file ", MIDASFile, " does not exist")
  }
}

# Read the MIDAS excel into a MIDASlist
cat("Processing MIDAS file...", "\n")
MIDASlist <- readMIDASExcel(MIDASFile, cell_types)
cat("Done.", "\n")

################# BoolNet setup ######################
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

# Try to match all treatments to rxncon nodes. If this fails, remove that treatment and all
# places where that treatment was used
removeTreatments <- numeric()
for (i in 1:length(MIDASlist$namesCues)) {
  treatment <- MIDASlist$treatmentDefs[MIDASlist$treatmentDefs$name == MIDASlist$namesCues[i],] # Get the treatment

  # Try mapping regex to nodes
  for (j in 1:length(treatment$regex[[1]])) {
    if(!any(grepl(paste0("^", treatment$regex[[1]][j], "(_.*--0|_.*-\\{0\\}|)$"), initStates$name))) { # If no nodes matched
      warning("Treatment ", treatment$name, " (node(s): ", paste0(treatment$nodes[[1]][j], collapse=", "),
              ") could not be matched to nodes in the rxncon system, removing data points involving this treatment.")

      # Remove rows from MIDAS list that used this treatment
      cueRows <- which(as.logical(MIDASlist$valueCues[,i]))
      if (length(cueRows) > 0) {
        MIDASlist$valueCues <- MIDASlist$valueCues[-cueRows,]
        MIDASlist$valueSignals <- MIDASlist$valueSignals[-cueRows,,]
        MIDASlist$valueVariances <- MIDASlist$valueVariances[-cueRows,,]
      }

      # Add this treatment to list of treatments to remove from treatmentDefs/namesCues
      removeTreatments <- c(removeTreatments, i)
    }
  }
}
if(length(removeTreatments) > 0) {
  MIDASlist$valueCues <- MIDASlist$valueCues[,-removeTreatments]
  MIDASlist$treatmentDefs <- MIDASlist$treatmentDefs[!(MIDASlist$treatmentDefs$name %in% MIDASlist$namesCues[removeTreatments]),]
  MIDASlist$namesCues <- MIDASlist$namesCues[-removeTreatments]
}

# Do the same as above but for measurements
removeMeasurements <- numeric()
for (i in 1:length(MIDASlist$namesSignals)) {
  if(!any(grepl(MIDASlist$regexSignals[i], symbolMapping$name))) {
    warning("Measurement ", MIDASlist$namesSignals[i], " could not be matched to a node in the rxncon system, ignoring.")
    removeMeasurements <- c(removeMeasurements, i)
  }
}
# Remove the measurements tagged for removal
if(length(removeMeasurements) > 0) {
  MIDASlist$valueSignals <- MIDASlist$valueSignals[,-removeMeasurements,]
  MIDASlist$valueVariances <- MIDASlist$valueVariances[,-removeMeasurements,]
  MIDASlist$namesSignals <- MIDASlist$namesSignals[-removeMeasurements]
  MIDASlist$regexSignals <- MIDASlist$regexSignals[-removeMeasurements]
}

# Get all ligand and mutant nodes/components
ligNodes <- character()
for (i in 1:nrow(MIDASlist$treatmentDefs)) {
  if (grepl("(Stimulus|Mutant)", MIDASlist$treatmentDefs$type[i])) {
    if (!is.na(MIDASlist$treatmentDefs$regex[[i]])) ligNodes <- c(ligNodes, MIDASlist$treatmentDefs$regex[[i]])
  }
}

# Set them to be off in initial states 
for (lig in ligNodes) {
  initStates$state[grepl(paste0("^", lig, "(_.*--0|_.*-\\{0\\}|)$"), initStates$name)] <- 0
}

# Map namesSignals to rxncon nodes
signalMapping <- list()
for (i in 1:length(MIDASlist$namesSignals)) {
  signalNodes <- symbolMapping$name[grepl(MIDASlist$regexSignals[i], symbolMapping$name)]
  cat("Signal", MIDASlist$namesSignals[i], "mapped to node(s)", paste0(signalNodes, collapse=", "), "\n")
  signalMapping <- c(signalMapping, list(signalNodes))
}

####################### Simulate #########################
# Create MIDAS list to store simulation results
simMIDASlist <- MIDASlist
simMIDASlist$timeSignals <- c(0, 1)
emptyArray <- array(dim = c(dim(MIDASlist$valueSignals)[1], dim(MIDASlist$valueSignals)[2], 2), # Empty array with same dimensions as original, except only 2 timepoints
                    dimnames = list(dimnames(MIDASlist$valueSignals)[[1]], dimnames(MIDASlist$valueSignals)[[2]], c(0, 1)))
simMIDASlist$valueSignals <- emptyArray
simMIDASlist$valueVariances <- emptyArray

# Simulation rounds is equal to number of cue combinations
rounds <- nrow(MIDASlist$valueCues)

# Set up progress bar
cat("Simulating", rounds, "different experimental conditions...", "\n")

# Set up variables to store simulation results
scores <- data.frame(init=numeric(rounds), final=numeric(rounds))

# Simulate
for (i in 1:rounds) {
  # Get all the cues applied in the round
  roundCues <- MIDASlist$namesCues[as.logical(MIDASlist$valueCues[i,])]
  
  # For each cue, sort it into the correct type
  roundLigs <- character()
  roundKOs <- character()
  roundInhibs <- character()
  roundMutants <- character()
  for (cue in roundCues) {
    cueType <- MIDASlist$treatmentDefs$type[MIDASlist$treatmentDefs$name == cue]
    cueRegex <- MIDASlist$treatmentDefs$regex[MIDASlist$treatmentDefs$name == cue][[1]]
    if (cueType == "Stimulus") roundLigs <- c(roundLigs, cueRegex)
    if (cueType == "KO") roundKOs <- c(roundKOs, cueRegex)
    if (cueType == "Inhibitor") roundInhibs <- c(roundInhibs, cueRegex)
    if (cueType == "Mutant") roundMutants <- c(roundMutants, cueRegex)
  }
  
  cat(paste0("Round ", i, " of ", rounds, ": simulating with KOs: ", paste0(roundKOs, collapse=", "), "; inhibitors: ",
             paste0(roundInhibs, collapse=", "), "; ligands: ", paste0(roundLigs, collapse=", "), "; mutants: ", paste0(roundMutants, collapse=", ")), "\n")
  
  # Apply KOs
  KOnetwork <- network
  if (length(roundKOs) > 0) {
    KOnetwork <- fixedNetwork(KOnetwork, roundKOs, symbolMapping$name, symbolMapping$ID, 0)
  }
  
  # Add mutants
  if (length(roundMutants) > 0) {
    KOnetwork <- fixedNetwork(KOnetwork, roundMutants, symbolMapping$name, symbolMapping$ID, 1)
  }
  
  # Initial simulation to neutral state (t = 0)
  neutralAttr <- getPathAndAttractor(KOnetwork, initStates$state, symbolMapping$name)$attractor
  
  # Inhibit inhibitor nodes
  if (length(roundInhibs) > 0) {
    inhibNetwork <- fixedNetwork(KOnetwork, roundInhibs, symbolMapping$name, symbolMapping$ID, 0)
  } else {
    inhibNetwork <- KOnetwork
  }
  
  # If pretreatment desired, simulated to inhibited attractor
  if (opt$pretreat) {
    neutralAttr <- getPathAndAttractor(inhibNetwork, initStates$state, symbolMapping$name)$attractor
  }
  
  # Use last neutral state as a new initial state
  newInitStates <- initStates
  newInitStates$state <- neutralAttr[,ncol(neutralAttr)]

  # Add ligands
  if (length(roundLigs) > 0) {
    for (lig in roundLigs) {
      newInitStates$state[grepl(paste0(lig, "(_.*--0$|_.*-\\{0\\}|$)"), initStates$name)] <- 1
    }
  }
  
  # Simulate to final state
  finalAttr <- getPathAndAttractor(inhibNetwork, newInitStates$state, symbolMapping$name)$attractor
  
  # Save simulation data to simMIDASlist
  for (j in 1:length(MIDASlist$namesSignals)) {
    # Average all the nodes for a signal
    neutralAttrMean <- mean(rowMeans(neutralAttr[signalMapping[[j]],,drop=F]))
    finalAttrMean <- mean(rowMeans(finalAttr[signalMapping[[j]],,drop=F]))
    
    # Put these values in the corresponding spot for simMIDASlist
    simMIDASlist$valueSignals[i,j,1] <- neutralAttrMean
    simMIDASlist$valueSignals[i,j,2] <- finalAttrMean
  }
}

cat("Simulations complete after", rounds, "rounds.", "\n")

################# Simulation normalizataion ##############
# Normalize simulation results if requested
# Divides all simulation results by the maximum value for that output
if (opt$normalize) {
  for (i in 1:ncol(simMIDASlist$valueSignals)) {
    maxVal <- max(simMIDASlist$valueSignals[,i,])
    if (maxVal > 0) {
      simMIDASlist$valueSignals[,i,] <- simMIDASlist$valueSignals[,i,] / maxVal
    }
  }
}

###################### Binning ###########################
# Bin the data if binning was requested (i.e. --bin != 0)
if (opt$bin != 0) {
  # Create bins and bin names/times
  bins <- list(c(-Inf, 0), c(0,3), c(3, 30), c(30, 300), c(300, Inf)) # These are the bins; for bin c(a, b) timepoints a < t <= b will be averaged
  binNames <- sapply(bins, paste0, collapse=" < t <= ")
  binTimes <- numeric(length(bins)) # This is for plotting purposes
  for (i in 1:length(bins)) {
    bin <- bins[[i]]
    if (bin[1] == -Inf) { # If bin goes from negative infinity to a point, set that point as the time 
      binTimes[i] <- bin[2]
    } else if (bin[2] == Inf) { # If bin goes from point to infinity, set point as the time
      binTimes[i] <- bin[1]
    } else { # Otherwise, set time as mean of bin points
      binTimes[i] <- mean(bin)
    }
  }
  binnedData <- array(dim = c(nrow(MIDASlist$valueCues), length(MIDASlist$namesSignals), length(bins)), # Create array to put binned data together
                      dimnames = list(1:nrow(MIDASlist$valueCues), MIDASlist$namesSignals, binNames))
  
  # Create MIDASlist to store the binned data
  binMIDASlist <- MIDASlist
  binMIDASlist$valueVariances <- binnedData
  binMIDASlist$timeSignals <- binTimes
  
  cat("Binning data using bins:", paste0(binNames, collapse=", "), "\n")
  cat(paste0("Bin #", opt$bin, " (", binNames[opt$bin+1], ") will be used for scoring."), "\n")
  for (i in 1:length(bins)) {
    bin <- bins[[i]]
    binTimes <- MIDASlist$timeSignals > bin[1] & MIDASlist$timeSignals <= bin[2] # Get all the time signals for the bin
    binArray <- MIDASlist$valueSignals[,,binTimes, drop=F] # Get the data of the bin
    binnedData[,,i] <- apply(binArray, 1:2, mean, na.rm = TRUE) # Average and store to binnedData
  }
  binnedData[is.nan(binnedData)] <- NA
  binMIDASlist$valueSignals <- binnedData
  simMIDASlist$timeSignals[2] <- binMIDASlist$timeSignals[opt$bin+1] # Set the simulation time to the bin time
  
  # Save only the bins being used for comparison to avgMIDASlist
  avgMIDASlist <- MIDASlist
  avgMIDASlist$timeSignals <- binMIDASlist$timeSignals[c(1, opt$bin+1)]
  avgMIDASlist$valueVariances <- binMIDASlist$valueVariances[,,c(1, opt$bin+1)]
  avgMIDASlist$valueSignals <- binMIDASlist$valueSignals[,,c(1, opt$bin+1)]
} else { # If binning not requested, average all data after t = 0 into a single bin
  cat("Averaging all data after t = 0 for use in scoring.", "\n")
  avgMIDASlist <- simMIDASlist
  avgMIDASlist$valueSignals[,,1] <- MIDASlist$valueSignals[,,1]
  avgMIDASlist$valueSignals[,,2] <- apply(MIDASlist$valueSignals[,,2:length(MIDASlist$timeSignals),drop=FALSE], 1:2, mean, na.rm = TRUE)
  avgMIDASlist$valueSignals[is.nan(avgMIDASlist$valueSignals)] <- NA # Turn all NaNs into NAs
}

####################### Error calculation ###########################
# Calculate squared errors and create new MIDASlist containing them
initErr <- (avgMIDASlist$valueSignals[,,1] - simMIDASlist$valueSignals[,,1])^2
finalErr <- (avgMIDASlist$valueSignals[,,2] - simMIDASlist$valueSignals[,,2])^2
avgErr <- apply(array(c(unlist(initErr), unlist(finalErr)), dim=c(nrow(initErr), ncol(initErr), 2)), 1:2, mean, na.rm = TRUE)
avgErr[is.nan(avgErr)] <- NA # Turn all NaNs into NAs

cat("Mean square error for initial timepoint: ", mean(initErr, na.rm = TRUE), "\n")
cat("Mean square error for final timepoint: ", mean(finalErr, na.rm = TRUE), "\n")
cat("Mean square error for all data:", mean(c(unlist(initErr), unlist(finalErr)), na.rm = TRUE), "\n")

##################### Plot and save data ##############################
# Plot setup
cat("Generating plots... ")
if (is.na(opt$height)) {
  height = 0.9 * nrow(MIDASlist$valueCues)
} else {
  height = opt$height
}
if (is.na(opt$width)) {
  width = 1.0 * (length(MIDASlist$namesSignals) + (sum(colSums(MIDASlist$valueCues) > 0)/4))
} else {
  width = opt$width
}

# Experimental data plot
p <- plotMIDAS(MIDASlist)
ggsave(paste0(outPath, "Experimental_MIDAS.pdf"), p, height = height, width = width)

# Either binned plot or averaged plot
if (opt$bin != 0) {
  p <- plotMIDAS(binMIDASlist)
  ggsave(paste0(outPath, "Experimental_Binned_MIDAS.pdf"), p, height = height, width = width)
} else {
  p <- plotMIDAS(avgMIDASlist)
  ggsave(paste0(outPath, "Experimental_Averaged_MIDAS.pdf"), p, height = height, width = width)
}

# Simulation data plot
p <- plotMIDAS(simMIDASlist)
ggsave(paste0(outPath, "Simulation_MIDAS.pdf"), p, height = height, width = width)

# Comparison plot
p <- plotMIDASComp(avgMIDASlist, simMIDASlist, avgErr)
ggsave(paste0(outPath, "Comparison_MIDAS.pdf"), p, height = height, width = width)
cat("Done.", "\n")


# Save MIDAS lists in RDS files
cat("Saving data... ")
saveRDS(MIDASlist, file = paste0(outPath, "Experimental_MIDAS.rds"))
if (opt$bin != 0) saveRDS(binMIDASlist, file = paste0(outPath, "Experimental_Binned_MIDAS.rds"))
if (opt$bin == 0) saveRDS(avgMIDASlist, file = paste0(outPath, "Experimental_Averaged_MIDAS.rds"))
saveRDS(simMIDASlist, file = paste0(outPath, "Simulation_MIDAS.rds"))

cat("Done.", "\n")
