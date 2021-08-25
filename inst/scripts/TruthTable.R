#!/usr/bin/env Rscript

#################################################
# TruthTable.R
# Adrian C
#
# Script to generate a "truth table" for a module
# with all combinations of inhibitors and stimuli.
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
suppressMessages(library(egg))
library(kboolnet)

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
default <- list(modules="", out="./out/", minQuality=0, file=NA, driveFile=NA,
                inputInhibs=NA, inputStimuli=NA, outputs=NA, initState=NA)
opt <- setDefaults(config, default)

# Create out dir if it does not exist
if (!dir.exists(opt$out)) {
  dir.create(opt$out)
}

# Normalize paths
outPath       <- paste0(normalizePath(opt$out), "/")

# Parse modules option to a list
modules <- opt$modules

# Same for inputs options
if (is.na(opt$inputInhibs) & is.na(opt$inputStimuli)) {
  stop("Please provide input nodes to be toggled in simulation rounds")
}
inputInhibs <- opt$inputInhibs
inputStimuli <- opt$inputStimuli

# Same for outputs
if (is.na(opt$outputs)) {
  stop("Please provide output nodes to be toggled in simulation rounds")
}
outputs <- opt$outputs

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
cat("Extracting modules... ")
callExtractModules(masterFile, modulesFile, modules, minQuality)

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

###################### Simulation setup ##################
# Create list of all possible combinations of inputs
numInputs <- length(inputStimuli) + length(inputInhibs)
rounds <- 2^numInputs
combinations <- matrix(nrow = rounds, ncol = numInputs)
for (i in 0:(rounds-1)) {
  combinations[i + 1,] <- as.logical(intToBits(i))[1:numInputs]
}
combinationsStimuli <- combinations[,1:length(inputStimuli), drop=FALSE]
if(length(inputInhibs) > 0) {
  combinationsInhibs <- combinations[,(length(inputStimuli)+1):ncol(combinations), drop=FALSE]
}

# Make sure regexes match something in network
inputStimuliRegex <- escapeRegex(inputStimuli)
for (i in 1:length(inputStimuli)) {
  matches <- symbolMapping$name[grepl(inputStimuliRegex[i], symbolMapping$name)]
  if(length(matches) == 0) {
    stop("Input ", inputStimuli[i], " could not be matched to a rxncon node.")
  } else {
    cat("Input", inputStimuli[i], "was matched to rxncon nodes", paste0(matches, collapse=", "), "\n")
  }
}

if (length(inputInhibs) > 0) {
  inputInhibsRegex <- escapeRegex(inputInhibs)
  for (i in 1:length(inputInhibs)) {
    matches <- symbolMapping$name[grepl(inputInhibsRegex[i], symbolMapping$name)]
    if(length(matches) == 0) {
      stop("Input ", inputInhibs[i], " could not be matched to a rxncon node.")
    } else {
      cat("Input", inputInhibs[i], "was matched to rxncon nodes", paste0(matches, collapse=", "), "\n")
    }
  }
}

outputsRegex <- escapeRegex(outputs)
for(i in 1:length(outputs)) {
  matches <- symbolMapping$name[grepl(outputsRegex[i], symbolMapping$name)]
  if(length(matches) == 0) {
    stop("Output ", outputs[i], " could not be matched to a rxncon node.")
  } else {
    cat("Output", outputs[i], "was matched to rxncon nodes", paste0(matches, collapse=", "), "\n")
  }
}

####################### Simulate #########################
# Create matrix to store simulation results
results <- matrix(nrow = rounds, ncol = length(outputs), dimnames = list(1:rounds, outputs))
attractors <- list()

# Simulation rounds
cat("Simulating", rounds, "activation/inhibition combinations of", ncol(combinations), "inputs...", "\n")
pb <- txtProgressBar(min = 0, max = rounds, style = 3)
for (i in 1:rounds) {
  roundNetwork <- network

  # Fix stimulus input nodes present in round on
  onNodes <- character()
  for (regex in inputStimuliRegex[combinationsStimuli[i,]]) {
    onNodes <- c(onNodes, symbolMapping$ID[grepl(regex, symbolMapping$name)])
  }
  roundNetwork <- fixGenes(roundNetwork, onNodes, 1)

  # Fix stimulus nodes not present in round off
  offNodes <- character()
  for (regex in inputStimuliRegex[!combinationsStimuli[i,]]) {
    offNodes <- c(offNodes, symbolMapping$ID[grepl(regex, symbolMapping$name)])
  }

  # Fix inhibited nodes present in round off
  if (length(inputInhibs > 0)) {
    for (regex in inputInhibsRegex[combinationsInhibs[i,]]) {
      offNodes <- c(offNodes, symbolMapping$ID[grepl(regex, symbolMapping$name)])
    }
  }
  roundNetwork <- fixGenes(roundNetwork, offNodes, 0)

  # Simulate
  attractor <- getPathAndAttractor(roundNetwork, initStates$state, initStates$name)$attractor

  # Save data
  for (j in 1:length(outputs)) {
    outNodes <- symbolMapping$name[grepl(outputsRegex[j], symbolMapping$name)]
    results[i, outputs[j]] <- mean(as.matrix(attractor[outNodes,,drop=FALSE]))
  }
  attractors[[i]] <- attractor

  setTxtProgressBar(pb, i)
}
close(pb)

##################### Plot ##############################
cat("Generating plots... ")
# Try and simplify conditions that result in the same outputs
results <- as.data.frame(results)
combinations <- as.data.frame(combinations)
newCombinations <- combinations[0,,drop=F]
newResults <- results[0,,drop=F]
for (i in 1:ncol(combinations)) { # For each input
  # Get the combinations in which which input is present
  inputPresent <- which(combinations[,i])

  # If there is no variation in the output value when the input is present
  identicalCols <- apply(results[inputPresent,], 2, function(x) {(max(x) - min(x)) == 0})
  if (all(identicalCols)) {
    # Append it to newCombinations and newResults
    newResults[nrow(newResults) + 1,] <- results[inputPresent[1],,]
    newCombinations[nrow(newCombinations) + 1,] <- rep("*", ncol(combinations)) # Set all inputs to *
    newCombinations[nrow(newCombinations), i] <- TRUE # Except for the input that's identical
  }
}

# Replace redundant combinations with simplified version
if (nrow(newCombinations) > 0) {
  for (i in 1:nrow(newCombinations)) {
    delete <- combinations[,which(newCombinations[i,] == TRUE)] # Get the redundant rows to delete
    combinations <- combinations[!delete,,drop=F] # Remove the redundant combinations
    results <- results[!delete,,drop=F] # Remove the redundant results
  }
  combinations <- rbind(combinations, newCombinations)
  results <- rbind(results, newResults)
}

# Gather the dataframes so they're in plottable format
resultsGather <- results
resultsGather$num <- 1:nrow(results)
resultsGather <- gather(resultsGather, "signal", "value", -num)
resultsGather$color[resultsGather$value > 0.5] <- "white"
resultsGather$color[resultsGather$value <= 0.5] <- "black"
resultsGather$signal <- factor(resultsGather$signal, levels = outputs) # This puts the outputs in the order inputted

# Do the same for inputs
if (length(inputInhibs) > 0) {
  colnames(combinations) <- c(inputStimuli, inputInhibs)
} else {
  colnames(combinations) <- inputStimuli
}
combinationsGather <- combinations
combinationsGather$num <- 1:nrow(combinations)
combinationsGather <- gather(combinationsGather, "input", "value", -num)

# Add type information to inputs
combinationsGather$value[combinationsGather$value == FALSE] = "None"
combinationsGather$value[(combinationsGather$input %in% inputStimuli) & combinationsGather$value == TRUE] = "Stimulus"
combinationsGather$value[(combinationsGather$input %in% inputInhibs) & combinationsGather$value == TRUE] = "Inhibitor"
combinationsGather$input <- factor(combinationsGather$input, levels = c(inputStimuli, inputInhibs)) # Puts inputs in order inputted

# Inputs plot
inputColors <- c("*"="white", "None"="white", "Inhibitor"="red2", "Stimulus"="green3")
inputPlot <- ggplot(combinationsGather, aes(x=.5, y=.5)) +
  facet_grid(cols=vars(input), rows=vars(num)) +
  geom_tile(aes(fill=value), colour="black", size=1) + # Create the tiles (colour and size affect borders)
  scale_fill_manual(values=inputColors, drop=FALSE) + # Color the tiles appropriately
  scale_x_continuous(limits=c(0,1), expand=c(0,0)) + scale_y_reverse(limits=c(1,0), expand=c(0,0)) + # Reverse the ordering and set proper scales
  theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), # Remove axis labels
        legend.position = "none", panel.background = element_blank(), strip.text.y = element_blank(), # Remove legend, background, and y labels
        strip.background = element_blank(), strip.text.x = element_text(angle=90, hjust=0), # Turn cue labels 90 degrees
        plot.margin = unit(c(0,5.5,0,0), "pt"), plot.background = element_blank()) # Remove margins, these are set by data plot
if(nrow(combinationsGather[combinationsGather$value == "*",]) > 0) { # If necessary, add stars to represent redundant inputs
  inputPlot <- inputPlot + geom_text(data=combinationsGather[combinationsGather$value == "*",], label="*", size=8, vjust=.7)
}

# Results plot
resultsPlot <- ggplot(resultsGather, aes(x=.5, y=.5)) +
  facet_grid(cols=vars(signal), rows=vars(num)) +
  geom_tile(aes(fill=value), colour="black", size=1) + # Create the tiles (colour and size affect borders)
  geom_text(aes(color=color, label=format(value, nsmall=2, digits=2))) + # Label the values on the plot
  scale_color_manual(values=c("black"="black", "white"="white"), guide="none") + # Color the labels
  scale_fill_gradient(limits=c(0,1), low="white", high="steelblue") + # Color the tiles appropriately
  scale_x_continuous(limits=c(0,1), expand=c(0,0)) + scale_y_reverse(limits=c(1,0), expand=c(0,0)) + # Reverse the ordering and set proper scales
  theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), # Remove axis labels
        panel.background = element_blank(), strip.text.y = element_blank(), legend.position = "none", # Remove legend, background, and y labels
        strip.background = element_blank(), strip.text.x = element_text(angle=90, hjust=0), # Turn cue labels 90 degrees
        plot.margin = unit(c(0,5.5,0,0), "pt"), plot.background = element_blank()) # Remove margins, these are set by data plot

widths = c(ncol(combinations), ncol(results)*3)
p <- suppressWarnings(suppressMessages(ggarrange(inputPlot, resultsPlot, widths = widths, nrow = 1, ncol = 2)))

height = (.7 * nrow(results)) + 1
width =  (0.25 * (ncol(combinations) + ncol(results) * 3)) + 1
ggsave(plot=p, device="pdf", paste0(outPath, "truth_table.pdf"), height = height, width = width)

cat("Done.", "\n")

################### Save data ###############
cat("Saving data... ")

# Add stimulus/inhibitors data to results
results$num <- NULL
for (i in 1:nrow(results)) {
  results$stimuli[i] <- paste0(inputStimuli[combinationsStimuli[i,]], collapse = ", ")
  if (length(inputInhibs) > 0) {
    results$inhibitors[i] <- paste0(inputInhibs[combinationsInhibs[i,]], collapse = ", ")
  }
}

# Write to csv
write.csv(results, file = paste0(outPath, "truth_table.csv"), row.names = FALSE)

cat("Done.", "\n")
