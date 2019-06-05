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

################# Argument parsing #################
# Get commandline args
option_list = list(
  make_option("--config", action="store", default=NA, type="character",
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
  make_option("--out", action="store", default=NA, type="character",
              help="Folder to which output files will be written [default: ./out/]"),
  make_option("--ligands", action="store", default=NA, type="character",
              help="Comma-separated rxncon name(s) of ligand node(s) to be toggled in verification simulation"),
  make_option("--rounds", action="store", default=NA, type="integer",
              help="Maximum number of ligand/no-ligand simulation pairs that should be run (default: 20)")
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
default <- list(modules="", out="./out/", minQuality=0, ligands=NA, file=NA, driveFile=NA, rounds=20)
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
suppressMessages(source(paste0(kboolnetPath, "functions/extractModules.R")))
suppressMessages(source(paste0(kboolnetPath, "functions/plotPath.R")))
suppressMessages(source(paste0(kboolnetPath, "functions/unbindLigand.R")))
suppressMessages(source(paste0(kboolnetPath, "functions/compMatrix.R")))

# Parse modules option to a list
modules <- trimws(strsplit(opt$modules, ",")[[1]])

# Same for ligands option
if (is.na(opt$ligands)) {
  stop("Please provide ligand(s) to be toggled in simulation rounds")
}
ligands <- trimws(strsplit(opt$ligands, ",")[[1]])

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

maxrounds <- opt$maxrounds
if (opt$maxrounds < 2) {
  stop("maxrounds must be >= 2")
}

################ Load and process rxncon file ###################

# If GDrive file provided
if (!(is.na(opt$driveFile))) {
  cat("Downloading rxncon file from Google Drive...", "\n")
  
  if (grepl("^https?:\\/\\/", opt$driveFile)) { # If URL provided
    # Deactivate authentication and try and access file 
    drive_auth_config(active = FALSE)
    
    # Try accessing file without authentication
    gDriveID <- as_id(tryCatch({
      drive_get(id = as_id(opt$driveFile))$id[1]
    }, error = function(e) {
      # If access fails, ask user if they want to try again with authentication
      repeat {
        i <- readline(prompt = "Public Google Drive file not found. Try again with authentication? (y)es/(n)o: ")
        
        if (grepl("^(?:yes|y)$", i, ignore.case = TRUE)) {
          drive_auth_config(active = TRUE) # Activate authentication
          
          return(tryCatch({ # Try to find file again
            drive_get(id = as_id(opt$driveFile))$id[1]
            
          }, error = function(e){ # If it fails again
            stop("Google Drive file not found. Either the file does not exist or you do not have permission to view it.")
          }))
        } else if (grepl("^(?:no|n)$", i, ignore.case = TRUE)){
          stop("Please provide an input file for rxncon.")
        }
      }
    }))
  } else { # Search for the file in Drive if name provided
    # Activate Google Drive authentication, required for searching in a user's drive
    drive_auth_config(active = TRUE)
    
    gDriveID <- as_id(drive_find(pattern = opt$driveFile, type = "spreadsheet")$id[1])
    
    if(is.na(gDriveID)) { # If file does not exist
      stop("File not found in Google Drive.")
    }
  }
  
  # Download file
  masterFile  <- paste0(outPath, "master.xlsx")
  drive_download(gDriveID, path = masterFile, type = "xlsx", overwrite = T) # Download the file
  
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
# Lists to store simulation results
ligAttr   <- list()
ligPath   <- list()
noLigAttr <- list()
noLigPath <- list()

# Matrices to store comparison results
scoreLig   <- matrix(nrow = maxrounds, ncol = maxrounds)
scoreNoLig <- matrix(nrow = maxrounds, ncol = maxrounds)

# Simulation rounds
cat("Starting simulations...", "\n")
rounds <- 0
for (i in 1:maxrounds) {
  rounds <- rounds + 1
  # Remove ligands 
  for (lig in 1:length(ligands)){
    initStates$state <- unbindLigand(ligands[lig], initStates$name, initStates$state) # Unbind ligand from complexes
    initStates$state[grepl(paste0(ligands[lig], "_.*--0$"), initStates$name)] <- 0 # Remove unbound ligand
  }
  
  # Simulate w/out ligand
  noLigPath[[i]]    <- getPathToAttractor(network, initStates$state) %>% t()
  rownames(noLigPath[[i]]) <- initStates$name
  initStates$state  <- noLigPath[[i]][,ncol(noLigPath[[i]])] # Get only attractor of no lig path
  noLigAttr[[i]]    <- getPathToAttractor(network, initStates$state) %>% t()
  rownames(noLigAttr[[i]]) <- initStates$name
  
  # Add ligands
  for (lig in 1:length(ligands)){
    initStates$state[grepl(paste0(ligands[lig], "_.*--0$"), initStates$name)] <- 1
  }
  
  # Simulate w/ ligand
  ligPath[[i]]      <- getPathToAttractor(network, initStates$state) %>% t()
  rownames(ligPath[[i]]) <- initStates$name
  initStates$state  <- ligPath[[i]][,ncol(ligPath[[i]])] # Get only attractor of no lig path
  ligAttr[[i]]      <- getPathToAttractor(network, initStates$state) %>% t()
  rownames(ligAttr[[i]]) <- initStates$name
  
  # Compare this round of simulation to previous rounds
  for (j in 1:rounds) {
    scoreLig[j,rounds] <- compMatrix(ligAttr[[rounds]], ligAttr[[j]])
    scoreNoLig[j,rounds] <- compMatrix(noLigAttr[[rounds]], noLigAttr[[j]])
  }
  
  # If both ligand and no-ligand simulations are identical to a previous round, then we know a
  # "meta-attractor" has been found and we can stop simulation
  if (rounds != 1 & any(scoreLig[1:rounds-1,rounds] == 1) & any(scoreNoLig[1:rounds-1,rounds] == 1)) {
    cat("Meta-attractor found! Stopping simulation.", "\n")
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
# Create output subdirs if they do not exist
if (!dir.exists(paste0(outPath, "lig/"))) dir.create(paste0(outPath, "lig/"))
if (!dir.exists(paste0(outPath, "lig/attractor/"))) dir.create(paste0(outPath, "lig/attractor/"))
if (!dir.exists(paste0(outPath, "lig/path/"))) dir.create(paste0(outPath, "lig/path/"))
if (!dir.exists(paste0(outPath, "nolig/"))) dir.create(paste0(outPath, "nolig/"))
if (!dir.exists(paste0(outPath, "nolig/attractor/"))) dir.create(paste0(outPath, "nolig/attractor/"))
if (!dir.exists(paste0(outPath, "nolig/path/"))) dir.create(paste0(outPath, "nolig/path/"))

# Write comparison data 
write.csv(scoreLig, file = paste0(outPath, "lig/simulation_comparison.csv"))
write.csv(scoreNoLig, file = paste0(outPath, "nolig/simulation_comparison.csv"))

# Plot the attractor comparison matrices
scoreLigPlot   <- plotCompMat(scoreLig)
scoreNoLigPlot <- plotCompMat(scoreNoLig)
suppressMessages(ggsave(paste0(outPath, "lig/simulation_comparison.pdf"), device="pdf", plot=scoreLigPlot))
suppressMessages(ggsave(paste0(outPath, "nolig/simulation_comparison.pdf"), device="pdf", plot=scoreNoLigPlot))

for (i in 1:rounds) {
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