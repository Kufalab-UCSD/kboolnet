#!/usr/bin/env Rscript

##############################################################
# AnimatePath.R
# Adrian C
#
# Script to generate an animation from an attractor/path csv
# and a network open in a Cytoscape session.
#
# Requires imagemagick and ffmpeg commandline tools
##############################################################

############## Library loading ###############################
options(stringsAsFactors = F)
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(RCy3))
library(kboolnet)

############### Argument parsing and setup ###################
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
default <- list(file=NA, out="./animation.gif", textsize=50, zoom=200)
opt <- setDefaults(config, default)

# Stop if no file provided
if (is.na(opt$file)) {
  stop("Please provide the file with data values to be plotted")
}

# Check if out is an mp4
if (!grepl("\\.gif$", opt$out)) {
  opt$out <- paste0(opt$out, ".gif")
}

# Create a temp folder to put the frames in
tmpdir <- paste0(tempdir(), "/frames/")
if (dir.exists(tmpdir)) {
  unlink(tmpdir, recursive = TRUE)
}
dir.create(tmpdir)

################# Check if Cytoscape environment is prepared ##################
# Check if Cytoscape running
tryCatch({
  invisible(cytoscapePing())
}, error=function(e){
  print(e)
  stop("Couldn't connect to Cytoscape. Please ensure that Cytoscape is running.")
})

# Check if there is a network loaded
tryCatch({
  invisible(getNetworkName())
}, error=function(e){
  print(e)
  stop("No loaded networks detected. Please ensure you have loaded a network.")
})

# Set correct visual style
if(!("rxncon2animation" %in% getVisualStyleNames())) {
  stop("rxncon2animation style not found. Please ensure you have loaded it before running.")
} else {
  invisible(setVisualStyle("rxncon2animation"))
}

################# Path loading/processing ##################
# Load path and set first col as row names
path <- read.csv(opt$file, header=TRUE)
rownames(path) <- path[,1]
path <- path[,2:ncol(path), drop=F]

# Try and match components to path data
nodeTypes <- getTableColumns(columns=c("rxnconID","type")) # Get types of nodes
rownames(nodeTypes) <- nodeTypes$rxnconID
components <- nodeTypes$rxnconID[nodeTypes$type == "component"] # Get list of only components
for (component in components) {
  if (component %in% rownames(path)) { # If the component is already in the path data, skip it
    next()
  }

  componentRows <- path[grepl(paste0(component, "_.*--0"), rownames(path)),,drop=F] # Get all the unbound state rows for that component
  path[component,] <- as.numeric(colSums(componentRows) > 0) # Create a new row for the component that = 1 if any unbound state = 1
}

################ Frame generation #####################
cat("Generating frames... DO NOT TOUCH YOUR COMPUTER WHILE DOING THIS OR THE GIF WILL NOT GENERATE PROPERLY.", "\n")

invisible(loadTableData(data.frame(rxnconID = nodeTypes$rxnconID, val = numeric(nrow(nodeTypes))), data.key.column = "rxnconID", table.key.column = "rxnconID")) # Load 0 to "val" column to reset edge routing
formatStr <- paste0("%0", nchar(toString(ncol(path))), "d.png") # Set formatting string for filenames

# Loop over each timepoint to create the frames
pb <- txtProgressBar(min=0, max=ncol(path), style=3) # Progress bar
for (i in 1:ncol(path)) {
  # Get the timepoint data
  cytoscapeTable <- path[,i,drop=F]
  colnames(cytoscapeTable) <- c("val")

  # Load it into Cytoscape
  invisible(loadTableData(cytoscapeTable, table.key.column = "rxnconID"))

  # Write the frame to the temporary directory
  exportImage(filename = paste0(tmpdir, sprintf(formatStr, i-1)), zoom = opt$zoom)

  setTxtProgressBar(pb, i) # Update progress bar
}
close(pb)

############### gif generation ########################
# Add frame number labels
cat("Labeling frames...", "\n")
frames <- dir(tmpdir)
pb <- txtProgressBar(min=0, max=length(frames), style=3) # Progress bar
for (i in 1:length(frames)) {
  frameFile <- normalizePath(paste0(tmpdir, frames[i]), mustWork = FALSE)
  if(.Platform$OS.type == "unix") {
    system2(command = "convert", args = c(frameFile, "-trim +repage -fill black -pointsize 50 -gravity NorthWest -annotate +10+10 '%t'", addQuotes(frameFile)))
  } else {
    system2(command = "magick", args = c("convert", frameFile, "-trim +repage -fill black -pointsize 50 -gravity NorthWest -annotate +10+10 '%t'", addQuotes(frameFile)))
  }
  setTxtProgressBar(pb, i)
}
close(pb)

# Use imagemagick to stitch to a gif
cat("Converting to video...", "\n")
if(.Platform$OS.type == "unix") {
  system2(command = "convert", args = c("-limit memory 2GiB -delay 20 -loop 0", addQuotes(normalizePath(paste0(tmpdir, '/*'), mustWork = FALSE)), addQuotes(opt$out)))
} else {
  system2(command = "magick", args = c("convert -delay 20 -loop 0", addQuotes(normalizePath(paste0(tmpdir, '/*'), mustWork = FALSE)), addQuotes(opt$out)))
}

cat("Animation written to file", opt$out, "\n")
