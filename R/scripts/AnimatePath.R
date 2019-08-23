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

############### Argument parsing and setup ###################
# Get commandline args
option_list = list(
  make_option(c("--file", "-f"), action="store", default=NA, type="character",
              help="Name of csv file containing path/attractor to be animated"),
  make_option(c("--out", "-o"), action="store", default="./animation.mp4", type="character",
              help="Name of mp4 file to which animation should be written. [default: animation.mp4]"),
  make_option(c("--textsize", "-t"), action="store", default=50, type="numeric",
              help="Size of frame number label. [default: 50]"),
  make_option(c("--zoom", "z"), action="store", default=200, type="numeric",
              help="Zoom factor for Cytoscape render. [default: 200].")
)
opt <- parse_args(OptionParser(option_list=option_list))

# Stop if no file provided
if (is.na(opt$file)) {
  stop("Please provide the file with data values to be plotted")
}

# Check if out is an mp4
if (!grepl("\\.mp4$", opt$out)) {
  opt$out <- paste0(opt$out, ".mp4")
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
nodeTypes <- getTableColumns(columns=c("name","type")) # Get types of nodes
rownames(nodeTypes) <- nodeTypes$name
components <- nodeTypes$name[nodeTypes$type == "component"] # Get list of only components
for (component in components) {
  if (component %in% rownames(path)) { # If the component is already in the path data, skip it
    next()
  }
  
  componentRows <- path[grepl(paste0(component, "_.*--0"), rownames(path)),,drop=F] # Get all the unbound state rows for that component
  path[component,] <- as.numeric(colSums(componentRows) > 0) # Create a new row for the component that = 1 if any unbound state = 1
}

################ Frame generation #####################
cat("Generating frames... DO NOT TOUCH YOUR COMPUTER WHILE DOING THIS OR THE GIF WILL NOT GENERATE PROPERLY.", "\n")

invisible(loadTableData(data.frame(row.names = nodeTypes$name, val = numeric(nrow(nodeTypes))))) # Load 0 to "val" column to reset edge routing
formatStr <- paste0("%0", nchar(toString(ncol(path))), "d.png") # Set formatting string for filenames

# Loop over each timepoint to create the frames
pb <- txtProgressBar(min=1, max=ncol(path), style=3) # Progress bar
for (i in 1:ncol(path)) {
  # Get the timepoint data
  cytoscapeTable <- path[,i,drop=F]
  colnames(cytoscapeTable) <- c("val")
  
  # Load it into Cytoscape
  invisible(loadTableData(cytoscapeTable))
  
  # Write the frame to the temporary directory
  exportImage(filename = paste0(tmpdir, sprintf(formatStr, i-1)), zoom = opt$zoom)
  
  setTxtProgressBar(pb, i) # Update progress bar
}
close(pb)

############### gif generation ########################
# Add frame number labels
cat("Labeling frames...", "\n")
frames <- dir(tmpdir)
pb <- txtProgressBar(min=1, max=length(frames), style=3) # Progress bar
for (i in 1:length(frames)) {
  frameFile <- paste0(tmpdir, frames[i])
  system2(command = "convert", args = c(frameFile, "-trim +repage -fill black -pointsize 50 -gravity NorthWest -annotate +10+10 '%t'", frameFile))
  setTxtProgressBar(pb, i)
}
close(pb)
 
# Use ffmpeg to stitch to a video
cat("Converting to video...", "\n")
system2(command = "ffmpeg", args = c("-framerate 5 -r 5 -loglevel panic -i", paste0(tmpdir, formatStr), opt$out))

cat("Animation written to file", opt$out, "\n")
