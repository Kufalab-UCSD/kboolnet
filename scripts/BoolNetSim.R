#!/usr/bin/env Rscript

options(stringsAsFactors = F)
invisible(library(BoolNet))
invisible(library(ggplot2))
invisible(library(tidyr))
invisible(library(optparse))

#############################################################
# Disclaimer
# 
# This is a basic script that simulates a bipartite Boolean 
# network created from a rxncon model. It requires three files:
#   model.boolnet
#   model_initial_vals.csv
#   model_symbols.csv
# 
# The script runs two simulations and plots the trajectories.
# The script creates five new files:
#   model.pdf
#   model_trajectory_first.csv
#   model_attractor.pdf
#   model_trajectory_second.csv
#   model_new_attractor.csv
# These five files will be overwritten each time you run the 
# script.
# No guarantee.
#############################################################


#############################################################
############# 1 NETWORK LOADING AND SIMULATION ##############

# accept network as commandline argument
option_list = list(
  make_option("--model", action="store", default=NA, type='character',
              help="Prefix of model to be simulated")
)
opt = parse_args(OptionParser(option_list=option_list))

# Provide network name
if(!is.na(opt$model)) {
  filePrefix <- opt$model
} else {
  stop("Please provide the name of the model.")
}


# Load Network as symbolic
network <- loadNetwork(paste0(filePrefix, ".boolnet"), symbolic=TRUE)

# Load symbols (data frame) and remove spaces
symbolMapping     <- read.csv(paste0(filePrefix, "_symbols.csv"), header=F)
symbolMapping$V2  <- gsub("[[:space:]]", "", symbolMapping$V2)
symbolMapping     <- symbolMapping[,1:2]

# Define initial states 
initStates        <- read.csv(paste0(filePrefix, '_initial_vals.csv'), header=F)
initStates        <- initStates[,1:2]

# Merge Symbols and initial state
init_symbols              <- cbind(symbolMapping, initStates)
colnames(init_symbols)[4] <- "initStates"

# First simulation
initState <- init_symbols[,4]
path      <- getPathToAttractor(network,initState)
t_path    <- t(path)
write.table(t_path,file=(paste0(filePrefix, '_trajectory_first_simulation.csv')), sep=",", col.names=F)

# Write data in format readable by Cytoscape
colnames(t_path)    <- sprintf("t%03d", 0:(ncol(t_path)-1))
rownames(t_path) <- symbolMapping[,2]
write.table(t_path,file=(paste0(filePrefix, '_trajectory_cytoscape.csv')), sep=",")

### Second simulation ###
new_initState <- path[nrow(path),]
new_path      <- getPathToAttractor(network,new_initState)
t_new_path    <- t(new_path)
write.table(t_new_path,file=(paste0(filePrefix, '_trajectory_second_simulation.csv')), sep=",", col.names=F)

new_attractor <- t_new_path[,ncol(t_new_path)]
write.table(new_attractor,file=(paste0(filePrefix, '_new_attractor.csv')), sep=",", col.names=F)


#############################################################
############# 2 DATA REARRANGEMENT AND PLOTTING #############

############# Sort both simulations #########################

# Assignments for first simulation
path_length       <- nrow(path)
path_mat          <- as.matrix(path)
symbolVec         <- as.character(symbolMapping[,2])
names(symbolVec)  <- as.character(symbolMapping[,1])
symbolVec         <- symbolVec[colnames(path_mat)]

colnames(path_mat)  <- symbolVec
path_mat            <- t(path_mat)

# Assignments for second simulation
new_path_length   <- nrow(new_path)
new_path_mat      <- as.matrix(new_path)
symbolVec         <- as.character(symbolMapping[,2])
names(symbolVec)  <- as.character(symbolMapping[,1])
symbolVec         <- symbolVec[colnames(new_path_mat)]

colnames(new_path_mat)  <- symbolVec
new_path_mat            <- t(new_path_mat)

both_path_mat <- cbind(path_mat, new_path_mat)

# Sort both in decreasing activation order
for(i in ncol(both_path_mat):1){
  both_path_mat  <- both_path_mat[order(both_path_mat[,i], decreasing = T), ]
}

path_mat     <- both_path_mat[,1:ncol(path_mat)]
new_path_mat <- both_path_mat[,-(1:ncol(path_mat))]

############# First simulation  #############################


# Set attractor values to 2 for different coloring later on
attractor_mat                                             <- path_mat[,-(1:(ncol(path_mat) - (ncol(new_path_mat)-1)))]
attractor_mat[attractor_mat == 1]                         <- 2
path_mat[,-(1:(ncol(path_mat) - (ncol(new_path_mat)-1)))] <- attractor_mat

path_mat2   <- data.frame(path_mat)

# Replace t..0 with t000
colnames(path_mat2) <- sprintf("t%03d", 0:(path_length-1))

# Add columnvector with symbols
path_mat2$symbols <- rownames(path_mat2)

# Rearranges data frame
path_mat2         <- path_mat2 %>% gather(t,value, 1:path_length)
path_mat2$symbols <- factor(path_mat2$symbols, levels = rownames(path_mat)[dim(path_mat)[1]:1])
path_mat2$t       <- factor(path_mat2$t)


# Plot first simulation 
p                 <- ggplot(path_mat2, aes(t, symbols)) + geom_tile(aes(fill = value), colour = "white") + scale_fill_gradient2(low = "white", mid = "steelblue", high = "red", midpoint = 1)
base_size         <- 8
p                 <- p + theme_grey(base_size = base_size) + labs(x = "", y = "") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + coord_fixed(ratio=0.5) + theme(legend.position = "none")
p
ggsave(paste0(filePrefix, ".pdf"), plot = last_plot(), height=(5 + length(levels(p$data$symbols)) * 0.25), width=(10 + length(levels(p$data$t)) * 0.45), scale = 1, units = "cm")

#############################################################

############# Second simulation   ###########################

new_path_mat2 <- data.frame(new_path_mat)

# Replace t..0 with t000
colnames(new_path_mat2) <- sprintf("t%03d", 0:(new_path_length-1))

# Add columnvector with symbols
new_path_mat2$symbols   <- rownames(new_path_mat2)

# Rearrange data frame
new_path_mat2         <- new_path_mat2 %>% gather(t,value, 1:new_path_length)
new_path_mat2$symbols <- factor(new_path_mat2$symbols, levels = rownames(new_path_mat)[dim(new_path_mat)[1]:1])
new_path_mat2$t       <- factor(new_path_mat2$t)

# Plot second simulation 
p2                    <- ggplot(new_path_mat2, aes(t, symbols)) + geom_tile(aes(fill = value), colour = "white") + scale_fill_gradient(low = "white", high = "steelblue")
base_size             <- 8
p2                    <- p2 + theme_grey(base_size = base_size) + labs(x = "", y = "") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + coord_fixed(ratio=0.5) + theme(legend.position = "none")
p2
ggsave(paste0(filePrefix, "_attractor.pdf"), plot = last_plot(), height=(5 + length(levels(p$data$symbols)) * 0.25), width=(10 + length(levels(p$data$t)) * 0.45), scale = 1, units = "cm")

#############################################################

