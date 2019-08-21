##############################################
# inhibitedNetwork.R
# Adrian C
#
# Create an inhibited BoolNet network based on
# a list of rxncon components or states to be
# inhibited.
#
# Dependencies: BoolNet
##############################################

inhibitedNetwork <- function(network, inhibList, names, symbols) {
  inhibNodes <- character()
  for (i in 1:length(inhib)) {
    # Make sure the inhibited nodes/components exist
    if (!(any(grepl(paste0("^", inhib[i], "(_.*--0|_.*-\\{0\\}|)$"), names)))) {
      stop("No neutral state found for inhibited node/component ", inhib[i], ". Please verify that ", inhib[i], " is a valid node/component in the rxncon system.")
    } 
    
    # Find which nodes the inhibition corresponds to
    inhibRegex <- paste0("(", c(paste0(inhib[i], "_.*--.*"), paste0(".*--", inhib[i], "_.*"),
                              paste0("^", inhib[i], "$"), paste0(inhib[i], "_.*-\\{.*\\}")), ")", collapse="|")
    inhibMatch <- grepl(inhibRegex, names)
    inhibNodes <- c(inhibNodes, symbols[inhibMatch])
    
    cat("Inhibitor", inhib[i], "matched to node(s)", paste0(names[inhibMatch], collapse=", "), "\n")
  }
  
  # Fix inhib nodes off
  network <- fixGenes(network, inhibNodes, 0)
  
  return(network)
}