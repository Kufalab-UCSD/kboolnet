##############################################
# fixedNetwork.R
# Adrian C
#
# Create a fixed BoolNet network based on
# a list of rxncon components or states to be
# inhibited.
#
# Dependencies: BoolNet
##############################################

fixedNetwork <- function(network, inhib, names, symbols, effect = 0, regex = FALSE) {
  if (!effect %in% c(0, 1)) {
    stop("Invalid effect on the network. Must be either 0 to fix off or 1 to fix on.")
  }
  inhibNodes <- character()
  for (i in 1:length(inhib)) {
    # Try to match to a node
    if (!regex) {
      # Make sure the inhibited nodes/components exist
      if (!(any(grepl(paste0("^", inhib[i], "(_.*--0|_.*-\\{0\\}|)$"), names)))) {
        stop("No neutral state found for inhibited node/component ", inhib[i], ". Please verify that ", inhib[i], " is a valid node/component in the rxncon system.")
      }

      # Find which nodes the inhibition corresponds to
      inhibRegex <- paste0("(", c(paste0(inhib[i], "_.*--.*"), paste0(".*--", inhib[i], "_.*"),
                                paste0("^", inhib[i], "$"), paste0(inhib[i], "_.*-\\{.*\\}")), ")", collapse="|")
    } else {
      inhibRegex <- escapeRegex(inhib[i])
    }

    inhibMatch <- grepl(inhibRegex, names)
    inhibNodes <- c(inhibNodes, symbols[inhibMatch])

    if (effect == 0) {
      cat("Inhibitor", inhib[i], "matched to node(s)", paste0(names[inhibMatch], collapse=", "), "\n")
    } else if (effect == 1) {
      cat("Activator", inhib[i], "matched to node(s)", paste0(names[inhibMatch], collapse=", "), "\n")
    }
  }

  # Fix inhib nodes off
  network <- BoolNet::fixGenes(network, inhibNodes, effect)

  return(network)
}
