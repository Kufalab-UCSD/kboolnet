####################################################
# removeLigand.R
# Adrian C
#
# Takes initial states of a rxncon BoolNet network and removes all
# occurences of a ligand from it. Ensures any  components the ligand 
# was complexed to are still present in unbound form.
#
# Args:
#   - ligand: rxncon component name to be removed
#   - names: character vector of names of rxncon nodes
#   - states: int vector of states
#
# Returns new int vector of states
#
# Depends on: BoolNet, dplyr
####################################################

unbindLigand <- function(ligand, names, states) {
  # If ligand isnt a component, throw an error
  if (grepl("_", ligand)) {
    stop("Ligand ", ligand, " is not a component.")
  }
  
  # If ligand doesn't exist in names, throw an error
  if (!any(grepl(paste0(ligand, "_"), names))) {
    stop("Ligand ", ligand, " does not exist in network.")
  }

  # Get all rxncon nodes that contain the ligand and are active
  ligandNodes   <- names[grepl(paste0(ligand, "_"), names) & states == 1]
  
  # Handling of complexes ligand has formed - we have to make sure unbound forms of whatever ligand was bound to are still present
  complexNames <- ligandNodes[grepl("--", ligandNodes) & !grepl("--0$", ligandNodes)] # Get all complexes the ligand is in
  if (length(complexNames) == 0) return(states) # If no complexes, quit
  ligandComplex <- data.frame(name = character(length(complexNames)),
                              compA = character(length(complexNames)),
                              compADomain = character(length(complexNames)),
                              compB = character(length(complexNames)),
                              compBDomain = character(length(complexNames)),
                              bindRxn = character(length(complexNames)),
                              unbindRxn = character(length(complexNames)),
                              unbound = character(length(complexNames)))
  ligandComplex$name <- complexNames
  
  # Extract complex info
  ligandComplex$compA <- lapply(ligandComplex$name, function(x) { # Get name of first component
    y <- strsplit(x, "--")[[1]][1] # Get component + domain name
    y <- gsub("_.*", "", y) # Delete everything after component name
    return(y)
  })
  ligandComplex$compADomain <- lapply(ligandComplex$name, function(x) { # Get domain of first component
    y <- strsplit(x, "--")[[1]][1] # Get component + domain name
    y <- gsub(".*\\[", "", y) # Remove everything before and after domain name
    y <- gsub("\\].*", "", y)
    return(y)
  })
  ligandComplex$compB <- lapply(ligandComplex$name, function(x) { # Get name of second component
    y <- strsplit(x, "--")[[1]][2] # Get component + domain name
    y <- gsub("_.*", "", y) # Delete everything after component name
    return(y)
  })
  ligandComplex$compBDomain <- lapply(ligandComplex$name, function(x) { # Get domain of second component
    y <- strsplit(x, "--")[[1]][2] # Get component + domain name
    y <- gsub(".*\\[", "", y) # Remove everything before and after domain name
    y <- gsub("\\].*", "", y)
    return(y)
  })
  
  # # Get reactions that form/dissociate complexes
  # for(i in 1:length(complexNames)) {
  #   # Get complex-forming reaction, tyring both component orderings
  #   bindRxn <- names[grepl(paste0(ligandComplex$compA, ".*_.*\\+_", ligandComplex$compB), names)]
  #   if (length(bindRxn) == 0) {
  #     bindRxn <- names[grepl(paste0(ligandComplex$compB, ".*_.*\\+_", ligandComplex$compA), names)]
  #   }
  #   ligandComplex$bindRxn[i] <- bindRxn
  #   
  #   # Same for dissociation reaction
  #   unbindRxn <- names[grepl(paste0(ligandComplex$compA, ".*_.*-_", ligandComplex$compB), names)]
  #   if (length(unbindRxn) == 0) {
  #     unbindRxn <- names[grepl(paste0(ligandComplex$compB, ".*_.*-_", ligandComplex$compA), names)]
  #   }
  #   ligandComplex$unbindRxn[i] <- unbindRxn
  # }
  
  # Reconstruct unbound forms of non-ligand components
  ligandComplex[ligandComplex$compA == ligand] <- ligandComplex[ligandComplex$compA == ligand] %>%
                                                  mutate(unbound = paste0(ligandComplex$compB, "_[", ligandComplex$compBDomain, "]--0"))
  ligandComplex[ligandComplex$compB == ligand] <- ligandComplex[ligandComplex$compB == ligand] %>%
                                                  mutate(unbound = paste0(ligandComplex$compA, "_[", ligandComplex$compADomain, "]--0"))
  
  # Remove any reactions that could result in complex formation/dissociation
  # states[names %in% ligandComplex$bindRxn] <- 0
  # states[names %in% ligandComplex$unbindRxn] <- 0

  # Remove all ligand-containing complexes or reactions
  states[names %in% complexNames] <- 0
  
  # Make sure unbound forms of non-ligand components are present
  states[names %in% ligandComplex$unbound] <- 1
  
  return(states)
}