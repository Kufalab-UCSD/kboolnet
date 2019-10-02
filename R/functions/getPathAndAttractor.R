########################################
# getPathAndAttractor
# Adrian C
#
# Small function that returns both the
# trajectory and attractor for a network
# given the network and an initial state
# vector.
#
# Dependencies: BoolNet
#
# Value:
# A list with two entries: path and attractor
########################################

getPathAndAttractor <- function(network, states, names=NULL) {
  path <- t(getPathToAttractor(network, states)) # Simulate the path
  rownames(path) <- names
  states <- path[,ncol(path)] # Use last path state as new start for attractor
  attractor <- t(getPathToAttractor(network, states)) # Get the attractor
  attractor <- attractor[,1:(ncol(attractor)-1), drop=FALSE] # Remove last repeated column from attractor
  rownames(attractor) <- names
  
  return(list("path" = path, "attractor" = attractor))
}
