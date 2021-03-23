#########################################
# calculateSensitivity.R
# Adrian C
#
# Takes two networks, generates a random
# set of initial states for them, then
# generates score based on how different
# the resulting attractors are.
#
# Requires: getPathAndAttractor,
# compMatrix
#
# Value: Sensitvity score
#########################################

findDifferenceForInitState <- function(state, network1, network2) {
  return(attractorDistance(getPathAndAttractor(network1, state)$attractor, getPathAndAttractor(network2, state)$attractor))
}

calculateSensitivity <- function(network1, network2, rounds) {
  if (any(network1$genes != network2$genes)) {
    stop("network1 and network2 must have the same set of genes!")
  }
  
  # Generate random initial states
  numGenes <- length(network1$genes)
  initStates <- matrix(sample(c(0,1), replace=TRUE, size=rounds*numGenes), nrow = rounds, ncol = numGenes)
  
  # Calculate the score
  score <- mean(apply(initStates, 1, findDifferenceForInitState, network1 = network1, network2 = network2))
  return(score)
}
