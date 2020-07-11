# ###################################################
# attractorDistance.R
# Adrian C
# 
# New version of compMatrix that uses the attractor
# similarity scoring algorithm described in RMut
# (https://doi.org/10.1371/journal.pone.0213736)
#
# Note: Matrices must have the same number of rows.
#
# Args:
#   - mat1: First matrix
#   - mat2: Second matrix
# 
# Returns: numerical similarity value
# 
# Depends on: BoolNet, numbers
# ###################################################

hamming <- function(vec1, vec2) {
  sum(as.numeric(vec1 != vec2))
}

attractorDistance <- function(attr1, attr2) {
  lengthLCM <- numbers::LCM(ncol(attr1), ncol(attr2))
  lengthGCD <- numbers::GCD(ncol(attr1), ncol(attr2))
  netSize <- nrow(attr1)
  
  scores <- sapply(0:(lengthGCD-1), function(shift) {
    return((1/(lengthLCM*netSize)) * sum(sapply(0:(lengthLCM-1), function(i) {
      return(hamming(attr1[,((shift + i) %% ncol(attr1)) + 1], attr2[,(i %% ncol(attr2)) + 1]))
    })))
  })
  
  return(min(scores))
}
