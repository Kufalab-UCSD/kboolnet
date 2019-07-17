# ###################################################
# compMatrix.R
# Adrian C
# 
# Takes two binary matrices, aligns them as closely
# as possible column-wise, and then computes the average
# simple matching coefficient across all columns.
#
# Note: Matrices must have the same number of rows. If
#   the number of columns does not match, extra columns
#   will not be included in score calculation.
#
# Args:
#   - mat1: First matrix
#   - mat2: Second matrix
# 
# Returns: numerical similarty value
# 
# Depends on: BoolNet, dplyr
# ###################################################

compMatrix <- function(mat1, mat2) {
  mat1 <- as.matrix(mat1)
  mat2 <- as.matrix(mat2)
  if (nrow(mat1) != nrow(mat2)) {
    stop("Number of matrix rows must match")
  }
  
  # If there is a mismatch in number of columns, add empty columns to the shorter matrix
  if (ncol(mat1) < ncol(mat2)) {
    mat1 <- cbind(mat1, matrix(nrow=nrow(mat1), ncol=ncol(mat2)-ncol(mat1)))
  } else if (ncol(mat2) < ncol(mat1)) {
    mat2 <- cbind(mat2, matrix(nrow=nrow(mat2), ncol=ncol(mat1)-ncol(mat2)))
  }
  
  # Align and score the matrices
  scores <- numeric()
  order  <- 1:ncol(mat1)
  for (i in 1:ncol(mat1)) {
    # Apply new ordering to matrix 2 and calculate the similarity (matching/total)
    vec1 <- as.vector(mat1)
    vec2 <- as.vector(mat2[,order])
    scores[i] <- sum(vec1 == vec2, na.rm=T)/length(vec1)
    
    # If perfect score, we can stop there
    if (scores[i] == 1) return(scores[i])
    
    # Create ordering for next round
    order <- c((ncol(mat1) - i + 1):ncol(mat1), 1:(ncol(mat1) - i))
  }
  
  # Return the highest score
  return(max(scores))
}
