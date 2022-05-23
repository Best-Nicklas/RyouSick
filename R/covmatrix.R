#' Co-variance matrix
#' 
#' This function creates the theoretical co-variance matrix from the 
#' multivariate normal distribution which is used to model the liabilities.
#' 
#' @param h2 Heritability parameter.
#' @param sib Number of siblings.
#' @return A co-variance matrix, created from the value of h2 and 
#' the amount of sib. 
#' 
#' 
#' @export
#' 
covmatrix <- function(h2, n_sib = 0) {
  cov <- matrix(h2/2, 4 + n_sib, 4 + n_sib)
  diag(cov) <- 1
  cov[3,4] <- cov[4,3] <- 0
  cov[1:2, 1] <- cov[1, 1:2] <- h2
  
  return(cov)
}