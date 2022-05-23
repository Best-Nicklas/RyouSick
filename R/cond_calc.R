#' placeholder
#' 
#' This function calculates the conditional mu and sigma given a co-variance 
#' matrix. 
#' @param i Index of person 
#' @param covmatrix A co-variance matrix 
#' @return The conditional values of mu and sigma
#' @export
cond_calc <- function(i, covmatrix) {
  s11 <- covmatrix[i, i]
  s12 <- covmatrix[i, -i]
  s21 <- covmatrix[-i, i]
  s22 <- covmatrix[-i, -i]
  
  new_mu <- s12 %*% solve(s22) 
  new_sigma <- s11 - (s12 %*% solve(s22) %*% s21)
  return(list(mu = new_mu, sigma = new_sigma))
}
