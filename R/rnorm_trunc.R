#' placeholder
#' 
#' Helper function to calculate threshold for truncated distribution 
#' 
#' @param n placeholder
#' @param range placeholder
#' @param mu placeholder
#' @param sigma placeholder
#' @return The threshold value used to get the truncated distribution
#' @export
#' 

rnorm_trunc <- function(n, range, mu, sigma) {
  
  lower <- pnorm(min(range), mu, sigma)
  upper <- pnorm(max(range), mu, sigma)
  
  u <- runif(n, lower, upper)
  
  return(qnorm(u, mu, sigma))
}