#' Description 
#' 
#' GWAX
#' 
#' @param child placeholder
#' @return placeholder
#' @export 
#' 

GWAX <- function(child, include = rows_along(child$genotypes)) {
  
  p1_Status <- child$FAM$p1_Status
  p2_Status <- child$FAM$p2_Status
  child_status <- child$FAM$Status
  FBM <- child$genotypes
  
  #Creates a vector of the proxy statuses for the child
  x <- (child_status == 1 | p1_Status == 1 | p2_Status == 1) + 0
  
  GWAS(child, x, include = include)
}

