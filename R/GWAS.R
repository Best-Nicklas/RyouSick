#' GWAS
#' Calculates a linear regression
#' @param person A .rds file with an FBM.code256 and accompanying FAM and MAP.
#' @param y Vector of regressands to regress on 
#' @param include Vector of rows to use in regression. Used with cross-validation. Default uses all rows.
#' @return A data.frame with slopes of each regression, standard errors of each slope, t-scores associated with each slope and P-values of each slope.
#' @export 
#' 
#' 
GWAS <- function(person, y, include = rows_along(person$genotypes)) {
  FBM <- person$genotypes
  
  #Uses function from bigSNPr package to do regression on FBM
  regr <- big_univLinReg(FBM, y[include], ind.train = include)
  #Adds column with P-values
  regr$p.value <- predict(regr, log10 = FALSE)
  
  return(data.frame(regr))
}
