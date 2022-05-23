#' Simulation of a disease 
#' 
#' This function is used to simulate a disease parameters. 
#' 
#' @param n_snp Amount of SNPs
#' @param prevalence Prevalence of disease in simulated population
#' @param h2 Herability parameter 
#' @param causal List of predetermined causal SNPs. Leave empty to get random causal SNPs. 
#' @param causal_n Amount of causal SNPs. Default value is 10 procent of total SNPs 
#' @param maf List of predetermined Minor Ale Frequency 
#' @param maf_low placeholder
#' @param maf_high placeholder
#' @param seed placeholder
#' @return placeholder
#'
#' 


sim_disease <- function(n_snp, prevalence, h2, causal = NULL, causal_n = round(n_snp / 10), maf = NULL, maf_low = 0.01, maf_high = 0.49, seed = NULL){
  if (!is.null(seed)) {set.seed(seed)}
  
  #Calculates a vector with causal SNPs at random positions if a causal vector not given
  if (is.null(causal)) {
    causal <- numeric(n_snp)
    causal[sample(1:n_snp, max(1, causal_n))] <- 1
  }
  
  #Calculates MAF if not given
  if (is.null(maf))  {
    maf <- runif(n_snp, maf_low, maf_high)
    
  }
  #Calcualtes betas for causal SNPs
  c <- sum(causal)
  beta <- ifelse(causal == 1, rnorm(n_snp, 0, sqrt(h2 / c)), 0)
  
  return(list(MAF = maf, BETA = beta, CAUSAL = causal, H2 = h2, PREVALENCE = prevalence, N_SNP = n_snp))
}
