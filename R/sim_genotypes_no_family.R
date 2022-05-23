#' Description - Simulation with no family
#' @param n placeholder
#' @param disease placeholder
#' @param path placeholder
#' @param overwrite placeholder
#' @param n_blocks placeholder
#' @return placeholder
#' 

sim_genotypes_no_family <- function(n, disease, path, overwrite = T, n_blocks = 20) {
  
  # Load disease information
  cols <- disease$N_SNP
  MAF <- disease$MAF
  beta <- disease$BETA
  causal <- disease$CAUSAL
  prevalence <- disease$PREVALENCE
  h2 <- disease$H2
  
  # Calculate normalization constants 
  norm_const <- calc_normalization_consts(MAF, causal)
  mu <- norm_const$mu
  sigma <- norm_const$sigma
  
  #Checks if a FBM with the given name and dimensions exists, else creates one
  FBM <- verifyRds(path, overwrite, n, cols)
  
  # Prepare correct block indexes
  blocks <-  round(seq(0, n, length = n_blocks + 1))
  
  # Inserts values into FBM and calculate genetic liabilities in block sizes
  g_liabs <- future_lapply(1:(length(blocks) - 1), function(i) {  
    b_start <- blocks[i] + 1
    b_end <- blocks[i + 1]
    b_size <- (b_end - b_start + 1)
    
    FBM_temp <- matrix(rbinom(cols * b_size, 2, MAF),
                       nrow = b_size,
                       byrow = T)
    
    FBM$genotypes[b_start:b_end, ] <- FBM_temp
    
    g_block <- calc_gliab(FBM_temp, beta, mu, sigma)
    
  }, future.seed = T) %>% do.call("rbind", .) %>% as.numeric()
  
  # Saves liability and status information as well as SNP information in Rds
  threshold <- qnorm(prevalence, lower.tail = F)
  FBM$FAM$Genetic_Liability <- g_liabs
  FBM$FAM$Full_Liability <- g_liabs + rnorm(n, 0, sqrt(1 - h2))
  FBM$FAM$Status <- (FBM$FAM$Full_Liability > threshold) + 0
  
  FBM$MAP$MAF <- MAF
  FBM$MAP$BETA  <- beta
  
  return(FBM)
}
