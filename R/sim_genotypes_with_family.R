#' Description - 
#' 
#' Simulation with family (Fixed and non fixed)
#' @param n placeholder
#' @param disease placeholder
#' @param path placeholder
#' @param n_sibs placeholder
#' @param overwrite placeholder
#' @param n_blocks placeholder
#' @return  placeholder
#' 

sim_genotypes_with_family <- function(n, disease, path, n_sibs = NULL, overwrite = T, n_blocks = 20) {
  # Load disease information
  cols <- disease$N_SNP
  MAF <- disease$MAF
  beta <- disease$BETA
  causal <- disease$CAUSAL
  prevalence <- disease$PREVALENCE
  h2 <- disease$H2
  
  # Create or find FBM file to fill with child genotypes
  FBM <- verifyRds(path, overwrite, n, cols)
  
  # Create vector with number of sibs for each child
  if (!is.null(n_sibs)) {
    sibs_pr_child <- if (length(n_sibs) == 1) 
    {rep(n_sibs,n)} else 
      sample(min(n_sibs):max(n_sibs), n, replace = T)
  }
  
  # Calculate normalization constants 
  norm_const <- calc_normalization_consts(MAF, causal)
  mu <- norm_const$mu
  sigma <- norm_const$sigma
  
  #prepare for block calculations of genetic liabilities 
  blocks <-  round(seq(0, n, length = n_blocks + 1))
  g_liabs <- future_lapply(1:(length(blocks) - 1), function(i) {  
    b_start <- blocks[i] + 1
    b_end <- blocks[i + 1]
    b_size <- (b_end - b_start + 1)
    
    # simulate parent 1 genotypes and genetic liabilities 
    p1 <- matrix(rbinom(cols * b_size, 2, MAF),
                 nrow = b_size,
                 byrow = T)
    p1_gliab <- calc_gliab(p1, beta, mu, sigma)
    
    # simulate parent 2 genotypes and genetic liabilities 
    p2 <- matrix(rbinom(cols * b_size, 2, MAF),
                 nrow = b_size,
                 byrow = T)
    p2_gliab <- calc_gliab(p2, beta, mu, sigma)
    
    # Simulate child genotypes and genetic liabilities - store genotypes in FBM
    child <- child_gen(p1, p2)
    FBM$genotypes[b_start:b_end, ] <- child
    child_gliab <- calc_gliab(child, beta, mu, sigma)
    
    # Generate sibs for each child and calculate their genetic liabilities
    sibs_gliab <- vector(mode = "list", b_size)
    if (!is.null(n_sibs)) {
      for (s in 1:b_size) {
        sibs <- child_gen(matrix(rep(p1[s, ], sibs_pr_child[s]), sibs_pr_child[s], cols, byrow = T),
                          matrix(rep(p2[s, ], sibs_pr_child[s]), sibs_pr_child[s], cols, byrow = T))
        
        sibs_gliab[[s]] <- calc_gliab(sibs, beta, mu, sigma)
      }
    }
    
    tibble(child_gliab, p1_gliab, p2_gliab, sibs_gliab)
    
  }, future.seed = T) %>% do.call("bind_rows", .)
  
  # Calculate full liabilities/status and insert in rds file object
  threshold <- qnorm(prevalence, lower.tail = F)
  FBM$FAM$Genetic_Liability <- g_liabs$child_gliab
  FBM$FAM$Full_Liability <- g_liabs$child_gliab + rnorm(n, 0, sqrt(1 - h2))
  FBM$FAM$Status <- (FBM$FAM$Full_Liability > threshold) + 0
  FBM$FAM$p1_Status <- (g_liabs$p1_gliab + rnorm(n, 0, sqrt(1 - h2)) > threshold) + 0
  FBM$FAM$p2_Status <- (g_liabs$p2_gliab + rnorm(n, 0, sqrt(1 - h2)) > threshold) + 0
  
  FBM$MAP$MAF <- MAF
  FBM$MAP$BETA  <- beta
  
  if (!is.null(n_sibs)) {
    sibs_Full_Liability <- purrr::map(g_liabs$sibs_gliab, .f = ~ 
                                        {if ( is.null(.x)) NULL 
                                          else .x + rnorm(length(.x), 0, sqrt(1 - h2))})
    FBM$FAM$sibs_Status <- purrr::map(sibs_Full_Liability, .f = ~ 
                                        {if ( is.null(.x)) NULL
                                          else(.x > threshold) + 0})
  }
  
  return(FBM)
  
  
}