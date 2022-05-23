#' Create .rds file.
#' 
#' This function is used create a file-backed matrices and save it as an .rds file with related information. 
#' 
#' @param nrow Amount of rows.
#' @param ncol Amount of columns.
#' @param path Path where to save and find the file. 
#' 
#' @return An rds object containing a list with a file-backed matrix, FAM and MAP information in a list.
#' 
#' @export
#'


createRds <-  function(path, nrow, ncol) {
  G = FBM.code256(nrow = nrow, # number of rows
                  ncol = ncol, # number of columns
                  code = c(0L, 1L, 2L, rep(NA_integer_, 256 - 3)),
                  backingfile = path) 
  
  obj.bigsnp = list(
    genotypes = G,
    MAP = tibble(SNP_ID = 1:ncol),
    FAM = tibble(ID = 1:nrow))
  snp_save(obj.bigsnp)
  
  
  return(obj.bigsnp)
}

