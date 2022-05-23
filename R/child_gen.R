#' Description 
#' 
#' wdwdwd
#' 
#' @param p1 placeholder
#' @param p2 placeholder
#' @return wd34
#' 

child_gen = function(p1, p2){
  #Finds the positions at which both parents have 1
  ph1 = which(p1 == 1, arr.ind = T) %>% as_tibble()
  ph2 = which(p2 == 1, arr.ind = T) %>% as_tibble()
  ind_11 = bind_rows(inner_join(ph1, ph2, by = c("row", "col")), 
                     inner_join(ph1, ph2, by = c("row", "col"))) %>% distinct()
  
  #calculates the avg genotypes for the parents 
  temp <- (p1 + p2) / 2
  
  # samples from 0,1 on positions where avg 0.5
  temp[temp == 0.5] <- sample(0:1,length(temp[temp == 0.5]), replace = TRUE)
  
  # samples from 1,2 on positions where avg 1.5
  temp[temp == 1.5] <- sample(1:2,length(temp[temp == 1.5]), replace = TRUE)
  
  # samples from 0,1,2 on positions where both parents have 1
  temp[as.matrix(ind_11)] <- sample(0:2, size = nrow(ind_11), replace = T)  
  return(temp)
}