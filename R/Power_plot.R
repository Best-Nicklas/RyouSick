#' Power plot
#' 
#' This function plots the power of 
#' 
#' @param x 
#' @return placeholder
#' @export 
Power_plot <- function(x) {
  powerplot <- x %>% mutate (causal_snp = p.value < 5e-3) %>%
    arrange(abs(estim)) %>%
    mutate (cpower = cumsum(causal_snp))
  
  return(powerplot %>% ggplot(aes(x = estim, y = cpower)) + geom_line())
}