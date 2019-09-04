#' theme_ROP
#' This function is a ggplot2 theme function
#'
#' @details This function is the ROP theme for ggplot2, it utilization is similar to any other ggplot2 theme function
#'
#' @author Geremy Clair
#' @export
#'
theme_ROP <- function(){
  theme_bw() %+replace%
    theme(panel.background  = element_blank(),
          plot.background = element_blank(),
          legend.background = element_rect(fill="transparent", colour=NA),
          legend.key = element_rect(fill="transparent", colour=NA),
          axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 14, face = "bold",hjust = 0.5))}
