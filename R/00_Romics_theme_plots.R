#' theme_ROP()
#' @description This function is a ggplot2 theme function
#' @details This function is the ROP theme for ggplot2, it utilization is similar to any other ggplot2 theme function
#' @author Geremy Clair
#' @export
theme_ROP <- function(){
  ggplot2::theme_bw() +
    ggplot2::theme(panel.background  = ggplot2::element_blank(),
          plot.background = ggplot2::element_blank(),
          legend.background = ggplot2::element_rect(fill="transparent", colour=NA),
          legend.key = ggplot2::element_rect(fill="transparent", colour=NA),
          axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(size = 14, face = "bold",hjust = 0.5))}

#' theme_map()
#' @description This function is a ggplot2 theme function for the spatial maps that can be created using romics
#' @details This function is the ROP theme for ggplot2, it utilization is similar to any other ggplot2 theme function
#' @author Geremy Clair
#' @export

theme_map<- function(){
  font <- "Helvetica"   #assign font family up front

  ggplot2::theme_minimal() +    #replace elements we want to change

    ggplot2::theme(

      #grid elements
      panel.background = ggplot2::element_rect(fill = 'black', color = 'gray50'),
      panel.grid.major = ggplot2::element_line(color = 'white', linetype = 'dotted',linewidth = 0.2),
      panel.grid.minor = ggplot2::element_line(color = 'white', linetype = 'dotted',linewidth = 0.2),
      #panel.background
      #panel.grid.major = ggplot2::element_blank(),    #strip major gridlines
      #panel.grid.minor = ggplot2::element_blank(),    #strip minor gridlines
      #axis.ticks = ggplot2::element_blank(),          #strip axis ticks

      #since theme_minimal() already strips axis lines,
      #we don't need to do that again

      #text elements
      plot.title = ggplot2::element_text(             #title
        family = font,            #set font family
        size = 16,                #set font size
        face = 'bold',            #bold typeface
        hjust = 0,                #left align
        vjust = 2),               #raise slightly

      aspect.ratio=1,

      plot.subtitle = ggplot2::element_text(          #subtitle
        family = font,            #font family
        size = 14),               #font size

      plot.caption = ggplot2::element_text(           #caption
        family = font,            #font family
        size = 9,                 #font size
        hjust = 1),               #right align

      axis.title = ggplot2::element_text(             #axis titles
        family = font,            #font family
        face = 'bold',            #bold typeface
        size = 12),               #font size

      axis.text = ggplot2::element_text(              #axis text
        family = font,            #axis famuly
        size = 12),                #font size

      axis.text.x = ggplot2::element_text(            #margin for axis text
        margin=ggplot2::margin(5, b = 10))


    )
}
