#' distribBoxplot
#' This function will plot the boxplot of the sample distribution using ggplot2
#'
#' @param romics_object has to be a log transformed romics_object created using romicsCreateObject() and transformed using the function log2transform() or log10transform()
#'
#' @details create a ggplot2 graphic output displaying the normalized or not intensites (or logged intensities) within each sample. ggplot2 methods can be used to visually adjust the plot.
#'
#' @return a plot generated using ggplot2 is generated with this function.
#'
#' @author Geremy Clair
#' @export
#'
distribBoxplot <- function(romics_object){
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}

  data_graph<-melt(romics_object$data,id.vars=0)
  fill_graph<-romics_object$colors
  fill_graph<-fill_graph[!is.na(data_graph[,2])]
  data_graph<-data_graph[!is.na(data_graph[,2]),]
  data_graph<-cbind(data_graph,fill_graph)
  colors<-as.character(t(romics_object$metadata[romicsFactorNames(romics_object)=="colors_romics",]))
  factor_graph<-data.frame(romics_object$metadata[romicsFactorNames(romics_object)==romics_object$main_factor,])
  factor_graph<-data.frame(variable=colnames(factor_graph),factor=as.character(t(factor_graph[1,])))
  data_graph<-merge(data_graph,factor_graph,by="variable")


p<-ggplot(aes(x=variable, y=value),data=data_graph)+
  geom_boxplot(aes(fill=factor),alpha=0.8)+
  scale_fill_manual(values=levels(factor(data_graph$fill_graph, levels = unique(data_graph$fill_graph))))+
  ggtitle("Boxplot of the data distribution within each sample")+
  theme_ROP()+
  xlab("Sample")

if(romicsLogCheck(romics_object)==T)
{if(sum(grepl("log2transform\\(",romics_object$steps))>0){
  p<-p+ylab("Log2(intensities)")}else{p<-p+ylab("Log10(intensities)")}}else{p<-p+ylab("Intensities")}

return(p)
}

#' distribViolin
#' This function will plot the Violin plot of the sample distribution using ggplot2
#'
#' @param romics_object has to be a log transformed romics_object created using romicsCreateObject() and transformed using the function log2transform() or log10transform()
#'
#' @details create a ggplot2 graphic output displaying the normalized or not intensites (or logged intensities) within each sample. ggplot2 methods can be used to visually adjust the plot.
#'
#' @return a plot generated using ggplot2 is generated with this function.
#'
#' @author Geremy Clair
#' @export
#'
distribViolin <- function(romics_object){
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}

  data_graph<-melt(romics_object$data,id.vars=0)
  fill_graph<-romics_object$colors
  fill_graph<-fill_graph[!is.na(data_graph[,2])]
  data_graph<-data_graph[!is.na(data_graph[,2]),]
  data_graph<-cbind(data_graph,fill_graph)
  colors<-as.character(t(romics_object$metadata[romicsFactorNames(romics_object)=="colors_romics",]))
  factor_graph<-data.frame(romics_object$metadata[romicsFactorNames(romics_object)==romics_object$main_factor,])
  factor_graph<-data.frame(variable=colnames(factor_graph),factor=as.character(t(factor_graph[1,])))
  data_graph<-merge(data_graph,factor_graph,by="variable")


  p<-ggplot(aes(x=variable, y=value),data=data_graph)+
    geom_violin(aes(fill=factor),alpha=0.8)+
    scale_fill_manual(values=levels(factor(data_graph$fill_graph, levels = unique(data_graph$fill_graph))))+
    ggtitle("Boxplot of the data distribution within each sample")+
    theme_ROP()+
    xlab("Sample")

  if(romicsLogCheck(romics_object)==T)
  {if(sum(grepl("log2transform\\(",romics_object$steps))>0){
    p<-p+ylab("Log2(intensities)")}else{p<-p+ylab("Log10(intensities)")}}else{p<-p+ylab("Intensities")}

  return(p)
}

#' distribJitter
#' This function will plot the jitter plot of the sample distribution using ggplot2
#'
#' @param romics_object has to be a log transformed romics_object created using romicsCreateObject() and transformed using the function log2transform() or log10transform()
#'
#' @details create a ggplot2 graphic output displaying the normalized or not intensites (or logged intensities) within each sample. ggplot2 methods can be used to visually adjust the plot.
#'
#' @return a plot generated using ggplot2 is generated with this function.
#'
#' @author Geremy Clair
#' @export
#'
distribJitter <- function(romics_object){
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}

  data_graph<-melt(romics_object$data,id.vars=0)
  fill_graph<-romics_object$colors
  fill_graph<-fill_graph[!is.na(data_graph[,2])]
  data_graph<-data_graph[!is.na(data_graph[,2]),]
  data_graph<-cbind(data_graph,fill_graph)
  colors<-as.character(t(romics_object$metadata[romicsFactorNames(romics_object)=="colors_romics",]))
  factor_graph<-data.frame(romics_object$metadata[romicsFactorNames(romics_object)==romics_object$main_factor,])
  factor_graph<-data.frame(variable=colnames(factor_graph),factor=as.character(t(factor_graph[1,])))
  data_graph<-merge(data_graph,factor_graph,by="variable")


  p<-ggplot(aes(x=variable, y=value,fill=factor, colours=fill_graph),data=data_graph)+
    geom_jitter(aes(color=factor),alpha=0.8)+
    scale_color_manual(values=levels(factor(data_graph$fill_graph, levels = unique(data_graph$fill_graph))))+
    ggtitle("Jitter plot of the data distribution within each sample")+
    theme_ROP()+
    xlab("Sample")

  if(romicsLogCheck(romics_object)==T)
  {if(sum(grepl("log2transform\\(",romics_object$steps))>0){
    p<-p+ylab("Log2(intensities)")}else{p<-p+ylab("Log10(intensities)")}}else{p<-p+ylab("Intensities")}

  return(p)
}

#' distribSina
#' This function will plot the Sina plot of the sample distribution using ggplot2
#'
#' @param romics_object has to be a log transformed romics_object created using romicsCreateObject() and transformed using the function log2transform() or log10transform()
#'
#' @details create a ggplot2 graphic output displaying the normalized or not intensites (or logged intensities) within each sample. ggplot2 methods can be used to visually adjust the plot.
#'
#' @return a plot generated using ggplot2 is generated with this function.
#'
#' @author Geremy Clair
#' @export
#'
distribSina <- function(romics_object){
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}

  data_graph<-melt(romics_object$data,id.vars=0)
  fill_graph<-romics_object$colors
  fill_graph<-fill_graph[!is.na(data_graph[,2])]
  data_graph<-data_graph[!is.na(data_graph[,2]),]
  data_graph<-cbind(data_graph,fill_graph)
  colors<-as.character(t(romics_object$metadata[romicsFactorNames(romics_object)=="colors_romics",]))
  factor_graph<-data.frame(romics_object$metadata[romicsFactorNames(romics_object)==romics_object$main_factor,])
  factor_graph<-data.frame(variable=colnames(factor_graph),factor=as.character(t(factor_graph[1,])))
  data_graph<-merge(data_graph,factor_graph,by="variable")


  p<-ggplot(aes(x=variable, y=value,fill=factor, colours=fill_graph),data=data_graph)+
    geom_sina(aes(color=factor),alpha=0.8)+
    scale_color_manual(values=levels(factor(data_graph$fill_graph, levels = unique(data_graph$fill_graph))))+
    ggtitle("Sina plot of the data distribution within each sample")+
    theme_ROP()+
    xlab("Sample")

  if(romicsLogCheck(romics_object)==T)
  {if(sum(grepl("log2transform\\(",romics_object$steps))>0){
    p<-p+ylab("Log2(intensities)")}else{p<-p+ylab("Log10(intensities)")}}else{p<-p+ylab("Intensities")}

  return(p)

}

#' distribRidges
#' This function will plot the Density ridges plot of the sample distribution using ggplot2
#'
#' @param romics_object has to be a log transformed romics_object created using romicsCreateObject() and transformed using the function log2transform() or log10transform()
#'
#' @details create a ggplot2 graphic output displaying the normalized or not intensites (or logged intensities) within each sample. ggplot2 methods can be used to visually adjust the plot.
#'
#' @return a plot generated using ggplot2 is generated with this function.
#'
#' @author Geremy Clair
#' @export
#'
distribRidges <- function(romics_object){
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}

  data_graph<-melt(romics_object$data,id.vars=0)
  fill_graph<-romics_object$colors
  fill_graph<-fill_graph[!is.na(data_graph[,2])]
  data_graph<-data_graph[!is.na(data_graph[,2]),]
  data_graph<-cbind(data_graph,fill_graph)
  colors<-as.character(t(romics_object$metadata[romicsFactorNames(romics_object)=="colors_romics",]))
  factor_graph<-data.frame(romics_object$metadata[romicsFactorNames(romics_object)==romics_object$main_factor,])
  factor_graph<-data.frame(variable=colnames(factor_graph),factor=as.character(t(factor_graph[1,])))
  data_graph<-merge(data_graph,factor_graph,by="variable")


  p<-ggplot(aes(x=value, y=variable),data=data_graph)+
    geom_density_ridges(scale=1.5,aes(fill=factor),alpha=0.8,panel_scaling = T)+
    scale_fill_manual(values=levels(factor(data_graph$fill_graph, levels = unique(data_graph$fill_graph))))+
    ggtitle("Density ridges plot of the data distribution within each sample")+
    theme_ROP()+
    ylab("Sample")

  if(romicsLogCheck(romics_object)==T)
  {if(sum(grepl("log2transform\\(",romics_object$steps))>0){
    p<-p+xlab("Log2(intensities)")}else{p<-p+xlab("Log10(intensities)")}}else{p<-p+xlab("Intensities")}

  return(p)
}

#' distribHistogram
#' This function will display a grid.arranged histogram plot showing the distribution of the data in each sample
#'
#' @param romics_object has to be a log transformed romics_object created using romicsCreateObject() and transformed using the function log2transform() or log10transform()
#' @param bin has to be a numerical value indicating the width of the frequency bins to use for the visualization
#' @param scale_y has to be a numerical/double vector of length 2 indicating the y_limits of each graph. if too low those will be automatically adjusted. by default  scale_y=c(0,100) will be used
#' @param scale_x has to be a numerical/double vector of length 2 indicating the x_limits of each graph. if too low those will be automatically adjusted. by default  scale_x=c(0,100) will be used
#' @param col has to be a double vector of lenght 1 indicating in how many columns the graphics will be displayed.
#'
#' @details plots a complex graphic output displaying the data distribution within each sample
#'
#' @return a plot generated using ggplot2 is generated with this function.
#'
#' @author Geremy Clair
#' @export
#'
distribHistogram<- function(romics_object, bin=1, scale_y=c(0,100), scale_x=c(-10,10), col=3){
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(bin)){bin<-1}
  if(missing(col)){col<-3}
  if(!missing(scale_x) & !(is.numeric(scale_x) | is.double(scale_x)) & length(scale_x)!=2) {stop("scale_x has to be a numerical or double vector of lenght 2")}
  if(!missing(scale_y) & !(is.numeric(scale_y) | is.double(scale_y)) & length(scale_y)!=2) {stop("scale_y has to be a numerical or double vector of lenght 2")}
  if(missing(scale_x)){scale_x=c(-10,10)}
  if(missing(scale_y)){scale_y=c(0,100)}

  data<-as.matrix(romics_object$data)
  fill<-as.character(t(romics_object$metadata[grepl("colors_romics",rownames(romics_object$metadata)),]))

  #Check X bins and modify accordingly

  if(!missing(scale_x) & scale_x[2]< max(data, na.rm = TRUE)){
    scale_x[2]<-ceiling(max(data, na.rm = TRUE))
    warning("Your scale_x maximum value was too low, it was adjusted")
  }
  if(!missing(scale_x) & scale_x[1]> min(data, na.rm = TRUE)){
    scale_x[1]<- floor(min(data, na.rm = TRUE))
    warning("Your scale_x minimum value was too high, it was adjusted")
  }

  #check Y bins and modify accodingly (only top value)
  #estimate the break histogram
  breaks<- seq(scale_x[1],scale_x[2],by=bin)
  #create frequency tables

  freq_table <-data.frame(hist(data[!is.na(data[1]),1],breaks,plot = FALSE)$counts)
  for (i in 2:ncol(data)){
    freq_table<- cbind(freq_table,hist(data[!is.na(data[i]),i],breaks,plot = FALSE)$counts)
  }

  colnames(freq_table)<- colnames(data)
  freq_table$bins<-breaks[1:nrow(freq_table)]

  if(!missing(scale_y) & scale_y[1]> 0){
    warning("Your scale_y minimum was higher than 0, it was set at 0")
    scale_y[1]<-0
  }

  if(!missing(scale_y) & scale_y[2]<max(freq_table)){
    scale_y[2]<- max(freq_table)
    warning("Your scale_y maximum was too low, it was adjusted to allow to visualize all the data")
  }

  #generate the plots
  myplots <- list()
  for (i in 1:ncol(data))
    local({
      i <- i
      p <- ggplot()+aes(data[,i])+
        geom_histogram(binwidth=bin, fill=fill[i],alpha=I(.8))+
        ggtitle(names(data)[i])+theme_ROP()+
        scale_y_continuous(limits=scale_y)+
        scale_x_continuous(limits=scale_x)+
        ylab("Frequency")+
        ggtitle(colnames(romics_object$data)[i])

        if(romicsLogCheck(romics_object)==T)
        {if(sum(grepl("log2transform\\(",romics_object$steps))>0){
          p<-p+xlab("Log2(intensities)")}else{p<-p+xlab("Log10(intensities)")}}else{p<-p+xlab("Intensities")}

        myplots[[i]] <<- p
    })
  #generate the text of the plotting function
  plot_hist<-"grid.arrange("
  for (i in 1:(length(myplots))){
    plot_hist<-paste0(plot_hist,"myplots[[",i,"]],")
  }
  plot_hist<-paste0(plot_hist,"ncol=",col,")")

  #run the script stored in plot_hist
  options(warn = -1)
  p<-eval(parse(text = plot_hist))
  options(warn = 1)

}
#' distribHistogram
#' This function will display a grid.arranged histogram plot showing the distribution of the data in each sample
#'
#' @param romics_object has to be a log transformed romics_object created using romicsCreateObject() and transformed using the function log2transform() or log10transform()
#' @param bin has to be a numerical value indicating the width of the frequency bins to use for the visualization
#'
#' @details plots a complex graphic output displaying the data distribution within each sample
#'
#' @return a plot generated using ggplot2 is generated with this function.
#'
#' @author Geremy Clair
#' @export
#'
distribHistogramGlobal<-function(romics_object,bin=1){
  if(missing(romics_object)){stop("romics_object is missing")}
  if(class(romics_object)!="romics_object"){stop("your romics_object was not created using the function romicsCreateObject")}
  if(romics_object$steps[1]!="romics_object"){stop("romics_object is not in the appropriate format")}
  if(missing(bin)){bin=1}

  melted_data<-melt(romics_object$data)
  info_data<-summary(melted_data$value)

  p<- ggplot(melted_data, aes(x = value),fill="gray") +
    geom_histogram(position="identity", alpha=0.8,binwidth = bin)+
    xlab("data distribution")+
    theme_ROP()+
    scale_fill_manual(values=c("gray"))+
    ylab("Frequency")+
    ggtitle("Distribution Frequency of the data")

  if(romicsLogCheck(romics_object)==T)
  {if(sum(grepl("log2transform\\(",romics_object$steps))>0){
    p<-p+xlab("Log2(intensities)")}else{p<-p+xlab("Log10(intensities)")}}else{p<-p+xlab("Intensities")}

  return(p)
  }
