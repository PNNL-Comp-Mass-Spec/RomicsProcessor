
#' distribBoxplot()
#' @description Plots the boxplot of the sample distribution using ggplot2. The colors used for the plotting will correspond to the main_factor of the romics_object.
#' @param romics_object has to be a romics_object
#' @param aggregate_by has to be either NA (by default, in this case a Violin will be generated for each sample), or a factor from a romics_object, the list of factor can be found using the function romicsFactorNames()
#' @details create a ggplot2 graphic output displaying the normalized or not intensites (or logged intensities) within each sample. ggplot2 methods can be used to visually adjust the plot.
#' @return a plot generated using ggplot2 is generated with this function.
#' @author Geremy Clair
#' @export
distribBoxplot <- function(romics_object,aggregate_by="factor"){
  if(!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")}
  if(missing(aggregate_by)){
    if(ncol(romics_object$data)>100){
      warning("Over 100 samples were found in the romics_object,the object was aggregated by the main factor. If you want to see the distribution for each column use <NA> for the parameter <aggregate_by>")
      aggregate_by<-romics_object$main_factor
      }else{aggregate_by<-NA}
    }
  if(!is.na(aggregate_by) & !aggregate_by %in% romicsFactorNames(romics_object)){
    stop("<aggregate_by> has to either be NA or a factor of your romics_object, list of factors can be found using the function romicsFactorNames().")
  }

  if(!is.na(aggregate_by)){
    original_object<-romics_object
    romics_object<-romicsChangeFactor(romics_object,main_factor = aggregate_by)}

  data_graph<-reshape2::melt(romics_object$data,id.vars=0)
  fill_graph<-romics_object$colors
  fill_graph<-fill_graph[!is.na(data_graph[,2])]
  data_graph<-data_graph[!is.na(data_graph[,2]),]
  data_graph<-cbind(data_graph,fill_graph)
  colors<-as.character(t(romics_object$metadata[romicsFactorNames(romics_object)=="colors_romics",]))
  factor_graph<-data.frame(romics_object$metadata[romicsFactorNames(romics_object)==romics_object$main_factor,])
  factor_graph<-data.frame(variable=colnames(factor_graph),factor=as.character(t(factor_graph[1,])))
  data_graph<-merge(data_graph,factor_graph,by="variable")

  if(is.na(aggregate_by)){
p<-ggplot(aes(x=variable, y=value),data=data_graph)+
  geom_boxplot(aes(fill=factor),alpha=0.8)+
  scale_fill_manual(values=levels(factor(data_graph$fill_graph, levels = unique(data_graph$fill_graph))))+
  ggtitle("Boxplot of the data distribution within each sample")+
  theme_ROP()+
  xlab("Sample")}else{
p<-ggplot(aes(x=factor, y=value),data=data_graph)+
  geom_boxplot(aes(fill=factor),alpha=0.8)+
  scale_fill_manual(values=levels(factor(data_graph$fill_graph, levels = unique(data_graph$fill_graph))))+
  ggtitle("Boxplot of the data distribution within each factor level")+
  theme_ROP()+
  xlab("factor")
romics_object<-original_object
  }

if(romicsLogCheck(romics_object)==T)
{if(sum(grepl("log2transform\\(",romics_object$steps))>0){
  p<-p+ylab("Log2(intensities)")}else{p<-p+ylab("Log10(intensities)")}}else{p<-p+ylab("Intensities")}

return(p)
}

#' distribViolin()
#' @description Plots the violin plot of the sample distribution using ggplot2. The colors used for the plotting will correspond to the main_factor of the romics_object.
#' @param romics_object has to be a romics_object
#' @param aggregate_by has to be either NA (by default, in this case a Violin will be generated for each sample), or a factor from a romics_object, the list of factor can be found using the function romicsFactorNames()
#' @details create a ggplot2 graphic output displaying the normalized or not intensites (or logged intensities) within each sample. ggplot2 methods can be used to visually adjust the plot.
#' @return a plot generated using ggplot2 is generated with this function.
#' @author Geremy Clair
#' @export
distribViolin <- function(romics_object,aggregate_by="factor"){
  if(!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")}
  if(missing(aggregate_by)){
    if(ncol(romics_object$data)>100){
      warning("Over 100 samples were found in the romics_object,the object was aggregated by the main factor. If you want to see the distribution for each column use <NA> for the parameter <aggregate_by>")
      aggregate_by<-romics_object$main_factor
    }else{aggregate_by<-NA}
  }
  if(!is.na(aggregate_by) & !aggregate_by %in% romicsFactorNames(romics_object)){
    stop("<aggregate_by> has to either be NA or a factor of your romics_object, list of factors can be found using the function romicsFactorNames().")
  }

  if(!is.na(aggregate_by)){
    original_object<-romics_object
    romics_object<-romicsChangeFactor(romics_object,main_factor = aggregate_by)}

  data_graph<-reshape2::melt(romics_object$data,id.vars=0)
  fill_graph<-romics_object$colors
  fill_graph<-fill_graph[!is.na(data_graph[,2])]
  data_graph<-data_graph[!is.na(data_graph[,2]),]
  data_graph<-cbind(data_graph,fill_graph)
  colors<-as.character(t(romics_object$metadata[romicsFactorNames(romics_object)=="colors_romics",]))
  factor_graph<-data.frame(romics_object$metadata[romicsFactorNames(romics_object)==romics_object$main_factor,])
  factor_graph<-data.frame(variable=colnames(factor_graph),factor=as.character(t(factor_graph[1,])))
  data_graph<-merge(data_graph,factor_graph,by="variable")

  if(is.na(aggregate_by)){
    p<-ggplot(aes(x=variable, y=value),data=data_graph)+
      geom_violin(aes(fill=factor),alpha=0.8)+
      scale_fill_manual(values=levels(factor(data_graph$fill_graph, levels = unique(data_graph$fill_graph))))+
      ggtitle("Boxplot of the data distribution within each sample")+
      theme_ROP()+
      xlab("Sample")}else{
        p<-ggplot(aes(x=factor, y=value),data=data_graph)+
          geom_violin(aes(fill=factor),alpha=0.8)+
          scale_fill_manual(values=levels(factor(data_graph$fill_graph, levels = unique(data_graph$fill_graph))))+
          ggtitle("Boxplot of the data distribution within each factor level")+
          theme_ROP()+
          xlab("factor")
        romics_object<-original_object
      }

  if(romicsLogCheck(romics_object)==T)
  {if(sum(grepl("log2transform\\(",romics_object$steps))>0){
    p<-p+ylab("Log2(intensities)")}else{p<-p+ylab("Log10(intensities)")}}else{p<-p+ylab("Intensities")}

  return(p)
  }

#' distribJitter()
#' @description Plots the Jitter plot of the sample distribution using ggplot2. The colors used for the plotting will correspond to the main_factor of the romics_object.
#' @param romics_object has to be a romics_object
#' @param aggregate_by has to be either NA (by default, in this case a Violin will be generated for each sample), or a factor from a romics_object, the list of factor can be found using the function romicsFactorNames()
#' @details create a ggplot2 graphic output displaying the normalized or not intensites (or logged intensities) within each sample. ggplot2 methods can be used to visually adjust the plot.
#' @return a plot generated using ggplot2 is generated with this function.
#' @author Geremy Clair
#' @export
distribJitter <- function(romics_object,aggregate_by="factor"){
  if(!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")}
  if(missing(aggregate_by)){
    if(ncol(romics_object$data)>100){
      warning("Over 100 samples were found in the romics_object,the object was aggregated by the main factor. If you want to see the distribution for each column use <NA> for the parameter <aggregate_by>")
      aggregate_by<-romics_object$main_factor
    }else{aggregate_by<-NA}
  }
  if(!is.na(aggregate_by) & !aggregate_by %in% romicsFactorNames(romics_object)){
    stop("<aggregate_by> has to either be NA or a factor of your romics_object, list of factors can be found using the function romicsFactorNames().")
  }

  if(!is.na(aggregate_by)){
    original_object<-romics_object
    romics_object<-romicsChangeFactor(romics_object,main_factor = aggregate_by)}

  data_graph<-reshape2::melt(romics_object$data,id.vars=0)
  fill_graph<-romics_object$colors
  fill_graph<-fill_graph[!is.na(data_graph[,2])]
  data_graph<-data_graph[!is.na(data_graph[,2]),]
  data_graph<-cbind(data_graph,fill_graph)
  colors<-as.character(t(romics_object$metadata[romicsFactorNames(romics_object)=="colors_romics",]))
  factor_graph<-data.frame(romics_object$metadata[romicsFactorNames(romics_object)==romics_object$main_factor,])
  factor_graph<-data.frame(variable=colnames(factor_graph),factor=as.character(t(factor_graph[1,])))
  data_graph<-merge(data_graph,factor_graph,by="variable")

  if(is.na(aggregate_by)){
    p<-ggplot(aes(x=variable, y=value),data=data_graph)+
      geom_jitter(aes(fill=factor),alpha=0.8)+
      scale_fill_manual(values=levels(factor(data_graph$fill_graph, levels = unique(data_graph$fill_graph))))+
      ggtitle("Boxplot of the data distribution within each sample")+
      theme_ROP()+
      xlab("Sample")}else{
        p<-ggplot(aes(x=factor, y=value),data=data_graph)+
          geom_jitter(aes(fill=factor),alpha=0.8)+
          scale_fill_manual(values=levels(factor(data_graph$fill_graph, levels = unique(data_graph$fill_graph))))+
          ggtitle("Boxplot of the data distribution within each factor level")+
          theme_ROP()+
          xlab("factor")
        romics_object<-original_object
      }

  if(romicsLogCheck(romics_object)==T)
  {if(sum(grepl("log2transform\\(",romics_object$steps))>0){
    p<-p+ylab("Log2(intensities)")}else{p<-p+ylab("Log10(intensities)")}}else{p<-p+ylab("Intensities")}

  return(p)
}

#' distribSina()
#' @description Plots the Sina plot of the sample distribution using ggplot2. The colors used for the plotting will correspond to the main_factor of the romics_object.
#' @param romics_object has to be a romics_object
#' @param aggregate_by has to be either NA (by default, in this case a Violin will be generated for each sample), or a factor from a romics_object, the list of factor can be found using the function romicsFactorNames()
#' @details create a ggplot2 graphic output displaying the normalized or not intensites (or logged intensities) within each sample. ggplot2 methods can be used to visually adjust the plot.
#' @return a plot generated using ggplot2 is generated with this function.
#' @author Geremy Clair
#' @export
distribSina <- function(romics_object){
  if(!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")}
  if(missing(aggregate_by)){
    if(ncol(romics_object$data)>100){
      warning("Over 100 samples were found in the romics_object,the object was aggregated by the main factor. If you want to see the distribution for each column use <NA> for the parameter <aggregate_by>")
      aggregate_by<-romics_object$main_factor
    }else{aggregate_by<-NA}
  }
  if(!is.na(aggregate_by) & !aggregate_by %in% romicsFactorNames(romics_object)){
    stop("<aggregate_by> has to either be NA or a factor of your romics_object, list of factors can be found using the function romicsFactorNames().")
  }

  if(!is.na(aggregate_by)){
    original_object<-romics_object
    romics_object<-romicsChangeFactor(romics_object,main_factor = aggregate_by)}

  data_graph<-reshape2::melt(romics_object$data,id.vars=0)
  fill_graph<-romics_object$colors
  fill_graph<-fill_graph[!is.na(data_graph[,2])]
  data_graph<-data_graph[!is.na(data_graph[,2]),]
  data_graph<-cbind(data_graph,fill_graph)
  colors<-as.character(t(romics_object$metadata[romicsFactorNames(romics_object)=="colors_romics",]))
  factor_graph<-data.frame(romics_object$metadata[romicsFactorNames(romics_object)==romics_object$main_factor,])
  factor_graph<-data.frame(variable=colnames(factor_graph),factor=as.character(t(factor_graph[1,])))
  data_graph<-merge(data_graph,factor_graph,by="variable")

  if(is.na(aggregate_by)){
    p<-ggplot(aes(x=variable, y=value),data=data_graph)+
      geom_sina(aes(fill=factor),alpha=0.8)+
      scale_fill_manual(values=levels(factor(data_graph$fill_graph, levels = unique(data_graph$fill_graph))))+
      ggtitle("Boxplot of the data distribution within each sample")+
      theme_ROP()+
      xlab("Sample")}else{
        p<-ggplot(aes(x=factor, y=value),data=data_graph)+
          geom_sina(aes(fill=factor),alpha=0.8)+
          scale_fill_manual(values=levels(factor(data_graph$fill_graph, levels = unique(data_graph$fill_graph))))+
          ggtitle("Boxplot of the data distribution within each factor level")+
          theme_ROP()+
          xlab("factor")
        romics_object<-original_object
      }

  if(romicsLogCheck(romics_object)==T)
  {if(sum(grepl("log2transform\\(",romics_object$steps))>0){
    p<-p+ylab("Log2(intensities)")}else{p<-p+ylab("Log10(intensities)")}}else{p<-p+ylab("Intensities")}

  return(p)
  }

#' distribRidges()
#' @description Plots the ridges plot of the sample distribution using ggplot2. The colors used for the plotting will correspond to the main_factor of the romics_object.
#' @param romics_object has to be a romics_object
#' @param aggregate_by has to be either NA (by default, in this case a Violin will be generated for each sample), or a factor from a romics_object, the list of factor can be found using the function romicsFactorNames()
#' @details create a ggplot2 graphic output displaying the normalized or not intensites (or logged intensities) within each sample. ggplot2 methods can be used to visually adjust the plot.
#' @return a plot generated using ggplot2 is generated with this function.
#' @author Geremy Clair
#' @export
distribRidges <- function(romics_object){
  if(!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")}
  if(missing(aggregate_by)){
    if(ncol(romics_object$data)>100){
      warning("Over 100 samples were found in the romics_object,the object was aggregated by the main factor. If you want to see the distribution for each column use <NA> for the parameter <aggregate_by>")
      aggregate_by<-romics_object$main_factor
    }else{aggregate_by<-NA}
  }
  if(!is.na(aggregate_by) & !aggregate_by %in% romicsFactorNames(romics_object)){
    stop("<aggregate_by> has to either be NA or a factor of your romics_object, list of factors can be found using the function romicsFactorNames().")
  }

  if(!is.na(aggregate_by)){
    original_object<-romics_object
    romics_object<-romicsChangeFactor(romics_object,main_factor = aggregate_by)}

  data_graph<-reshape2::melt(romics_object$data,id.vars=0)
  fill_graph<-romics_object$colors
  fill_graph<-fill_graph[!is.na(data_graph[,2])]
  data_graph<-data_graph[!is.na(data_graph[,2]),]
  data_graph<-cbind(data_graph,fill_graph)
  colors<-as.character(t(romics_object$metadata[romicsFactorNames(romics_object)=="colors_romics",]))
  factor_graph<-data.frame(romics_object$metadata[romicsFactorNames(romics_object)==romics_object$main_factor,])
  factor_graph<-data.frame(variable=colnames(factor_graph),factor=as.character(t(factor_graph[1,])))
  data_graph<-merge(data_graph,factor_graph,by="variable")

  if(is.na(aggregate_by)){
    p<-ggplot(aes(x=variable, y=value),data=data_graph)+
      geom_density_ridges(aes(fill=factor),alpha=0.8)+
      scale_fill_manual(values=levels(factor(data_graph$fill_graph, levels = unique(data_graph$fill_graph))))+
      ggtitle("Boxplot of the data distribution within each sample")+
      theme_ROP()+
      xlab("Sample")}else{
        p<-ggplot(aes(x=factor, y=value),data=data_graph)+
          geom_density_ridges(aes(fill=factor),alpha=0.8)+
          scale_fill_manual(values=levels(factor(data_graph$fill_graph, levels = unique(data_graph$fill_graph))))+
          ggtitle("Boxplot of the data distribution within each factor level")+
          theme_ROP()+
          xlab("factor")
        romics_object<-original_object
      }

  if(romicsLogCheck(romics_object)==T)
  {if(sum(grepl("log2transform\\(",romics_object$steps))>0){
    p<-p+ylab("Log2(intensities)")}else{p<-p+ylab("Log10(intensities)")}}else{p<-p+ylab("Intensities")}

  return(p)
}

#' distribHistogram()
#' @description Plots individual frequency plots showing the discribution of the data in each sample using grid.Extra and ggplots2. Please use the function distribHistogramGlobal() to plot the global distribution of the whole dataset.
#' @param romics_object has to be a log transformed romics_object created using romicsCreateObject() and transformed using the function log2transform() or log10transform()
#' @param bin has to be a numerical value indicating the width of the frequency bins to use for the visualization
#' @param scale_y has to be a numerical/double vector of length 2 indicating the y_limits of each graph. if too low those will be automatically adjusted. by default  scale_y=c(0,100) will be used
#' @param scale_x has to be a numerical/double vector of length 2 indicating the x_limits of each graph. if too low those will be automatically adjusted. by default  scale_x=c(0,100) will be used
#' @param col has to be a double vector of lenght 1 indicating in how many columns the graphics will be displayed.
#' @details plots a complex graphic output displaying the data distribution within each sample
#' @return a plot generated using ggplot2 is generated with this function.
#' @author Geremy Clair
#' @export
distribHistogram<- function(romics_object, bin=1, scale_y=c(0,100), scale_x=c(-10,10), col=3){
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
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

#' distribFeature()
#' @description Plots the distribution of a specific feature across all samples, with options to display quartiles and threshold lines.
#' @param romics_object A romics_object created using romicsCreateObject()
#' @param feature_name Character string. Name of the feature to plot. Must match a row name in the data exactly.
#' @param bin Numeric. Width of the histogram bins. Default: 0.1
#' @param show_quantiles Logical. Whether to display quartile lines (Q1, median, Q3) on the plot. Default: TRUE
#' @param threshold_value Numeric. Value at which to draw a vertical threshold line. If NULL, no threshold line is drawn. Default: NULL
#' @param threshold_color Character. Color of the threshold line. Default: "red"
#' @param threshold_linetype Character. Line type for threshold ("solid", "dashed", "dotted"). Default: "dashed"
#' @param threshold_size Numeric. Size/width of the threshold line. Default: 1
#' @param scale_x Numeric vector of length 2. X-axis limits c(min, max). If NULL, uses data range. Default: NULL
#' @param scale_y Numeric vector of length 2. Y-axis limits c(min, max). If NULL, auto-scales. Default: NULL
#' @param fill_color Character. Fill color for histogram. If NULL, uses sample colors from romics_object. Default: NULL
#' @param alpha Numeric. Transparency of histogram bars (0-1). Default: 0.7
#' @param title Character. Plot title. If NULL, uses feature name. Default: NULL
#' @param show_stats Logical. Whether to display summary statistics on the plot. Default: TRUE
#' @details This function creates a histogram showing the distribution of a specific feature across all samples.
#' Quartile lines show Q1 (25th percentile), median (50th percentile), and Q3 (75th percentile).
#' The threshold line can be used to visualize cutoff values or reference points.
#' Automatic log-transform detection adjusts axis labels appropriately.
#' @return A ggplot2 object
#' @examples
#' \dontrun{
#' # Basic feature distribution
#' distribFeature(romics_obj, feature_name = "Gene1")
#'
#' # With threshold line at value 5
#' distribFeature(romics_obj, feature_name = "Gene1", threshold_value = 5)
#'
#' # Custom styling with no quartiles
#' distribFeature(romics_obj, feature_name = "Gene1",
#'                show_quantiles = FALSE, fill_color = "blue",
#'                threshold_value = 2.5, threshold_color = "green")
#' }
#' @author Geremy Clair
#' @export
distribFeature <- function(romics_object,
                           feature_name = NULL,
                           bin = 0.1,
                           show_quantiles = TRUE,
                           threshold_value = NULL,
                           threshold_color = "red",
                           threshold_linetype = "dashed",
                           threshold_size = 1,
                           scale_x = NULL,
                           scale_y = NULL,
                           fill_color = NULL,
                           alpha = 0.7,
                           title = NULL,
                           show_stats = TRUE) {

  # Input validation
  if(!is.romicsObject(romics_object) || missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }

  if(is.null(feature_name) || missing(feature_name)) {
    stop("feature_name must be provided")
  }

  if(!is.character(feature_name) || length(feature_name) != 1) {
    stop("feature_name must be a single character string")
  }

  if(!feature_name %in% rownames(romics_object$data)) {
    stop(paste("Feature", feature_name, "not found in romics_object data."))
  }

  if(!is.numeric(bin) || bin <= 0) {
    stop("bin must be a positive numeric value")
  }

  if(!is.null(threshold_value) && !is.numeric(threshold_value)) {
    stop("threshold_value must be numeric or NULL")
  }

  if(!is.null(scale_x) && (length(scale_x) != 2 || !is.numeric(scale_x))) {
    stop("scale_x must be a numeric vector of length 2 or NULL")
  }

  if(!is.null(scale_y) && (length(scale_y) != 2 || !is.numeric(scale_y))) {
    stop("scale_y must be a numeric vector of length 2 or NULL")
  }

  # Extract feature data
  feature_data <- as.numeric(romics_object$data[feature_name, ])
  feature_data_clean <- feature_data[!is.na(feature_data)]

  if(length(feature_data_clean) == 0) {
    stop("No non-missing values found for the specified feature")
  }

  # Determine x-axis limits
  if(is.null(scale_x)) {
    data_range <- range(feature_data_clean)
    buffer <- (data_range[2] - data_range[1]) * 0.05  # 5% buffer
    scale_x <- c(data_range[1] - buffer, data_range[2] + buffer)
  } else {
    # Validate and adjust scale_x if needed
    if(scale_x[1] > min(feature_data_clean)) {
      scale_x[1] <- floor(min(feature_data_clean))
      warning("scale_x minimum was too high and was adjusted")
    }
    if(scale_x[2] < max(feature_data_clean)) {
      scale_x[2] <- ceiling(max(feature_data_clean))
      warning("scale_x maximum was too low and was adjusted")
    }
  }

  # Calculate statistics
  feature_stats <- list(
    mean = mean(feature_data_clean),
    median = median(feature_data_clean),
    q1 = quantile(feature_data_clean, 0.25),
    q3 = quantile(feature_data_clean, 0.75),
    min = min(feature_data_clean),
    max = max(feature_data_clean),
    sd = sd(feature_data_clean),
    n = length(feature_data_clean),
    n_missing = sum(is.na(feature_data))
  )

  # Determine fill color
  if(is.null(fill_color)) {
    # Try to get color from romics object
    if("colors_romics" %in% rownames(romics_object$metadata)) {
      sample_colors <- as.character(t(romics_object$metadata["colors_romics", ]))
      fill_color <- sample_colors[1]  # Use first sample's color
    } else {
      fill_color <- "steelblue"  # Default color
    }
  }

  # Determine plot title
  if(is.null(title)) {
    title <- paste("Distribution of", feature_name)
  }

  # Determine x-axis label based on log transformation
  if(romicsLogCheck(romics_object)) {
    if(sum(grepl("log2transform\\(", romics_object$steps)) > 0) {
      x_label <- "Log2(intensities)"
    } else {
      x_label <- "Log10(intensities)"
    }
  } else {
    x_label <- "Intensities"
  }

  # Create the base plot
  p <- ggplot2::ggplot() +
    ggplot2::aes(feature_data_clean) +
    ggplot2::geom_histogram(binwidth = bin, fill = fill_color, alpha = alpha, color = "white", size = 0.2) +
    ggplot2::scale_x_continuous(limits = scale_x) +
    ggplot2::labs(
      title = title,
      x = x_label,
      y = "Frequency"
    ) +
    theme_ROP()

  # Add y-axis limits if specified
  if(!is.null(scale_y)) {
    p <- p + ggplot2::scale_y_continuous(limits = scale_y)
  }

  # Add quartile lines if requested
  if(show_quantiles) {
    p <- p +
      ggplot2::geom_vline(xintercept = feature_stats$q1, color = "blue", linetype = "dotted", size = 0.8) +
      ggplot2::geom_vline(xintercept = feature_stats$median, color = "blue", linetype = "solid", size = 1) +
      ggplot2::geom_vline(xintercept = feature_stats$q3, color = "blue", linetype = "dotted", size = 0.8) +
      ggplot2::annotate("text", x = feature_stats$q1, y = Inf, label = "Q1",
                        vjust = 1.2, hjust = -0.1, color = "blue", size = 3) +
      ggplot2::annotate("text", x = feature_stats$median, y = Inf, label = "Median",
                        vjust = 1.2, hjust = -0.1, color = "blue", size = 3) +
      ggplot2::annotate("text", x = feature_stats$q3, y = Inf, label = "Q3",
                        vjust = 1.2, hjust = -0.1, color = "blue", size = 3)
  }

  # Add threshold line if specified
  if(!is.null(threshold_value)) {
    p <- p +
      ggplot2::geom_vline(xintercept = threshold_value, color = threshold_color,
                          linetype = threshold_linetype, size = threshold_size) +
      ggplot2::annotate("text", x = threshold_value, y = Inf,
                        label = paste("Threshold:", round(threshold_value, 3)),
                        vjust = 1.2, hjust = 1.1, color = threshold_color, size = 3.5, fontface = "bold")
  }

  # Add summary statistics text if requested
  if(show_stats) {
    stats_text <- paste0(
      "n = ", feature_stats$n,
      " (", feature_stats$n_missing, " missing)",
      "\nMean = ", round(feature_stats$mean, 3),
      "\nSD = ", round(feature_stats$sd, 3),
      "\nRange = [", round(feature_stats$min, 3), ", ", round(feature_stats$max, 3), "]"
    )

    p <- p +
      ggplot2::annotation_custom(
        ggplot2::ggplotGrob(
          ggplot2::ggplot() +
            ggplot2::annotate("text", x = 0, y = 0, label = stats_text, hjust = 0, vjust = 1, size = 3) +
            ggplot2::theme_void()
        ),
        xmin = scale_x[1] + (scale_x[2] - scale_x[1]) * 0.7,
        xmax = scale_x[2],
        ymin = -Inf,
        ymax = Inf
      )
  }

  # Print summary statistics
  message(paste("Feature:", feature_name))
  message(paste("Samples:", feature_stats$n, "(", feature_stats$n_missing, "missing)"))
  message(paste("Mean:", round(feature_stats$mean, 3)))
  message(paste("Median:", round(feature_stats$median, 3)))
  message(paste("Q1:", round(feature_stats$q1, 3), "| Q3:", round(feature_stats$q3, 3)))
  message(paste("Range: [", round(feature_stats$min, 3), ",", round(feature_stats$max, 3), "]"))
  if(!is.null(threshold_value)) {
    below_threshold <- sum(feature_data_clean < threshold_value, na.rm = TRUE)
    above_threshold <- sum(feature_data_clean >= threshold_value, na.rm = TRUE)
    message(paste("Threshold", threshold_value, ":", below_threshold, "below,", above_threshold, "above"))
  }

  return(p)
}

#' distribHistogramGlobal()
#' @description Plots the frequency plot showing the discribution of the data in the whole data layer using ggplots2. Please use the function distribHistogram() to plot the global distribution for each sample.
#' @param romics_object has to be a log transformed romics_object created using romicsCreateObject() and transformed using the function log2transform() or log10transform()
#' @param bin has to be a numerical value indicating the width of the frequency bins to use for the visualization
#' @details plots a complex graphic output displaying the data distribution within each sample
#' @return a plot generated using ggplot2 is generated with this function.
#' @author Geremy Clair
#' @export
distribHistogramGlobal<-function(romics_object,bin=1,show_quantiles=T){
  if(missing(romics_object)){stop("romics_object is missing")}
  if(!inherits(romics_object, "romics_object")){stop("your romics_object was not created using the function romicsCreateObject")}
  if(romics_object$steps[1]!="romics_object"){stop("romics_object is not in the appropriate format")}
  if(missing(bin)){bin=1}
  if(missing(show_quantiles)){show_quantiles=T}

  melted_data<-reshape2::melt(romics_object$data)
  q<-quantile(melted_data$value,na.rm = T)

  info_data<-summary(melted_data$value)

  p<- ggplot(melted_data, aes(x = value),fill="gray") +
    geom_histogram(position="identity", alpha=0.8,binwidth = bin)+
    xlab("data distribution")+
    theme_ROP()+
    scale_fill_manual(values=c("gray"))+
    ylab("Frequency")+
    ggtitle("Distribution Frequency of the data")

  if(show_quantiles==T){
    p<-p+geom_vline(xintercept = q[2],col="orange",linetype = "dashed")
    p<-p+annotate("text",x=q[2],y=0,label=paste0("Q1(",round(q[2],2),")"),col="orange",hjust=-0.2,vjust=-0.5,angle=90)
    p<-p+geom_vline(xintercept = q[3],col="red",linetype = "dashed")
    p<-p+annotate("text",x=q[3],y=0,label=paste0("Q2(",round(q[3],2),")"),col="red",hjust=-0.2,vjust=-0.5,angle=90)
    p<-p+geom_vline(xintercept = q[4],col="darkred",linetype = "dashed")
    p<-p+annotate("text",x=q[4],y=0,label=paste0("Q3(",round(q[4],2),")"),col="darkred",hjust=-0.2,vjust=-0.5,angle=90)

  }


  if(romicsLogCheck(romics_object)==T)
  {if(sum(grepl("log2transform\\(",romics_object$steps))>0){
    p<-p+xlab("Log2(intensities)")}else{p<-p+xlab("Log10(intensities)")}}else{p<-p+xlab("Intensities")}

  return(p)
}

#' romicsCountDetect()
#' @description This function will detect in how many samples each feature is detected and create a summary data.frame from these counts
#' @param romics_object has to be a log transformed romics_object created using romicsCreateObject() and transformed using the function log2transform() or log10transform()
#' @return a data frame with the counts of proteins detected in different samples
#' @author Geremy Clair
#' @export
romicsCountDetect<-function(romics_object){
  if(!is.romicsObject(romics_object)){stop("Your <romics_object> has to be an romics_object created using the function romicsCreateObject().")}
  t<-data.frame(count=rowSums((romics_object$missingdata)))
  t$count<-ncol(romics_object$missingdata)-t
  Detected_in<-as.factor(t(t$count))
  Detected_in<-as.data.frame(t(table(Detected_in)))[-1]
  Detected_in$Percent<-Detected_in$Freq/nrow(romics_object$missingdata)*100
  return(Detected_in)}

#' romicsCountDetectPlot()
#' @description This function will detect in how many samples each feature is detected and create a plot from these counts
#' @param romics_object has to be a log transformed romics_object created using romicsCreateObject() and transformed using the function log2transform() or log10transform()
#' @return a ggplot plot
#' @author Geremy Clair
#' @export
romicsCountDetectPlot<-function(romics_object,percent=TRUE){
  t<-romicsCountDetect(romics_object)
  t$Percent_detected<-round(t$Percent,2)
  if(percent==TRUE){
    p<-ggplot(t,aes(x=Detected_in,y=Percent_detected))+geom_bar(stat="identity")+theme_ROP()+geom_text(aes(label=Percent_detected), vjust=-0.3, color="gray20", size=3.5)
  }else{
    p<-ggplot(t,aes(x=Detected_in,y=Freq))+geom_bar(stat="identity")+theme_ROP()+geom_text(aes(label=Freq), vjust=-0.3, color="gray20", size=3.5)
  }
  return(p)
}
