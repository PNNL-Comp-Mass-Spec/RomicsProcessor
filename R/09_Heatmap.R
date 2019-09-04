#' romicsHeatmap()
#' This function will plot a scaled heatmap of the data layer from a romics_object. This data can (or not be filtered) based on the statistics layer of the romics_object.
#'
#' @param romics_object A romics_object created using romicsCreateObject().
#' @param color_palette Character vector of colors. By default the viridis(n=20) will be used
#' @param color_boundaries Numerical vector of length 2. Indicates the min and max of the color scale. By default will be c(-2,2)
#' @param sample_hclust Boolean. Indicates if hclust has to be done for the samples.
#' @param sample_hclust_method_dist Sample dist method to be used. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman" or "kendall".
#' @param sample_hclust_method_hclust Sample agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward", "single", "complete", "average", "mcquitty", "median" or "centroid".
#' @param variable_hclust Boolean. Indicates if hclust has to be done for the variables.
#' @param variable_hclust_method_dist Variable dist method to be used. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman" or "kendall".
#' @param variable_hclust_method_hclust Variable agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward", "single", "complete", "average", "mcquitty", "median" or "centroid".
#' @param variable_hclust_number Numerical of length 1. Indicates the number of clusters to be used for the coloring of the variable hclust
#' @param ANOVA_filter Either 'none', 'p' or 'padj'. Indicates if an the ANOVA filter has to be used to plot the Heatmap (only the proteins below the filter will be displayed on the heatmap)
#' @param p Numerical of length 1 indicating the value of the ANOVA_filter cutoff (anything below this value will be conserved).
#' @param statCol_filter A column contained in the statistical layer of the romics_object, the list of columns can be obtained by using the function romicsCalculatedStats().
#' @param statCol_text Character to indicate how this column should be filtered (e.g. '<=0.05','>0.05','==1', '==TRUE', '>2')
#' @param statCol_filter2 A column contained in the statistical layer of the romics_object, the list of columns can be obtained by using the function romicsCalculatedStats().
#' @param statCol_text2 Character to indicate how this column should be filtered (e.g. '<=0.05','>0.05','==1', '==TRUE', '>2')
#' @param statCol_filter3 A column contained in the statistical layer of the romics_object, the list of columns can be obtained by using the function romicsCalculatedStats().
#' @param statCol_text3 Character to indicate how this column should be filtered (e.g. '<=0.05','>0.05','==1', '==TRUE', '>2')
#' @param notecol (optional) character string specifying the color for cellnote text. Defaults to "black".
#' @param density.info character string indicating whether to superimpose a 'histogram', a 'density' plot, or no plot ('none') on the color-key.
#' @param trace has to be in c("column","row","both","none"). See trace() documentation.
#' @param LabRow character vectors with row labels to use; these default to rownames(x)
#' @param cexCol positive numbers, used as cex.axis in for the row or column axis labeling. The defaults currently only use number of columns.
#' @param margin numeric vector of length 2 containing the margins (see par(mar= *)) for column and row names, respectively.
#' @param key.title main title of the color key. If set to NA no title will be plotted.
#' @param key.xlab x axis label of the color key. If set to NA no label will be plotted.
#' @param ...
#' @details Create a customizable and filterable heatmap based on the romics_object statistics Layer. the ANOVA filter enables to restrict the variable displayed to be only the ones passing an ANOVA
#' @details 3 stat column filters (StatCol) can be set simultaneously to restrict the variable displayed. Each filter enable to sort based on a given column of the statistics layer (statCol_filter) of an romics_object (the list of columns can be obtained by using the function romicsCalculatedStats()) using a specific text (statCol_text) this text indicate what parameter should be used to filter this column (example: column has to be positive -> statCol_text= '>0'). Note that the ANOVA filter is applied first (if any) and then the filters are applied sequencially (first, then second, then third).
#' @return A heatmap generated using the gplots::heatmap.2() function. Subsequently any adjustment of the heatmap.2 can be performed as described in gplots::heatmap.2() documentation.
#'
#' @author Geremy Clair
#' @export
#'
romicsHeatmap<-function(romics_object,
                        color_palette =viridis(20),
                        color_boundaries =c(-2,2),
                        sample_hclust=TRUE,
                        sample_hclust_method_dist = "euclidean",
                        sample_hclust_method_hclust = "ward.D",
                        variable_hclust=TRUE,
                        variable_hclust_method_dist= "euclidean",
                        variable_hclust_method_hclust= "ward.D",
                        variable_hclust_number = 1,
                        ANOVA_filter="none",
                        p=0.05,
                        statCol_filter="none",
                        statCol_text="<=0.05",
                        statCol_filter2="none",
                        statCol_text2="<=0.05",
                        statCol_filter3="none",
                        statCol_text3="<=0.05",
                        notecol="black",
                        density.info="none",
                        trace = "none",
                        labRow = FALSE,
                        cexCol=1,
                        margins = c(15, 5),
                        key.title = "Scaled Heatmap",
                        key.xlab = "Z-scores",
                        ...){

  arguments<<-as.list(match.call())
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(color_palette)){color_palette =viridis::viridis(20)}
  if(missing(color_boundaries)){color_boundaries=c(-2,2)}
  if(!is.numeric(color_boundaries) | length(color_boundaries)!=2){"<color_boundaries> has to be numerical of length 2."}
  if(missing(sample_hclust) | !is.logical(sample_hclust)){sample_hclust=TRUE}
  if(missing(variable_hclust) | !is.logical(variable_hclust)){sample_hclust=TRUE}
  if(missing(sample_hclust_method_dist)){sample_hclust_method_dist<-"euclidean"}
  if(missing(variable_hclust_method_dist)){variable_hclust_method_dist<-"euclidean"}
  if(missing(sample_hclust_method_hclust)){sample_hclust_method_hclust<-"ward.D"}
  if(missing(variable_hclust_method_hclust)){variable_hclust_method_hclust<-"ward.D"}
  if(!sample_hclust_method_dist %in% c("euclidean","maximum","manhattan","canberra","binary","minkowski")){stop("The sample_hclust_method_dist used is not appropriate, please check the documentation of the function.")}
  if(!variable_hclust_method_dist %in% c("euclidean","maximum","manhattan","canberra","binary","minkowski")){stop("The variable_hclust_method_dist used is not appropriate, please check the documentation of the function.")}
  if(!sample_hclust_method_hclust %in% c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")){stop("The sample_hclust_method_hclust used is not appropriate, , please check the documentation of the function.")}
  if(!variable_hclust_method_hclust %in% c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")){stop("The variable_hclust_method_hclust used is not appropriate, please check the documentation of the function.")}
  if(missing(variable_hclust_number)){variable_hclust_number=1}
  if(!is.numeric(variable_hclust_number)&variable_hclust_number<0){stop("variable_hclust_number should be an integer higher than 0")}
  if(missing(ANOVA_filter)){ANOVA_filter="none"}
  if(!ANOVA_filter %in% c("none","p","padj")){stop("Your ANOVA filter should be either none, p or padj")}
  if(missing(p)){p<-0.05}
  if(missing(statCol_filter)){statCol_filter="none"}
  if(!statCol_filter %in% c("none",romicsCalculatedStats(romics_object))){stop("<statCol_filter> is not a calculated stat column of the romics_object. Use the function romicsCalculatedStats() to identify the usable columns.")}
  if(missing(statCol_filter2)){statCol_filter="none"}
  if(!statCol_filter2 %in% c("none",romicsCalculatedStats(romics_object))){stop("<statCol_filter2> is not a calculated stat column of the romics_object. Use the function romicsCalculatedStats() to identify the usable columns.")}
  if(missing(statCol_filter3)){statCol_filter="none"}
  if(!statCol_filter3 %in% c("none",romicsCalculatedStats(romics_object))){stop("<statCol_filter3> is not a calculated stat column of the romics_object. Use the function romicsCalculatedStats() to identify the usable columns.")}

  #extract the data
  data<-data.frame(romics_object$data)

  #Filter based on ANOVA
  if(ANOVA_filter=="p"){
    if(is.null(romics_object$statistics$ANOVA_p)){
      warning("The ANOVA_p has not been calculated, no filtering was applied")}else{
        data <- data[romics_object$statistics$ANOVA_p<p,]
      }}
  if(ANOVA_filter=="padj"){
    if(is.null(romics_object$statistics$ANOVA_padj)){
      warning("The ANOVA_padj has not been calculated, no filtering was applied")}else{
        data <- data[romics_object$statistics$ANOVA_padj<p,]
      }}

  #Filter based on statCol_filter
  if(statCol_filter!="none"){
    if(missing(statCol_text)){
      warning("The stat <statCol_text> was missing the stat column was not filtered")
     }else{
      text<-paste0(arguments[[2]],"$statistics$`",statCol_filter,"`",statCol_text)
      message("The rows passing the following test were conserved for the heatmap:")
      message(text)
      test<-eval(parse(text=text))
      if(!is.logical(test) | length(test)!= nrow(romics_object$statistics)){
        message("The result of following test is not a logical vector of the same length as the dataset :")
        message(test)
      }else{data <- data[test,]}
    }}

  if(statCol_filter2!="none"){
    if(missing(statCol_text2)){
      warning("the stat <statCol_text2> was missing the stat column was not filtered")
    }else{
      text<-paste0(arguments[[2]],"$statistics$`",statCol_filter2,"`",statCol_text2)
      message("The rows passing the following test were conserved for the heatmap:")
      message(text)
      test<-eval(parse(text=text))
      if(!is.logical(test) | length(test)!= nrow(romics_object$statistics)){
        message("The result of following test is not a logical vector of the same length as the dataset :")
        message(test)
      }else{data <- data[test,]}
    }}

  if(statCol_filter3!="none"){
    if(missing(statCol_text3)){
      warning("the stat <statCol_text3> was missing the stat column was not filtered")
    }else{
      text<-paste0(arguments[[2]],"$statistics$`",statCol_filter3,"`",statCol_text3)
      message("The rows passing the following test were conserved for the heatmap:")
      message(text)
      test<-eval(parse(text=text))
      if(!is.logical(test) | length(test)!= nrow(romics_object$statistics)){
        message("The result of following test is not a logical vector of the same length as the dataset :")
        message(test)
      }else{data <- data[test,]}
    }}

  #the variable_hclust_number as to be a round number
  variable_hclust_number<-as.integer(variable_hclust_number)

  #scale data
  scaled_data <- t(scale(t(data), center = TRUE, scale = TRUE))

  #sample dendrogram
  if(sample_hclust==TRUE){
    #extract the colors
    colors_dend <- as.character(t(romics_object$metadata[nrow(romics_object$metadata),]))
    #define the clustering method
    distance_samples<-dist(t(data), method=sample_hclust_method_dist)
    #define the hclust agglomeration method
    hc<- hclust(distance_samples, method=sample_hclust_method_hclust)
    #convert the hclust into a dendrogram
    sample_dd <- as.dendrogram(hc)
    #order the colors as the hclust
    colors_dend<-colors_dend[order=hc$order]
    #colors branches
    sample_dd <- color_branches(sample_dd,k = NULL, h = NULL,col=as.vector(colors_dend))
    } else {sample_dd=FALSE}


#variable hclust
  if(variable_hclust==TRUE){
    #define the clustering method
    dist_variable<-dist(scaled_data, method=variable_hclust_method_dist)
    #define the hclust agglomeration method
    variable_hc<- hclust(dist_variable, method=variable_hclust_method_hclust)
    #convert the hclust into a dendrogram
    variable_dd <- as.dendrogram(variable_hc)
    #colors sample_dd
    variable_dd <- color_branches(variable_dd, k=variable_hclust_number )
    variable_dd <- color_labels(variable_dd, k=variable_hclust_number )
  }else{variable_dd=FALSE}

  #make the color palette
  palette <- colorRampPalette(color_palette)(n = 299)
  #define the color palette breaks based on the minimum and maximum
  col_breaks = c(seq(color_boundaries[1],color_boundaries[2],length=300))

  if(missing(notecol)){notecol="black"} # change font color of cell labels to black
  if(missing(density.info)){density.info="none"} # turns off density plot inside color legend
  if(missing(trace)){trace = "none"}  # turns off trace lines inside the heat map
  if(missing(labRow)){labRow = FALSE} #remove row labels
  if(missing(cexCol)){cexCol=1} #size col text
  if(missing(margins)){(margins = c(15, 5))} #margin (bottom, right)
  if(missing(key.title)){key.title = "Scaled Heatmap"} #change title
  if(missing(key.xlab)){key.xlab = "Z-scores"} #change title

  #make heatmap
  heatmap.2(scaled_data,
            notecol=notecol,
            density.info=density.info,
            trace=trace,
            col=palette,
            breaks=col_breaks,
            Rowv = variable_dd,
            Colv = sample_dd,
            labRow = labRow,
            margins=margins,
            cexCol=cexCol,
            keysize = 1,
            key.title = key.title,
            key.xlab = key.xlab,
            ...)
}


#' romicsVariableHclust
#' This function will plot a hierarchical clustering for the variables. It will also add two new column in the statistical layer of the romics object indicating the order of the clustering and the clusters to which the variables belong.
#'
#' @param romics_object A romics_object created using romicsCreateObject().
#' @param clusters Numerical of length 1. Indicates the number of clusters to be used for the coloring of the variable hclust. 8 by default.
#' @param ANOVA_filter Either 'none', 'p' or 'padj'. Indicates if an the ANOVA filter has to be used to plot the Heatmap (only the proteins below the filter will be displayed on the heatmap)
#' @param p Numerical of length 1 indicating the value of the ANOVA_filter cutoff (anything below this value will be conserved).
#' @param method_dist Dist method to be used. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman" or "kendall".
#' @param method_hclust Agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward", "single", "complete", "average", "mcquitty", "median" or "centroid".
#' @param statCol_filter A column contained in the statistical layer of the romics_object, the list of columns can be obtained by using the function romicsCalculatedStats().
#' @param statCol_text Character to indicate how this column should be filtered (e.g. '<=0.05','>0.05','==1', '==TRUE', '>2')
#' @param statCol_filter2 A column contained in the statistical layer of the romics_object, the list of columns can be obtained by using the function romicsCalculatedStats().
#' @param statCol_text2 Character to indicate how this column should be filtered (e.g. '<=0.05','>0.05','==1', '==TRUE', '>2')
#' @param statCol_filter3 A column contained in the statistical layer of the romics_object, the list of columns can be obtained by using the function romicsCalculatedStats().
#' @param statCol_text3 Character to indicate how this column should be filtered (e.g. '<=0.05','>0.05','==1', '==TRUE', '>2')
#' @param plot Boolean indicating wether the clustering should be plotted or not.
#' @param scale_data Boolean indicating wether the data should be scaled or not.
#' @details Create a customizable and filterable hierarchical clustering based on the romics_object statistics Layer. The ANOVA filter enables to restrict the variable displayed to be only the ones passing an ANOVA
#' @details 3 stat column filters (StatCol) can be set simultaneously to restrict the variable displayed. Each filter enable to sort based on a given column of the statistics layer (statCol_filter) of an romics_object (the list of columns can be obtained by using the function romicsCalculatedStats()) using a specific text (statCol_text) this text indicate what parameter should be used to filter this column (example: column has to be positive -> statCol_text= '>0'). Note that the ANOVA filter is applied first (if any) and then the filters are applied sequencially (first, then second, then third).
#' @return return an Romics_object containing two new columns in its statistical layer : hclust_cluster, hclust_order. if those columns were pre-existing they'll be replaced
#'
#' @author Geremy Clair
#' @export
#'
romicsVariableHclust<-function(romics_object, clusters=8, ANOVA_filter="none", p=0.05, method_dist = "euclidean", method_hclust = "ward.D", plot=TRUE, scale_data=TRUE, statCol_filter="none", statCol_text="<=0.05", statCol_filter2="none", statCol_text2="<=0.05", statCol_filter3="none", statCol_text3="<=0.05"){
  arguments<-as.list(match.call())
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(method_dist)){sample_hclust_method_dist<-"euclidean"}
  if(missing(method_hclust)){sample_hclust_method_hclust<-"ward.D"}
  if(!method_dist %in% c("euclidean","maximum","manhattan","canberra","binary","minkowski")){stop("The sample_hclust_method_dist used is not appropriate, please check the documentation of the function.")}
  if(!method_hclust %in% c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")){stop("The sample_hclust_method_hclust used is not appropriate, , please check the documentation of the function.")}
  if(missing(clusters)){clusters=8}
  if(!is.numeric(clusters)&clusters<0){stop("variable_hclust_number should be an integer higher than 0")}
  if(missing(ANOVA_filter)){ANOVA_filter="none"}
  if(!ANOVA_filter %in% c("none","p","padj")){stop("Your ANOVA filter should be either none, p or padj")}
  if(missing(p)){p<-0.05}
  if(missing(statCol_filter)){statCol_filter="none"}
  if(!statCol_filter %in% c("none",romicsCalculatedStats(romics_object))){stop("<statCol_filter> is not a calculated stat column of the romics_object. Use the function romicsCalculatedStats() to identify the usable columns.")}
  if(missing(statCol_filter2)){statCol_filter="none"}
  if(!statCol_filter2 %in% c("none",romicsCalculatedStats(romics_object))){stop("<statCol_filter2> is not a calculated stat column of the romics_object. Use the function romicsCalculatedStats() to identify the usable columns.")}
  if(missing(statCol_filter3)){statCol_filter="none"}
  if(!statCol_filter3 %in% c("none",romicsCalculatedStats(romics_object))){stop("<statCol_filter3> is not a calculated stat column of the romics_object. Use the function romicsCalculatedStats() to identify the usable columns.")}
  if(missing(plot)){plot=TRUE}
  if(!is.logical(plot) | length(plot)!=1){stop("plot has to be either TRUE or FALSE")}
  if(missing(scale_data)){scale_data=TRUE}
  if(!is.logical(scale_data)| length(scale_data)!=1){stop("scale_data has to be either TRUE or FALSE")}

  data<-romics_object$data

  #Filter based on ANOVA
  if(ANOVA_filter=="p"){
    if(is.null(romics_object$statistics$ANOVA_p)){
      warning("The ANOVA_p has not been calculated, no filtering was applied")}else{
        data <- data[romics_object$statistics$ANOVA_p<p,]
      }}

  if(ANOVA_filter=="padj"){
    if(is.null(romics_object$statistics$ANOVA_p)){
      warning("The ANOVA_p has not been calculated, no filtering was applied")}else{
        data <- data[romics_object$statistics$ANOVA_padj<p,]
      }}

  #Filter based on statCol_filter
  if(statCol_filter!="none"){
    if(missing(statCol_text)){
      warning("the stat <statCol_text> was missing the stat column was not filtered")
    }else{
      text<-paste0(arguments[[2]],"$statistics$`",statCol_filter,"`",statCol_text)
      message("The rows passing the following test were conserved for the heatmap:")
      message(text)
      test<-eval(parse(text=text))
      if(!is.logical(test) | length(test)!= nrow(romics_object$statistics)){
        message("The result of following test is not a logical vector of the same length as the dataset :")
        message(test)
      }else{data <- data[test,]}
    }}

  if(statCol_filter2!="none"){
    if(missing(statCol_text2)){
      warning("the stat <statCol_text2> was missing the stat column was not filtered")
    }else{
      text<-paste0(arguments[[2]],"$statistics$`",statCol_filter2,"`",statCol_text2)
      message("The rows passing the following test were conserved for the heatmap:")
      message(text)
      test<-eval(parse(text=text))
      if(!is.logical(test) | length(test)!= nrow(romics_object$statistics)){
        message("The result of following test is not a logical vector of the same length as the dataset :")
        message(test)
      }else{data <- data[test,]}
    }}

  if(statCol_filter3!="none"){
    if(missing(statCol_text3)){
      warning("the stat <statCol_text3> was missing the stat column was not filtered")
    }else{
      text<-paste0(arguments[[2]],"$statistics$`",statCol_filter3,"`",statCol_text3)
      message("The rows passing the following test were conserved for the heatmap:")
      message(text)
      test<-eval(parse(text=text))
      if(!is.logical(test) | length(test)!= nrow(romics_object$statistics)){
        message("The result of following test is not a logical vector of the same length as the dataset :")
        message(test)
      }else{data <- data[test,]}
    }}

  #if data_scale == TRUE then Z score the data if not then keep the data as is.
  if (scale_data==TRUE){data<- t(scale(t(data), center = TRUE, scale = TRUE))}

  #define the clustering method
  dist_variable<-dist(data, method=method_dist)
  #define the hclust agglomeration method
  variable_hc<- hclust(dist_variable, method=method_hclust)
  #convert the hclust into a dendrogram
  variable_dd <- as.dendrogram(variable_hc)
  #colors sample_dd
  variable_dd <- color_branches(variable_dd, k=clusters )
  variable_dd <- color_labels(variable_dd, k=clusters )
  #calculate the clusters using cuttreee
  results_hc<- data.frame(cbind(hclust_clusters=cutree(variable_hc, k=clusters), hclust_order=variable_hc$order))

  #create an object containing matching same order as the original data
  columns_hc<-data.frame(matrix(ncol=2, nrow=nrow(romics_object$data)))
  rownames(columns_hc)<- rownames(romics_object$data)
  columns_hc<-results_hc[match(rownames(columns_hc),rownames(results_hc)),]
  rownames(columns_hc)<- rownames(romics_object$data)

  #if already exist replace if not create
  if("hclust_clusters" %in% colnames(romics_object$statistics)){
    romics_object$statistics<-romics_object$statistics[,colnames(romics_object$statistics)!="hclust_clusters"]
    warning("The hclust_clusters were previously calculated, the previous calculation was removed and replaced with the newer one")
  }

  if("hclust_order" %in% colnames(romics_object$statistics)){
    romics_object$statistics<-romics_object$statistics[,colnames(romics_object$statistics)!="hclust_order"]
    warning("The hclust_order were previously calculated, the previous calculation was removed and replaced with the newer one")
  }
  #add the columns to $statistics
  romics_object$statistics<-cbind(romics_object$statistics,columns_hc)
  print("The columns hclust_clusters and hclust_order were added to the statistics")

  if (plot==TRUE){
    plot(variable_dd)
  }

  #update Steps
  romics_object<-romicsUpdateSteps(romics_object,arguments)

  return(romics_object)
}
