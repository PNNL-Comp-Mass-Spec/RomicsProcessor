#' romics_outlier_eval()
#' @description Plots the the pmartR method for the evaluation of outliers (require the installation of the package pmartR). To remove the outlier below a certain treshold use the function Romics_outlier_remove()
#' @param romics_object A log transformed romics_object created using romicsCreateObject() and transformed using the function log2transform() or log10transform()
#' @param seed An integer of length 1, by default 42 will be used
#' @param metrics A character vector containing the following terms to indicate which parameters to use for the filtering of the data : 'Correlation', 'Proportion_Missing','MAD', 'Skewness'. By defaults all parameters will be used.
#' @param pvalue_threshold A numeric vector of lenght 1 indicating the pvalue threshold to be used.
#' @param label Either TRUE or FALSE to indicate if the labels have to be plotted.
#' @details This function requires the package 'pmartR' to be installed and loaded to be excecuted. It will calculate and plot the samples to be filtered out using the function Romics_outlier_eval().
#' @return This function will print the pmartR filtering details and will return 2 plots the first one is a scatter plot of the pvalue by log2(Robust Mahalanobis Distance) the second one is a scatter plot of the log2(Robust Mahalanobis Distance) per sample.
#' @references Matzke, M., Waters, K., Metz, T., Jacobs, J., Sims, A., Baric, R., Pounds, J., and Webb-Robertson, B.J. (2011), Improved quality control processing of peptide-centric LC-MS proteomics data. Bioinformatics. 27(20): 2866-2872.
#' @author Geremy Clair
#' @export
romicsOutlierEval<-function(romics_object, seed=42, metrics = c("Correlation", "Proportion_Missing", "MAD", "Skewness"), pvalue_threshold=0.01,label=TRUE){
  arguments<-as.list(match.call())
  if(!"pmartR" %in% rownames(installed.packages()) & !"package:pmartR"  %in% search()){stop("to run this function the package 'pmartR' has to be installed and loaded")}
  if(!"pmartR" %in% (.packages())){
    library("pmartR")
    print("pmartR was not loaded it was loaded to execute this function")
  }
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(metrics)){metrics = c("Correlation", "Proportion_Missing", "MAD", "Skewness")}
  if(missing(pvalue_threshold)){pvalue_threshold=0.01}
  if(missing(seed)){set.seed(Sys.time())}
  if(missing(label)){label=TRUE}

  set.seed(seed)
  Pmart_data<-romicsPmartR(romics_object)
  Pmart_data <- group_designation(omicsData = Pmart_data, main_effects = romics_object$main_factor, covariates = NULL)
  rmdfilter<-rmd_filter(Pmart_data, metrics=metrics)

  summary(rmdfilter, pvalue_threshold = pvalue_threshold)

  md_threshold<- max(rmdfilter$Log2.md[rmdfilter$pvalue>pvalue_threshold])

  plot1<-ggplot(rmdfilter,aes(x=pvalue,y=Log2.md))+geom_point(size=4, alpha=0.6,fill="gray75")+
    theme_ROP()+
    geom_hline(aes(yintercept=md_threshold),colour = "red")+
    ggtitle("pvalue vs. log2(Robust Mahalanobis Distance)")

  #force the order to stay the same
  rmdfilter$SampleID<-factor(rmdfilter$SampleID, levels = rmdfilter$SampleID[order(rmdfilter$Group)])
  rmdfilter<-rmdfilter[order(rmdfilter$Group),]

  colors<-unique(romics_object$colors)

  plot2<-ggplot(rmdfilter,aes(x=SampleID,y=Log2.md,color=Group))+
    geom_point(size=4, alpha= ifelse(rmdfilter$Log2.md > md_threshold, 0.8, 0.4))+
    geom_hline(aes(yintercept=md_threshold),colour = "red")+
    scale_color_manual(values=colors)+
    ylab("log2(Robust Mahalanobis Distance)")+
    theme_ROP()+
    ggtitle(paste0("Sample outlier result (p<",pvalue_threshold,")"))

  if (label==TRUE){plot2<-plot2+geom_text(label=rmdfilter$SampleID)}

  return(list(plot1,plot2))
}

#' Romics_outlier_remove()
#' @description Removes outliers below a certain pvalue threshold. The pmartR method for the removal of outliers is used (require the installation of the package pmartR). The evaluation of the outlier removal can be done prior to apply this function using the function Romics_outlier_eval().
#' @param romics_object A log transformed romics_object created using romicsCreateObject() and transformed using the function log2transform() or log10transform()
#' @param seed An integer of length 1, by default 42 will be used
#' @param metrics A character vector containing the following terms to indicate which parameters to use for the filtering of the data : 'Correlation', 'Proportion_Missing','MAD', 'Skewness'. By defaults all parameters will be used.
#' @param pvalue_threshold A numeric vector of lenght 1 indicating the pvalue threshold to be used.
#' @details This function requires the package 'pmartR' to be installed and loaded to be excecuted. It will calculate and plot the samples to be filtered out using the function Romics_outlier_eval().
#' @return This function will print the pmartR filtering details and will return 2 plots the first one is a scatter plot of the pvalue by log2(Robust Mahalanobis Distance) the second one is a scatter plot of the log2(Robust Mahalanobis Distance) per sample.
#' @references Matzke, M., Waters, K., Metz, T., Jacobs, J., Sims, A., Baric, R., Pounds, J., and Webb-Robertson, B.J. (2011), Improved quality control processing of peptide-centric LC-MS proteomics data. Bioinformatics. 27(20): 2866-2872.
#' @author Geremy Clair
#' @export
romicsOutlierRemove<-function(romics_object, seed=42, metrics = c("Correlation", "Proportion_Missing", "MAD", "Skewness"), pvalue_threshold=0.01){
  arguments<-as.list(match.call())
  if(!"pmartR" %in% rownames(installed.packages()) & !"package:pmartR"  %in% search()){stop("to run this function the package 'pmartR' has to be installed and loaded")}
  if(!"pmartR" %in% (.packages())){
    library("pmartR")
    print("pmartR was not loaded it was loaded to execute this function")
    }
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(metrics)){metrics = c("Correlation", "Proportion_Missing", "MAD", "Skewness")}
  if(missing(pvalue_threshold)){pvalue_threshold=0.01}
  if(missing(seed)){set.seed(Sys.time())}

  set.seed(seed)
  Pmart_data<-romicsPmartR(romics_object)
  Pmart_data <- group_designation(omicsData = Pmart_data, main_effects = romics_object$main_factor, covariates = NULL)
  rmdfilter<-rmd_filter(Pmart_data, metrics=metrics)

  rmdsummary<-summary(rmdfilter, pvalue_threshold = pvalue_threshold)

  romics_object$data<-romics_object$data[,!colnames(romics_object$data) %in% rmdsummary$filtered_samples]
  romics_object$metadata<-romics_object$metadata[,!colnames(romics_object$metadata) %in% rmdsummary$filtered_samples]
  romics_object$missingdata<-romics_object$missingdata[,!colnames(romics_object$missingdata) %in% rmdsummary$filtered_samples]
  romics_object<-romicsUpdateColor(romics_object)

  romics_object<-romicsUpdateSteps(romics_object,arguments)
  romics_object<-romicsAddDependency(romics_object,new_dependency = "pmartR")
  return(romics_object)
}
