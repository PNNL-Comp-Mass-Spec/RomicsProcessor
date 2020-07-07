#' log2transform()
#' @description log2-tranforms the romics_object data layer.
#' @param romics_object has to be an romics_object created using romicsCreateObject() that has not been previously log-transformed.
#' @details will log2 transform the romics object
#' @return This function returns the transformed romics_object with updated data layer
#' @author Geremy Clair
#' @export
#'
log2transform<-function(romics_object){
  arguments<-as.list(match.call())
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(romicsLogCheck(romics_object)){stop("Your romics_object was previously log transformed, see your romics_object$steps layer for more details")}

  romics_object$data<- data.frame(apply(romics_object$data,2,log2))

  romics_object<-romicsUpdateSteps(romics_object,arguments)
  return(romics_object)
}

#' log10transform()
#' @description log10-tranforms the romics_object data layer.
#' @param romics_object has to be an romics_object created using romicsCreateObject() that has not been previously log-transformed.
#' @details will log10 transform the romics object
#' @return This function returns the transformed romics_object with updated data layer
#' @author Geremy Clair
#' @export
log10transform<-function(romics_object){
  arguments<-as.list(match.call())
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(romicsLogCheck(romics_object)){stop("Your romics_object was previously log transformed, see the romics_object$steps layer for more details")}

  romics_object$data<- data.frame(apply(romics_object$data,2,log10))

  romics_object<-romicsUpdateSteps(romics_object,arguments)
  return(romics_object)
}

#' unlog2data()
#' @description Reverses the log2 tranformation of the romics_object data layer.
#' @param romics_object has to be an romics_object created using romicsCreateObject() that has not been previously log-transformed using log2transform()
#' @details will reverse the log2 transformation of the romics_object
#' @return This function returns the transformed romics_object with updated data layer
#' @author Geremy Clair
#' @export
unlog2data<-function(romics_object){
  arguments<-as.list(match.call())
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(!sum(grepl("log2transform\\(",romics_object$steps))>0){
    stop("Your romics_object was not previously log2 transformed using the function log2transform(). See the romics_object$steps layer for more details.")}
  if(sum(grepl("log2transform\\(",romics_object$steps))<=sum(grepl("unlog2data\\(",romics_object$steps))){
   stop("Your romics_object was already unlogged using the function unlog2data(). See the romics_object$steps layer for more details.")}

  romics_object$data<- 2^romics_object$data

  romics_object<-romicsUpdateSteps(romics_object,arguments)

  return(romics_object)

}

#' unlog10data()
#' @description Reverses the log10 tranformation of the romics_object data layer.
#' @param romics_object has to be an romics_object created using romicsCreateObject() that has not been previously log-transformed using log10transform()
#' @details will reverse the log10 transformation of the romics_object
#' @return This function returns the transformed romics_object with updated data layer
#' @author Geremy Clair
#' @export
unlog10data<-function(romics_object){
  arguments<-as.list(match.call())
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(!sum(grepl("log10transform\\(",romics_object$steps))>0){
    stop("Your romics_object was not previously log10 transformed using the function log10transform(). See the romics_object$steps layer for more details.")}
  if(sum(grepl("log10transform\\(",romics_object$steps))<=sum(grepl("unlog10data\\(",romics_object$steps))){
    stop("Your romics_object was already unlogged using the function unlog10data(). See the romics_object$steps layer for more details.")}

  romics_object$data<- 2^romics_object$data

  romics_object<-romicsUpdateSteps(romics_object,arguments)

  return(romics_object)
  }

#' MedianNormSample()
#' @description Normalizes the samples by their median. The median of the medians of all the samples is used as the alignment point.
#' @param romics_object has to be an romics_object created using romicsCreateObject() that has not been previously log-transformed using log10transform()
#' @details Median normalize within each sample the median of each sample median will be used as alignment point. If you waht to center the median at 0 please use the function medianCenterSample().
#' @return This function returns the transformed romics_object with updated data layer.
#' @author Geremy Clair
#' @export
medianNormSample<-function(romics_object){
  arguments<-as.list(match.call())
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}

  col_median<- apply(romics_object$data,2,function(x){median(x[!is.na(x)])})
  median_median<- median(col_median)

  if(romicsLogCheck(romics_object)==FALSE){
    for(i in 1:ncol(romics_object$data)){
      romics_object$data[,i]<- romics_object$data[,i]/col_median[i]*median_median
    }
  } else {
    for(i in 1:ncol(romics_object$data)){
      romics_object$data[,i]<- romics_object$data[,i]-col_median[i]+median_median
    }
  }

  romics_object<-romicsUpdateSteps(romics_object,arguments)

    return(romics_object)
}

#' medianCenterSample()
#' @description Normalizes the samples by their median. Zero is used as the median alignment center.
#' @param romics_object has to be an romics_object created using romicsCreateObject() that has not been previously log-transformed using log10transform()
#' @details Median normalize within each sample the median will be zero centered
#' @return This function returns the transformed romics_object with updated data layer.
#' @author Geremy Clair
#' @export
medianCenterSample<-function(romics_object){
  arguments<-as.list(match.call())
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}

  col_median<- apply(romics_object$data,2,function(x){median(x[!is.na(x)])})
  median_median<- median(col_median)

    if(romicsLogCheck(romics_object)==FALSE){
    for(i in 1:ncol(romics_object$data)){
      romics_object$data[,i]<- romics_object$data[,i]/col_median[i]
    }
  } else {
    for(i in 1:ncol(romics_object$data)){
      romics_object$data[,i]<- romics_object$data[,i]-col_median[i]
      }
  }

  romics_object<-romicsUpdateSteps(romics_object,arguments)

  #return object
  return(romics_object)
}

#' medianNormFactor()
#' @description Normalizes the samples by their median within a given factor. The median of the median within this factor will be used as factor-specific median center.
#' @param romics_object has to be an romics_object created using romicsCreateObject() that has not been previously log-transformed using log10transform()
#' @details median normalize the samples within a given factor the median of the median of this factor will be used as new factor median.
#' @return This function returns the transformed romics_object with updated data layer.
#' @author Geremy Clair
#' @export
medianNormFactor<-function(romics_object, main_factor= "factor"){
  if(missing(romics_object)){stop("romics_object is missing")}
  if(class(romics_object)!="romics_object"){stop("your romics_object was not created using the function romicsCreateObject")}
  if(romics_object$steps[1]!="romics_object"){stop("romics_object is not in the appropriate format")}
  if(missing(main_factor)){main_factor<-romics_object$main_factor}
  #import the data
  data<-romics_object$data

  #calculate the median for each column
  col_median<- apply(data,2,function(x){median(x[!is.na(x)])})

  #unique factors
  factor<- as.character(t(romics_object$metadata[rownames(romics_object$metadata)==romics_object$main_factor,]))
  unique_factor<-unique(factor)

  #calculate the median of the median within each factor
  median_factor<-as.numeric()

  for (i in 1:length(unique_factor)){
    median_factor[i]<-median(col_median[factor==unique_factor[i]])
  }

  #calculate the adjusting of the median for each sample (based on a given factor)
  adjust_median <- median_factor[match(factor,unique_factor)]

  #do the normalization using the - option if the data was log transformed or with the / if the data was not transformed
  if(sum(romics_object$steps=="log10transform()")+sum(romics_object$steps=="log2transform()")==0){
    print("Your data was median transformed")
    for(i in 1:ncol(data)){
      data[,i]<- data[,i]/col_median[i]*adjust_median[i]
    }
  } else {
    print("Your log transformed data was median transformed")
    for(i in 1:ncol(data)){
      data[,i]<- data[,i]-col_median[i]+adjust_median[i]
    }
  }

  #place this in your romics_object
    romics_object$data<-data
    #append the steps
    romics_object$steps<- c(romics_object$steps,paste0(gsub(" ","_",format(Sys.time(),"%b_%d_%Y_%X")),": the data was was median normalized for each column"))
    romics_object$steps<- c(romics_object$steps,"medianNormFactor()")
    #return object
    return(romics_object)

}
