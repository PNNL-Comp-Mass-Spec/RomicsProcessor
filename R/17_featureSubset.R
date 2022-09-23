#' featureSubset()
#' @description Replaces zeros in the romics_object by NA values
#' @param romics_object has to be an romics_object created using romicsCreateObject()
#' @param feature_list has to be a character vector containing features to be extracted
#' @details This function will create a new romics object in which only the features contained in the feature_list are contained.
#' @details The function supports partial matches and allows the features from the feature list not to be present in the romics_object
#' @return This function returns an romics_object
#' @author Geremy Clair
#' @export
featureSubset<-function(romics_object,feature_list=c("feature1","feature2","etc")){
  arguments<-as.list(match.call())
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("<romics_object> is missing or is not in the appropriate format")}
  if(!is.character(feature_list) & missing(feature_list)){stop("<feature_list> is not in the approriate character vector format or is missing")}

  data<-romics_object$data
  missingdata<-romics_object$missingdata
  if(!is.null(romics_object$statistics)){statistics<-romics_object$statistics}

  new_data<-data[-(1:nrow(data)),]
  new_missingdata<-missingdata[-(1:nrow(missingdata)),]
  if(!is.null(romics_object$statistics)){new_statistics<-statistics[-(1:nrow(statistics)),]}
  for(i in 1:length(feature_list)){
    datarow<-data[grepl(feature_list[i],rownames(data)),]
    new_data<-rbind(new_data,datarow)

    missingdatarow<-missingdata[grepl(feature_list[i],rownames(missingdata)),]
    new_missingdata<-rbind(new_missingdata,missingdatarow)

    if(!is.null(romics_object$statistics)){
      statisticsrow<-statistics[grepl(feature_list[i],rownames(statistics)),]
      new_statistics<-rbind(new_statistics,statisticsrow)
    }

    if(nrow(datarow)>1){
      print(paste("For the feature",feature_list[i],":",nrow(datarow),"elements of the data were matching, all matching elements were conserved."))}

  }
  romics_object$data<-new_data
  romics_object$missingdata<-new_missingdata
  if(!is.null(romics_object$statistics)){
    romics_object$statistics<-new_statistics
  }
  romics_object<-romicsUpdateSteps(romics_object,arguments)
  return(romics_object)
  }
