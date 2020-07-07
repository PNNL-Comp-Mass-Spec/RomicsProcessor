#' romicsBatchCorrection()
#' @description Performs the sva::ComBat() batch correction on the data layer of the romics_object. The data layer must not contain missing values and the factor utilized will be the one used for the correction.
#' @param romics_object A romics_object created using romicsCreateObject()
#' @param batch_factor has to be a factor contained in the romics_object that will serve as batch covariate To obtain the
#' @param method has to be either 'ComBat' or 'mean.only' to indicate how the ComBat function should be run.
#' @param ... parameters can be passed to sva::ComBat(), see sva::ComBat() documentation for more details.
#' @details This function is used to perform a ComBat batch correction on a romics_object. it can be performed using the ComBat method or using a mean.only method. sva::ComBat() documentation for more details.
#' @return This function returns a transformed romics_object.
#' @author Geremy Clair
#' @export
romicsBatchCorrection<-function(romics_object, batch_factor="factor", method="ComBat", ...){
  #run the line to record the user input arguments
  arguments<-as.list(match.call())

  #general checking
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(batch_factor)){stop("The <batch_factor> is missing the list of romics_factor can be obtained with the function romicsFactorNames()")}
  if(!batch_factor %in% rownames(romics_object$metadata)){stop("The <batch_factor> selected was not present in the romics_object. The list of romics_factor can be obtained with the function romicsFactorNames()")}
  if(missing(method)){method="ComBat"}
  if(!method %in% c("ComBat","mean.only")){stop("<method> has to be either 'ComBat' or 'mean.only'")}


  #verify if the package sva is installed if not indicate to the user that sva is required to run this function (This step would not be necessary if the romicsBatchCorrection was part of the sva package).
  if(!"sva" %in% rownames(installed.packages())){stop("This function requires the package 'sva' to be installed to work. Please install this package from bioConductor prior to run this function.")}
  #verify if the sva library is loaded if not load it and indicate to the user that the library was loaded (This step would not be necessary if the romicsBatchCorrection was part of the sva package).
  if(!"sva" %in% (.packages())){
    print("The package 'sva' was not loaded, it was loaded to execute this function.")
    library("sva")
    }

  #The batch_factor is extracted using the function romicsExtractFactor()
  f <- romicsExtractFactor(romics_object,batch_factor)
  #Sva::ComBat require the batch to be numerical the factor extracted is converted in double (see ComBat documentation for more details)
  f <- as.double(f)

  #Set the content of the mean.only based in the content of method
  if(method=="ComBat"){mean.only=FALSE}else{mean.only=TRUE}

  #print a message indicating that sva::ComBat() is applied
  print("A batch correction using the sva::ComBat method will be applied to the romics_object.")
  print(citation("sva"),bibtex =F)

  #Run the sva::ComBat on the data of the romics_object,
  romics_object$data<- sva::ComBat(as.matrix(romics_object$data),f,mean.only = mean.only,...)

  #For romics function it is necessary to update the steps
  romics_object<-romicsUpdateSteps(romics_object,arguments)

  #As sva is required for this function to run and in order to keep track of dependency used for the romics_object to be transformed
  romics_object<-romicsAddDependency(romics_object,new_dependency = "sva")

  return(romics_object)
}
