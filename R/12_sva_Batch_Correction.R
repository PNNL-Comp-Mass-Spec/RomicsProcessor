#This is a Wrapper for sva::combat batch correction
RomicsBatchCorrection<-function(romics_object, batch_factor= "factor" , method="method"){
  #general checking
  if(missing(romics_object)){stop("romics_object is missing")}
  if(class(romics_object)!="romics_object"){stop("your romics_object was not created using the function romicsCreateObject")}
  if(romics_object$steps[1]!="romics_object"){stop("romics_object is not in the appropriate format")}
  if(missing(batch_factor)){stop("Your batch factor is missing.")}
  if(missing(method)){method="ComBat"}
  if(method %in% c("ComBat","mean.only"))

  #seed setting
  #extract the data from the romics_object$data
  dat<-as.matrix(romics_object$data)

  #check if the batch_factor exists
  if(!batch_factor %in% rownames(romics_object$metadata)){
  print("Your batch factor has to be in the following list:")
  print(rownames(romics_object$metadata))
  stop()
  }

  #Extract the batch from the romics_object$metadata
  batch <- data.frame(t(romics_object$metadata))
  batch <- as.double(batch[,colnames(batch)==batch_factor])

  #if method = "ComBat" run the ComBat
  if(method == "ComBat"){
  corrected_dat<-sva::ComBat(dat, batch)
  }

  #if method="mean.only" run ComBat with mean.only
  if(method == "mean.only"){
    corrected_dat<-sva::ComBat(dat, batch, mean.only = TRUE)
  }

  #update the data
  romics_object$data<-corrected_dat

  #update the steps
  romics_object$steps<- c(romics_object$steps,paste0(gsub(" ","_",format(Sys.time(),"%b_%d_%Y_%X")),": the data was batch corrected"))
  romics_object$steps<- c(romics_object$steps,paste0('RomicsBatchCorrection(romics_object,batch_factor="',batch_factor,'",method="',method,'",',sep=""))

  return(romics_object)
}
