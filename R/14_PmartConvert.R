#' RomicsPmartR()
#' Converts romics_object to pmartR object.
#'
#' @param romics_object A romics_object created using the function romicsCreateObject()
#' @param type Has to be "lipidData","proData","pepData","metabData" to indicate what data type to use for the pmartR object
#'
#' @details This function converts an romics_object to a pmartR object
#' @return return the pmartR object
#'
#' @author Geremy Clair
#' @export
#'
RomicsPmartR<- function(romics_object, type="proData"){
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if (missing(type)){type="proData"}
  if (!type %in% c("lipidData","proData","pepData","metabData")){stop("your data should be one of the following: 'lipidData','proData','pepData','metabData'")}

e_data <- data.frame(cbind(Mass_Tag_ID = 1:nrow(romics_object$data),romics_object$data))
rownames(e_data)  <- NULL
f_data <- data.frame(cbind(SampleID=colnames(romics_object$metadata),t(romics_object$metadata)))
rownames(f_data)  <- NULL
e_meta <- data.frame(cbind(Mass_Tag_ID = e_data$Mass_Tag_ID,IDs=rownames(romics_object$data)))
if(sum(grepl("log10transform",romics_object$steps))+sum(grepl("log2transform",romics_object$steps))>0){data_scale = "log"} else {data_scale = "abundance"}
if(type=="pepData"){
  output <- as.pepData(e_data = e_data, f_data = f_data, e_meta = e_meta, edata_cname = "Mass_Tag_ID", fdata_cname = "SampleID", emeta_cname = "Mass_Tag_ID", data_scale = data_scale, check.names = FALSE)
  return(pepData)
  }
if(type=="proData"){
  output <- as.proData(e_data = e_data, f_data = f_data, e_meta = e_meta, edata_cname = "Mass_Tag_ID", fdata_cname = "SampleID", emeta_cname = "Mass_Tag_ID", data_scale = data_scale, check.names = FALSE)
  }
if(type=="lipidData"){
  output <- as.lipidData(e_data = e_data, f_data = f_data, e_meta = e_meta, edata_cname = "Mass_Tag_ID", fdata_cname = "SampleID", emeta_cname = "Mass_Tag_ID", data_scale = data_scale, check.names = FALSE)
}
if(type=="metabData"){
  output <- as.metabData(e_data = e_data, f_data = f_data, e_meta = e_meta, edata_cname = "Mass_Tag_ID", fdata_cname = "SampleID", emeta_cname = "Mass_Tag_ID", data_scale = data_scale, check.names = FALSE)
}
print(paste("An object of class", class(output), "was created", sep=" "))
return(output)
}
