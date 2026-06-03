#' romicsPmartR()
#' @description DEPRECATED: Use romicsToPmartR() instead (moved to 24_Format_Conversion.R)
#' Converts romics_object to pmartR object (requires the "pmartR' package to be installed).
#' @param romics_object A romics_object created using the function romicsCreateObject()
#' @param type Has to be "lipidData","proData","pepData","metabData" to indicate what data type to use for the pmartR object
#' @details This function is deprecated. Use romicsToPmartR() instead, which has been moved to the Format_Conversion module.
#' This wrapper maintains backward compatibility by calling the new function.
#' @return return the pmartR object
#' @author Geremy Clair
#' @export
romicsPmartR <- function(romics_object, type = "proData") {
  warning("romicsPmartR() is deprecated. Use romicsToPmartR() instead.")
  # Call the new function with renamed parameter
  romicsToPmartR(romics_object, data_type = type)
}
