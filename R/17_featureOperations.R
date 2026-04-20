#' romicsFeatureSubset()
#' @description This will create an romics_object subseted for a list of features only
#' @param romics_object has to be an romics_object created using romicsCreateObject()
#' @param feature_list has to be a character vector containing features to be extracted
#' @details This function will create a new romics object in which only the features contained in the feature_list are contained.
#' @details The function supports exact matches only, if some features are missing, the function will still work but only will return an object with the existing features.
#' @return This function returns a subseted romics_object
#' @author Geremy Clair
#' @export
romicsFeatureSubset<-function(romics_object,feature_list=c("feature1","feature2","etc")){
  arguments<-as.list(match.call())
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("<romics_object> is missing or is not in the appropriate format")}
  if(!is.character(feature_list) & missing(feature_list)){stop("<feature_list> is not in the approriate character vector format or is missing")}

  if(sum(!feature_list %in% rownames(romics_object$data)>0)){
    warning("Some of the feature in your <feature_list> were not present in the romics_object and were therefore not included:")
    print(feature_list[!feature_list %in% rownames(romics_object$data)>0])
  }

  feat_name<-romicsFeatureNames(romics_object)
  romics_object$data<-romics_object$data[feat_name %in% feature_list,]
  romics_object$missingdata<-romics_object$missingdata[feat_name %in% feature_list,]
  if(!is.null(romics_object$statistics)){
  romics_object$statistics<-romics_object$statistics[feat_name %in% feature_list,]}

  romics_object<-romicsUpdateSteps(romics_object,arguments)
  return(romics_object)
  }

#' romicsRemoveUnknown()
#' @description Remove list of unkown features from the dataset
#' @param romics_object has to be an romics_object created using romicsCreateObject()
#' @param unknown_identifier_txt string that indicates which features are unknown and are to be removed by default the string 'Unknown'
#' @details This function will create a new romics object in which only the unkown features are removed from the data, missingdata and statistics layers.
#' @details The function removes unkown features
#' @return This function returns an romics_object
#' @author Geremy Clair
#' @export
romicsRemoveUnknown<-function(romics_object,unknown_identifier_txt="Unknown"){
  arguments<-as.list(match.call())
  if(missing(unknown_identifier_txt)){unknown_identifier_txt="Unknown"}
  if(!is.character(unknown_identifier_txt)&length(unknown_identifier_txt)!=1){stop("Your <unknown_identifier_txt> has to be a string of lenght 1.")}
  if(!is.romicsObject(romics_object)){stop("Your <romics_object> was not in the appropriate format.")}
  romics_object$data<-romics_object$data[!grepl(unknown_identifier_txt,rownames(romics_object$data)),]
  romics_object$missingdata<-romics_object$missingdata[!grepl(unknown_identifier_txt,rownames(romics_object$missingdata)),]
  romics_object$statistics<-romics_object$statistics[!grepl(unknown_identifier_txt,rownames(romics_object$statistics)),]
  romics_object<-romicsUpdateSteps(romics_object,arguments)
}

#' romicsFeatureCutoff()
#' @description Applies cutoff thresholds to specified features in a romics_object. Values outside the specified range are set to NA and the missingdata layer is updated accordingly.
#' @param romics_object A romics_object created using romicsCreateObject().
#' @param feature_name Character string or vector. Name(s) of the feature(s) to apply cutoffs to. Must match row names in the data exactly.
#' @param low_cutoff Numeric. Lower threshold value. Values below this will be set to NA. If NULL, uses the minimum value of the feature. Default: NULL
#' @param high_cutoff Numeric. Upper threshold value. Values above this will be set to NA. If NULL, uses the maximum value of the feature. Default: NULL
#' @param apply_to_all Logical. If TRUE, applies the same cutoffs to all features in the dataset. Overrides feature_name. Default: FALSE
#' @details This function is useful for removing outliers or applying quality control thresholds to specific features.
#' When cutoffs are applied, the original data structure is preserved but values outside the range become missing values.
#' The missingdata layer is automatically updated to reflect the new missing values.
#' If low_cutoff or high_cutoff are NULL, they default to the min/max values of each feature respectively,
#' effectively meaning no cutoff is applied on that end.
#' @return Returns a romics_object with cutoffs applied and missingdata layer updated.
#' @author Geremy Clair
#' @export
romicsFeatureCutoff <- function(romics_object,
                                feature_name = NULL,
                                low_cutoff = NULL,
                                high_cutoff = NULL,
                                apply_to_all = FALSE) {

  arguments <- as.list(match.call())

  # Input validation
  if(!is.romicsObject(romics_object) || missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }

  if(!apply_to_all && (is.null(feature_name) || missing(feature_name))) {
    stop("feature_name must be provided when apply_to_all is FALSE")
  }

  if(!is.logical(apply_to_all)) {
    stop("apply_to_all must be TRUE or FALSE")
  }

  # Validate cutoff values if provided
  if(!is.null(low_cutoff) && (!is.numeric(low_cutoff) || length(low_cutoff) != 1)) {
    stop("low_cutoff must be a single numeric value or NULL")
  }

  if(!is.null(high_cutoff) && (!is.numeric(high_cutoff) || length(high_cutoff) != 1)) {
    stop("high_cutoff must be a single numeric value or NULL")
  }

  if(!is.null(low_cutoff) && !is.null(high_cutoff) && low_cutoff >= high_cutoff) {
    stop("low_cutoff must be less than high_cutoff")
  }

  # Determine which features to process
  if(apply_to_all) {
    features_to_process <- rownames(romics_object$data)
    message("Applying cutoffs to all ", length(features_to_process), " features")
  } else {
    # Validate feature names exist
    missing_features <- feature_name[!feature_name %in% rownames(romics_object$data)]
    if(length(missing_features) > 0) {
      stop(paste("Feature(s) not found in data:", paste(missing_features, collapse = ", ")))
    }
    features_to_process <- feature_name
    message("Applying cutoffs to ", length(features_to_process), " specified feature(s)")
  }

  # Initialize counters for reporting
  total_values_processed <- 0
  total_values_set_to_na <- 0

  # Process each feature
  for(feature in features_to_process) {
    # Extract feature data
    feature_data <- as.numeric(romics_object$data[feature, ])
    original_na_count <- sum(is.na(feature_data))

    # Determine actual cutoffs for this feature
    if(is.null(low_cutoff)) {
      actual_low <- min(feature_data, na.rm = TRUE)
    } else {
      actual_low <- low_cutoff
    }

    if(is.null(high_cutoff)) {
      actual_high <- max(feature_data, na.rm = TRUE)
    } else {
      actual_high <- high_cutoff
    }

    # Count values that will be affected
    values_below_low <- sum(feature_data < actual_low, na.rm = TRUE)
    values_above_high <- sum(feature_data > actual_high, na.rm = TRUE)
    values_affected <- values_below_low + values_above_high

    # Apply cutoffs
    feature_data[feature_data < actual_low | feature_data > actual_high] <- NA

    # Update the data
    romics_object$data[feature, ] <- feature_data

    # Update counters
    total_values_processed <- total_values_processed + length(feature_data)
    total_values_set_to_na <- total_values_set_to_na + values_affected

    # Feature-specific reporting for individual features
    if(!apply_to_all) {
      final_na_count <- sum(is.na(feature_data))
      new_na_count <- final_na_count - original_na_count

      message(paste("Feature:", feature))
      message(paste("  Cutoff range: [", actual_low, ", ", actual_high, "]"))
      message(paste("  Values set to NA:", new_na_count,
                    "(", values_below_low, " below low +", values_above_high, " above high)"))
      message(paste("  Total NA values now:", final_na_count, "/", length(feature_data)))
    }
  }

  # Update missingdata layer
  romics_object$missingdata <- is.na(romics_object$data)

  # Summary reporting
  message("Cutoff summary:")
  message(paste("  Total values processed:", total_values_processed))
  message(paste("  Total values set to NA:", total_values_set_to_na,
                "(", round(total_values_set_to_na/total_values_processed*100, 2), "%)"))

  # Update steps
  romics_object <- romicsUpdateSteps(romics_object, arguments)

  return(romics_object)
}
