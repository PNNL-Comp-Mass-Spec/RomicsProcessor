#' romicsBatchCorrection()
#' @description Performs the sva::ComBat() batch correction on the data layer of the romics_object. The data layer must not contain missing values and the factor utilized will be the one used for the correction.
#' @param romics_object A romics_object created using romicsCreateObject()
#' @param batch_factor A character string specifying the factor contained in the romics_object that will serve as batch covariate. To obtain the list of factors, use romicsFactorNames()
#' @param method A character string that must be either 'ComBat' or 'mean.only' to indicate how the ComBat function should be run.
#' @param verbose Logical. If TRUE, prints detailed information about the batch correction process. Default: TRUE
#' @param ... Additional parameters passed to sva::ComBat(). See sva::ComBat() documentation for more details.
#' @details This function performs ComBat batch correction on a romics_object. It can be performed using the ComBat method or using a mean.only method. The batch factor is automatically converted to the appropriate format for ComBat. See sva::ComBat() documentation for more details.
#' @return This function returns a transformed romics_object with batch-corrected data.
#' @author Geremy Clair
#' @export
romicsBatchCorrection <- function(romics_object, batch_factor, method = "ComBat", verbose = TRUE, ...) {

  # Record user input arguments
  arguments <- as.list(match.call())

  # Input validation
  if (!is.romicsObject(romics_object) || missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }

  if (missing(batch_factor)) {
    stop("The <batch_factor> is missing. The list of romics_factor can be obtained with the function romicsFactorNames()")
  }

  if (!batch_factor %in% rownames(romics_object$metadata)) {
    stop("The <batch_factor> selected was not present in the romics_object. The list of romics_factor can be obtained with the function romicsFactorNames()")
  }

  if (!method %in% c("ComBat", "mean.only")) {
    stop("<method> has to be either 'ComBat' or 'mean.only'")
  }

  # Check for sva package
  if (!requireNamespace("sva", quietly = TRUE)) {
    stop("This function requires the package 'sva' to be installed. Please install it from Bioconductor: BiocManager::install('sva')")
  }

  # Extract and process batch factor
  batch_raw <- romicsExtractFactor(romics_object, batch_factor)
  if (verbose) {
    message("Original batch factor:")
    print(table(batch_raw))
  }

  # Convert batch factor to appropriate format for ComBat
  if (is.character(batch_raw) || is.factor(batch_raw)) {
    # Convert to factor first, then to numeric
    batch_factor_levels <- as.factor(batch_raw)
    batch_numeric <- as.numeric(batch_factor_levels)

    if (verbose) {
      message("Batch factor after conversion:")
      batch_mapping <- data.frame(
        Original = levels(batch_factor_levels),
        Numeric = 1:nlevels(batch_factor_levels)
      )
      print(batch_mapping)
      print(table(batch_numeric))
    }
  } else if (is.numeric(batch_raw)) {
    batch_numeric <- batch_raw
    if (verbose) {
      message("Batch factor is already numeric:")
      print(table(batch_numeric))
    }
  } else {
    stop("Unable to convert batch factor to appropriate format for ComBat")
  }

  # Check if there are multiple batches
  n_batches <- length(unique(batch_numeric[!is.na(batch_numeric)]))
  if (n_batches < 2) {
    warning("Only ", n_batches, " unique batch(es) detected. Batch correction may not be meaningful.")
  }

  # Check for missing values in data
  if (any(is.na(romics_object$data))) {
    stop("The data layer contains missing values. ComBat requires complete data. Consider using imputation first.")
  }

  # Set mean.only parameter
  mean_only <- (method == "mean.only")

  if (verbose) {
    message("Running batch correction using sva::ComBat:")
    message("- Number of features: ", nrow(romics_object$data))
    message("- Number of samples: ", ncol(romics_object$data))
    message("- Number of batches: ", n_batches)
    message("- Method: ", method)
    message("- Mean only: ", mean_only)
  }

  # Print citation
  if (verbose) {
    message("Please cite the sva package:")
    print(citation("sva"), bibtex = FALSE)
  }

  # Store original data for comparison
  original_data <- romics_object$data

  # Apply ComBat correction
  tryCatch({
    corrected_data <- sva::ComBat(
      dat = as.matrix(romics_object$data),
      batch = batch_numeric,
      mean.only = mean_only,
      ...
    )

    # Update romics object
    romics_object$data <- corrected_data

    # Verify that correction occurred
    if (verbose) {
      data_diff <- sum(abs(original_data - corrected_data), na.rm = TRUE)
      message("Sum of absolute differences between original and corrected data: ", round(data_diff, 6))

      if (data_diff < 1e-10) {
        warning("The corrected data appears identical to original data. Check your batch factor and ComBat parameters.")
      } else {
        message("Batch correction completed successfully.")
      }
    }

  }, error = function(e) {
    stop("Error during ComBat correction: ", e$message)
  })

  # Update steps and dependencies
  romics_object <- romicsUpdateSteps(romics_object, arguments)
  romics_object <- romicsAddDependency(romics_object, new_dependency = "sva")

  return(romics_object)
}
