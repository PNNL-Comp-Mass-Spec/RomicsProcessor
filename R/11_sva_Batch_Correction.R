#' romicsBatchCorrection()
#' @description Performs batch correction on the data layer of the romics_object with options to use a reference region.
#' @param romics_object A romics_object created using romicsCreateObject()
#' @param batch_factor A character string specifying the factor for batch covariate.
#' @param method A character string that must be either 'ComBat', 'mean.only', or 'reference'. Default: "ComBat"
#' @param reference_filter A logical vector or function to identify the reference region. Required when method="reference".
#' @param verbose Logical. If TRUE, prints detailed information about the batch correction process. Default: TRUE
#' @param use_ref_params_only Logical. When TRUE, learns parameters only from reference region and applies to all data. Default: FALSE
#' @param ... Additional parameters passed to sva::ComBat().
#' @details
#' This function performs batch correction on a romics_object. When method="reference", it uses the reference
#' region to learn batch effect patterns and applies the same correction to all samples uniformly.
#'
#' For the reference-based method, you must provide a reference_filter that identifies the reference region.
#'
#' @return This function returns a transformed romics_object with batch-corrected data.
#' @author Modified by Brittney Gorman, with backward-compatible enhancements
#' @export
romicsBatchCorrection <- function(
  romics_object,
  batch_factor,
  method = "ComBat",
  reference_filter = NULL,
  verbose = TRUE,
  use_ref_params_only = FALSE,
  ...
) {
  # Record user input arguments
  arguments <- as.list(match.call())

  # Input validation
  if (!is.romicsObject(romics_object) || missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }

  if (missing(batch_factor)) {
    stop(
      "The <batch_factor> is missing. The list of romics_factor can be obtained with the function romicsFactorNames()"
    )
  }

  if (!batch_factor %in% rownames(romics_object$metadata)) {
    stop(
      "The <batch_factor> selected was not present in the romics_object. The list of romics_factor can be obtained with the function romicsFactorNames()"
    )
  }

  if (!method %in% c("ComBat", "mean.only", "reference")) {
    stop("<method> has to be either 'ComBat', 'mean.only', or 'reference'")
  }

  # Reference method validation
  if (method == "reference" && is.null(reference_filter)) {
    stop(
      "When method='reference', you must provide a reference_filter to identify the reference region"
    )
  }

  # Check for sva package
  if (!requireNamespace("sva", quietly = TRUE)) {
    stop(
      "This function requires the package 'sva' to be installed. Please install it from Bioconductor: BiocManager::install('sva')"
    )
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
    warning(
      "Only ",
      n_batches,
      " unique batch(es) detected. Batch correction may not be meaningful."
    )
  }

  # Check for missing values in data
  if (any(is.na(romics_object$data))) {
    stop(
      "The data layer contains missing values. ComBat requires complete data. Consider using imputation first."
    )
  }

  # Store original data for comparison
  original_data <- romics_object$data

  # Process reference region if needed
  if (method == "reference") {
    # Determine which samples are in the reference region
    if (is.logical(reference_filter)) {
      if (length(reference_filter) != ncol(romics_object$data)) {
        stop(
          "reference_filter logical vector length must match number of samples in romics_object"
        )
      }
      reference_samples <- reference_filter
    } else if (is.function(reference_filter)) {
      reference_samples <- reference_filter(romics_object)
      if (
        !is.logical(reference_samples) ||
          length(reference_samples) != ncol(romics_object$data)
      ) {
        stop(
          "reference_filter function must return a logical vector matching number of samples in romics_object"
        )
      }
    } else if (is.character(reference_filter)) {
      reference_samples <- colnames(romics_object$data) %in% reference_filter
    } else {
      stop(
        "reference_filter must be a logical vector, function, or character vector of sample names"
      )
    }

    # Validate reference samples
    if (sum(reference_samples) == 0) {
      stop("No reference samples identified from the reference_filter")
    }

    if (verbose) {
      message("Reference-based batch correction:")
      message("- Total samples: ", ncol(romics_object$data))
      message(
        "- Reference samples: ",
        sum(reference_samples),
        " (",
        round(sum(reference_samples) / ncol(romics_object$data) * 100, 1),
        "%)"
      )
    }

    # Get batch information for reference samples
    ref_batches <- batch_numeric[reference_samples]

    # Check if reference samples span all batches
    if (length(unique(ref_batches)) < n_batches) {
      warning(
        "Reference samples don't include all batches. Some batch effects may not be properly corrected."
      )
    }

    if (use_ref_params_only) {
      # Extract reference region data
      reference_data <- romics_object$data[, reference_samples, drop = FALSE]

      if (verbose) {
        message(
          "Learning batch effect parameters from reference region ONLY and applying to ALL data..."
        )
      }

      # Create a copy of the full data to manipulate
      full_data <- romics_object$data

      # Get the mean of each feature in the reference region for each batch
      batch_means <- list()
      batch_sds <- list()

      # For each batch, calculate the mean and SD of each feature in the reference region
      for (batch in unique(batch_numeric)) {
        batch_ref_samples <- reference_samples & (batch_numeric == batch)
        if (sum(batch_ref_samples) > 0) {
          batch_means[[as.character(batch)]] <- rowMeans(romics_object$data[,
            batch_ref_samples,
            drop = FALSE
          ])
          batch_sds[[as.character(batch)]] <- apply(
            romics_object$data[, batch_ref_samples, drop = FALSE],
            1,
            sd
          )
        } else {
          warning(
            "Batch ",
            batch,
            " has no reference samples. Using overall means."
          )
          batch_means[[as.character(batch)]] <- rowMeans(reference_data)
          batch_sds[[as.character(batch)]] <- apply(reference_data, 1, sd)
        }
      }

      # Calculate the grand mean of each feature across all reference samples
      grand_means <- rowMeans(reference_data)

      # Create a standardized version of all data
      standardized_data <- romics_object$data

      # For each batch, standardize the data
      for (batch in unique(batch_numeric)) {
        batch_samples <- batch_numeric == batch
        batch_mean <- batch_means[[as.character(batch)]]
        batch_sd <- batch_sds[[as.character(batch)]]

        # Standardize: (x - batch_mean) / batch_sd
        for (i in 1:nrow(standardized_data)) {
          if (batch_sd[i] > 0) {
            standardized_data[i, batch_samples] <-
              (romics_object$data[i, batch_samples] - batch_mean[i]) /
              batch_sd[i]
          } else {
            # If SD is zero, just center the data
            standardized_data[i, batch_samples] <-
              romics_object$data[i, batch_samples] - batch_mean[i]
          }
        }
      }

      # Rescale standardized data to the grand reference scale
      corrected_data <- standardized_data
      for (i in 1:nrow(corrected_data)) {
        # Calculate reference SD
        ref_sd <- sd(reference_data[i, ])
        if (is.na(ref_sd) || ref_sd == 0) {
          ref_sd <- 1
        }

        # Rescale: z-score * ref_sd + grand_mean
        corrected_data[i, ] <- standardized_data[i, ] * ref_sd + grand_means[i]
      }

      # Update romics object with the corrected data
      romics_object$data <- corrected_data
    } else {
      # Apply standard ComBat but prioritize the reference samples
      if (verbose) {
        message(
          "Running standard ComBat but using reference samples to guide the correction..."
        )
      }

      # Create design matrix that includes reference info
      mod <- matrix(1, ncol(romics_object$data), 1)
      # Add reference indicator as a covariate
      mod <- cbind(mod, reference_samples * 1) # Convert logical to numeric

      # Apply ComBat correction with reference covariate in the model
      corrected_data <- sva::ComBat(
        dat = as.matrix(romics_object$data),
        batch = batch_numeric,
        mod = mod, # Include reference status in the model
        mean.only = FALSE,
        ...
      )

      # Update romics object with corrected data
      romics_object$data <- corrected_data
    }
  } else {
    # Standard ComBat or mean.only method (backward compatible)
    if (verbose) {
      message("Running batch correction using sva::ComBat:")
      message("- Number of features: ", nrow(romics_object$data))
      message("- Number of samples: ", ncol(romics_object$data))
      message("- Number of batches: ", n_batches)
      message("- Method: ", method)
      message("- Mean only: ", method == "mean.only")
    }

    # Apply ComBat correction
    tryCatch(
      {
        corrected_data <- sva::ComBat(
          dat = as.matrix(romics_object$data),
          batch = batch_numeric,
          mean.only = (method == "mean.only"),
          ...
        )
        # Update romics object
        romics_object$data <- corrected_data
      },
      error = function(e) {
        stop("Error during ComBat correction: ", e$message)
      }
    )
  }

  # Print citation
  if (verbose) {
    message("Please cite the sva package:")
    print(citation("sva"), bibtex = FALSE)
  }

  # Verify that correction occurred
  if (verbose) {
    data_diff <- sum(abs(original_data - romics_object$data), na.rm = TRUE)
    message(
      "Sum of absolute differences between original and corrected data: ",
      round(data_diff, 6)
    )
    if (data_diff < 1e-10) {
      warning(
        "The corrected data appears identical to original data. Check your batch factor and ComBat parameters."
      )
    } else {
      message("Batch correction completed successfully.")
    }
  }

  # Update steps and dependencies
  romics_object <- romicsUpdateSteps(romics_object, arguments)
  romics_object <- romicsAddDependency(romics_object, new_dependency = "sva")

  return(romics_object)
}
