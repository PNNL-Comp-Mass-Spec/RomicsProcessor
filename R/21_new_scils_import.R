#' Extract Data from SCiLS Lab Files
#' @description Extracts intensity data, metadata, and feature IDs from a SCiLS '.slx' file.
#' The function efficiently processes large datasets common in mass spectrometry imaging.
#'
#' @param file Character string. Path to the SCiLS '.slx' file.
#' @param scils_executable Character string. Path to the SCiLS executable.
#'        Default: 'C:/Program Files/SCiLS/SCiLS Lab/scilsMsServer.exe'
#' @param port Numeric. Port to use for the connection with the SCiLS environment.
#'        Default: 8082
#' @param feature_list Character string. Name of the feature list to import.
#'        Default: "All Features"
#' @param normId Character string. Normalization ID to apply during import.
#'        Default: "" (no normalization)
#' @param metadata_fields Character vector. Names of metadata fields to extract.
#'        Default: NULL (extracts all available metadata fields)
#'
#' @details
#' This function extracts data from a SCiLS '.slx' file through the SCiLS Lab API.
#' It requires both the '.slx' file and corresponding '.sbd' file to be present in the same directory.
#'
#' The function is optimized for handling large datasets with potentially hundreds of thousands of
#' data points. It imports:
#' \itemize{
#'   \item Intensity matrix for the selected features
#'   \item Spot IDs and associated metadata
#'   \item Detailed information about each feature
#' }
#'
#' The necessary companion '.sbd' file must be in the same directory as the '.slx' file.
#'
#' @return A list with three components:
#' \describe{
#'   \item{data}{Data frame containing the intensity matrix with features as rows and samples as columns}
#'   \item{metadata}{Data frame containing spot IDs and other spot attributes}
#'   \item{IDs}{Data frame containing detailed information about each feature}
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage with default parameters
#' scils_data <- extractScils("path/to/your/file.slx")
#'
#' # Custom parameters
#' scils_data <- extractScils(
#'   file = "path/to/your/file.slx",
#'   scils_executable = "C:/Custom/Path/scilsMsServer.exe",
#'   port = 8083,
#'   feature_list = "My Custom Feature List",
#'   normId = "TIC"
#' )
#'
#' # Access the components
#' data <- scils_data$data
#' metadata <- scils_data$metadata
#' feature_ids <- scils_data$IDs
#' }
#'
#' @author Geremy Clair, Brittney Gorman
#' @export
extractScils <- function(file, scils_executable, port, feature_list, normId, metadata_fields = NULL) {
  # Input validation with detailed error messages
  if(missing(file)) {
    stop("Parameter 'file' is required: please provide the path to your SCiLS .slx file.")
  }

  if(!is.character(file) || length(file) != 1) {
    stop("Parameter 'file' must be a character string specifying the path to your SCiLS .slx file.")
  }

  if(!file.exists(file)) {
    stop(sprintf("File not found: '%s'. Please verify the file path is correct.", file))
  }

  # Check for the basedata file
  basedata_file <- gsub("\\.slx$", ".sbd", file)
  if(!file.exists(basedata_file)) {
    stop(sprintf("Required companion file not found: '%s'.\n  This file must be in the same directory as your .slx file.", basedata_file))
  }

  # Set default parameters if missing
  if(missing(scils_executable)) {
    scils_executable <- "C:/Program Files/SCiLS/SCiLS Lab/scilsMsServer.exe"
  }

  if(!file.exists(scils_executable)) {
    stop(sprintf("SCiLS executable not found: '%s'.\n  Please provide the correct path to scilsMsServer.exe.", scils_executable))
  }

  if(missing(port)) {
    port <- 8082
  } else if(!is.numeric(port) || length(port) != 1) {
    stop("Parameter 'port' must be a single numeric value.")
  }

  if(missing(feature_list)) {
    feature_list <- "All Features"
  }

  if(missing(normId)) {
    normId <- ""
  }

  # Ensure proper cleanup on exit or error
  ScilsFileEnv <- NULL
  on.exit({
    if(!is.null(ScilsFileEnv)) {
      message(sprintf("Closing SCiLS Lab session (%s).", Sys.time()))
      try(close(ScilsFileEnv), silent = TRUE)
    }
  })

  # Load the SCiLS environment and import the data
  tryCatch({
    ScilsFileEnv <- SCiLSLabOpenLocalSession(file, port = port, executable = scils_executable)
    message(sprintf("%s was opened to import the SCiLS data in R (%s).",
                    getServerVersion(ScilsFileEnv), Sys.time()))
  }, error = function(e) {
    stop(sprintf("Failed to open SCiLS session: %s", e$message))
  })

  # Load available feature lists
  tryCatch({
    Flists <- as.character(t(SCiLSLabClient::getFeatureLists(ScilsFileEnv)[, 2]))

    # Verify if the feature_list exists
    if(!feature_list %in% Flists) {
      stop(sprintf("Feature list '%s' not found. Available lists are: %s",
                   feature_list, paste(Flists, collapse = ", ")))
    }

    # Import the list of features
    selected_features <- SCiLSLabClient::getFeatures(
      ScilsFileEnv,
      name = feature_list,
      includeVisibleUserColumns = TRUE
    )
  }, error = function(e) {
    stop(sprintf("Failed to retrieve feature lists: %s", e$message))
  })

  # Import intensity matrix
  message(sprintf("Beginning data matrix importation (%s).", Sys.time()))
  tryCatch({
    data <- do.call(rbind, lapply(selected_features$id, function(x) {
      intensities <- SCiLSLabClient::getFeatureIntensities(
        ScilsFileEnv,
        x,
        regionId = 'Regions',
        normId = normId
      )$intensity
      as.numeric(intensities)  # Convert to numeric
    }))
    message(sprintf("Data matrix imported (%s).", Sys.time()))
  }, error = function(e) {
    stop(sprintf("Failed to import data matrix: %s", e$message))
  })

  # Get spot IDs
  spot_ids <- getRegionSpots(ScilsFileEnv, regionId = 'Regions')$spotId

  # Import metadata using the more efficient approach
  message(sprintf("Beginning metadata importation (%s).", Sys.time()))
  tryCatch({
    # Get region tree for metadata extraction
    regionTree <- getRegionTree(ScilsFileEnv)

    # Ensure regionTree is of the correct class
    if (!"RegionTree" %in% class(regionTree)) {
      stop("RegionTree argument is not of class 'RegionTree'")
    }

    # Flatten the region tree
    allRegions <- flattenRegionTree(regionTree)

    # Initialize spots data frame with coordinates
    spots_data <- allRegions[[1]]$spots
    base_fields <- c("spotId", "x", "y", "z")
    base_fields <- base_fields[base_fields %in% colnames(spots_data)]

    # Extract only base fields initially
    return_df <- spots_data[, base_fields, drop = FALSE]

    # Build attribute dictionary
    attributes <- list()
    for (region in allRegions) {
      # Skip regions with no attributes or only single attributes
      if (is.null(region$attributes) ||
          (is.character(region$attributes$name) && length(region$attributes$name) <= 1)) {
        next
      }

      # Process each attribute in the region
      for (att in seq_len(nrow(region$attributes))) {
        attribute_name <- region$attributes[att, ]$name

        # Skip "Date" attributes
        if (attribute_name == "Date") {
          next
        }

        attribute_level <- region$attributes[att, ]$value

        # Handle "Class" attribute specially
        if (attribute_name == "Class") {
          attribute_name <- attribute_level
          attribute_level <- "positive"
        }

        # Skip if metadata_fields is specified and this attribute isn't in it
        if (!is.null(metadata_fields) && !(attribute_name %in% metadata_fields)) {
          next
        }

        # Initialize attribute_name and levels if not already in the list
        if (!attribute_name %in% names(attributes)) {
          attributes[[attribute_name]] <- list()
        }

        if (!attribute_level %in% names(attributes[[attribute_name]])) {
          attributes[[attribute_name]][[attribute_level]] <- region$spots$spotId
        } else {
          # Append unique values using union
          attributes[[attribute_name]][[attribute_level]] <- union(
            attributes[[attribute_name]][[attribute_level]], region$spots$spotId
          )
        }
      }
    }

    # If metadata_fields is NULL, use all attributes found
    if (is.null(metadata_fields)) {
      metadata_fields <- names(attributes)
    } else {
      # Ensure we only process fields that actually exist
      metadata_fields <- intersect(metadata_fields, names(attributes))
    }

    # Vectorized function to determine attribute levels for a spot
    get_region_levels <- function(spot_id, attr) {
      levels <- names(attr)[sapply(attr, function(level_spots) {
        spot_id %in% level_spots
      })]
      if (length(levels) == 0) {
        return("unknown")
      } else {
        return(paste(levels, collapse = ","))
      }
    }

    # Process each metadata field
    for (attribute_name in metadata_fields) {
      if (!attribute_name %in% names(attributes)) {
        # Add column of "unknown" if attribute doesn't exist
        return_df[[attribute_name]] <- "unknown"
        next
      }

      attr <- attributes[[attribute_name]]

      # Use vectorized apply for efficiency
      return_df[[attribute_name]] <- vapply(
        return_df$spotId,
        FUN = get_region_levels,
        FUN.VALUE = character(1),
        attr = attr
      )
    }

    # Clean up any NA values
    return_df[is.na(return_df)] <- "unknown"

    # Clean up duplicate positive values
    for (col in colnames(return_df)) {
      return_df[[col]] <- gsub("positive,positive(,positive)*", "positive", return_df[[col]])
    }

    # Format for consistency with the original function
    metadata <- t(return_df)
    rownames(metadata) <- c(base_fields, metadata_fields)

    message(sprintf("Metadata imported (%s).", Sys.time()))
  }, error = function(e) {
    stop(sprintf("Failed to import metadata: %s", e$message))
  })



  # Process feature IDs
  message(sprintf("Beginning of IDs importation (%s).", Sys.time()))
  tryCatch({
    # Create feature names
    features <- makeFeatureNames(selected_features, useName = FALSE, digits = 5)

    # Handle the "m/z:NA" case
    na_rows <- rownames(features) == "m/z:NA"
    if (any(na_rows)) {
      features[na_rows, ] <- makeFeatureNames(selected_features[na_rows, ], useName = TRUE)
    }

    # Drop duplicated rows
    duplicated_mask <- duplicated(features) | (features == "m/z:NA")
    features <- features[!duplicated_mask]
    data <- data[!duplicated_mask, ]

    # Combine features with intensity matrix and set column names
    spot_id_cols <- paste0("spotId_", spot_ids)
    colnames(data) <- spot_id_cols
    data <- as.data.frame(cbind(features, data))

    # Process metadata to match data format
    colnames(metadata) <- spot_id_cols
    metadata <- cbind(data.frame(factor = rownames(metadata)), metadata)

    # Process IDs with deduplication
    IDs <- selected_features[!duplicated_mask, ]
    IDs <- cbind(features, IDs)
    IDs <- IDs[!duplicated(IDs$features), ]
    rownames(IDs) <- seq_len(nrow(IDs))
  }, error = function(e) {
    stop(sprintf("Failed to process feature IDs: %s", e$message))
  })

  # Create the result list
  result <- list(
    data = data,
    metadata = metadata,
    IDs = IDs
  )
  close(ScilsFileEnv)
  message(sprintf("SCiLS data has been successfully extracted (%s).", Sys.time()))

  return(result)
}
