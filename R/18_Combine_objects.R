#' combineRomicsObjects()
#' @description Combine 2 or more Romics_objects, containing the same <Samples>. It will remove the calculated stats ensure the log transformation are identical. For the combinations of multiple spatial objects use the function combineRomicsSpatialObjects()
#' @param ... A succession of romics_objects the object should have the exact same <sample> names (data column names).
#' @details This function will combine 2 or more Romics_objects, containing the same <Samples>. It will remove the calculated stats ensure the log transformation are identical. For the combinations of multiple spatial objects use the function combineRomicsSpatialObjects()
#' @return A combined romics_object containing the normalized data from multiple romics_objects
#' @author Geremy Clair
#' @export
combineRomicsObjects <- function(...) {
  romics_objects <- list(...)
  names(romics_objects) <- as.character(match.call())[-1]

  # Input validation
  for (i in 1:length(romics_objects)) {
    if (!is.romicsObject(romics_objects[[i]])) {
      stop(paste0("The object named <", names(romics_objects)[i], "> is not an romics_object."))
    }
  }

  # Initialize data frames
  og_data <- data.frame()
  norm_data <- data.frame()
  missing_data <- data.frame()
  combined_IDs <- data.frame()
  combined_original_metadata <- list()  # New: for combining original metadata

  # Get common properties from first object
  samples <- romicsSampleNames(romics_objects[[1]])
  f <- romics_objects[[1]]$main_factor
  log <- numeric()
  m <- cbind(f = romicsFactorNames(romics_objects[[1]]), romics_objects[[1]]$metadata)
  rownames(m) <- NULL

  # Check log transformation status
  for (i in 1:length(romics_objects)) {
    if (romicsLogCheck(romics_objects[[i]])) {
      if (sum(grepl("log2transform\\(", romics_objects[[i]]$steps)) -
          sum(grepl("unlog2data\\(", romics_objects[[i]]$steps)) == 1) {
        log <- c(log, 2)
      } else {
        log <- c(log, 10)
      }
    } else {
      log <- c(log, FALSE)
    }
  }
  names(log) <- names(romics_objects)

  # Validate log status consistency
  if (!sum(log == log[1]) == length(romics_objects)) {
    print("The log status of the different objects were the following:")
    print(log)
    stop("The log status was different in the different romics_object. Please ensure the different objects are logged identically.")
  }

  # Find common columns across all objects for each data type
  # For original_data
  all_original_cols <- list()
  for (i in 1:length(romics_objects)) {
    if (!is.null(romics_objects[[i]]$original_data)) {
      all_original_cols[[i]] <- colnames(romics_objects[[i]]$original_data)
    }
  }
  common_original_cols <- if (length(all_original_cols) > 0) {
    Reduce(intersect, all_original_cols)
  } else {
    NULL
  }

  # For missingdata
  all_missing_cols <- list()
  for (i in 1:length(romics_objects)) {
    if (!is.null(romics_objects[[i]]$missingdata)) {
      all_missing_cols[[i]] <- colnames(romics_objects[[i]]$missingdata)
    }
  }
  common_missing_cols <- if (length(all_missing_cols) > 0) {
    Reduce(intersect, all_missing_cols)
  } else {
    NULL
  }

  # Process each romics object
  for (i in 1:length(romics_objects)) {
    # Validate sample names consistency
    if (!sum(romicsSampleNames(romics_objects[[i]]) %in% samples) == length(samples)) {
      stop(paste0("The samples present in the romics_object are different; these cannot be combined."))
    }

    # Remove statistics layer if present
    if ("statistics" %in% names(romics_objects[[i]])) {
      message(paste0("The romics_object named <", names(romics_objects)[i],
                     "> had a statistics layer, it was removed."))
    }

    # Change factor to match main factor
    romics_objects[[i]] <- romicsChangeFactor(romics_objects[[i]], main_factor = f)

    # Combine original data (only common columns)
    if (!is.null(romics_objects[[i]]$original_data) && !is.null(common_original_cols)) {
      current_og_data <- romics_objects[[i]]$original_data[, common_original_cols, drop = FALSE]
      if (nrow(og_data) == 0) {
        og_data <- current_og_data
      } else {
        og_data <- rbind(og_data, current_og_data)
      }
    }

    # Process normalized data with mode suffix
    nd <- romics_objects[[i]]$data
    mode_suffix <- as.character(names(romics_objects)[i])
    nd <- cbind(IDs = paste0(rownames(nd), "@", mode_suffix), nd)
    rownames(nd) <- NULL

    if (nrow(norm_data) == 0) {
      norm_data <- nd
    } else {
      norm_data <- rbind(norm_data, nd)
    }

    # Combine missing data (only common columns)
    if (!is.null(romics_objects[[i]]$missingdata) && !is.null(common_missing_cols)) {
      current_missing_data <- romics_objects[[i]]$missingdata[, common_missing_cols, drop = FALSE]
      if (nrow(missing_data) == 0) {
        missing_data <- current_missing_data
      } else {
        missing_data <- rbind(missing_data, current_missing_data)
      }
    }

    # Store original metadata for each object
    if (!is.null(romics_objects[[i]]$metadata)) {
      combined_original_metadata[[names(romics_objects)[i]]] <- romics_objects[[i]]$metadata
    }

    # Process IDs layer if it exists
    if ("IDs" %in% names(romics_objects[[i]])) {
      current_IDs <- romics_objects[[i]]$IDs

      # Add mode suffix to rownames to match the data layer format
      modified_rownames <- paste0(rownames(current_IDs), "@", mode_suffix)

      # Create a data frame with modified IDs
      current_IDs_df <- data.frame(
        combined_ID = modified_rownames,
        current_IDs,
        stringsAsFactors = FALSE
      )

      # Add source information
      current_IDs_df$source_object <- mode_suffix

      # Combine with existing IDs
      if (nrow(combined_IDs) == 0) {
        combined_IDs <- current_IDs_df
      } else {
        # Ensure column names match for proper binding
        common_cols <- intersect(colnames(combined_IDs), colnames(current_IDs_df))

        # If columns don't match perfectly, we need to align them
        if (length(common_cols) < max(ncol(combined_IDs), ncol(current_IDs_df))) {
          # Get all unique column names
          all_cols <- unique(c(colnames(combined_IDs), colnames(current_IDs_df)))

          # Add missing columns with NA values
          for (col in all_cols) {
            if (!col %in% colnames(combined_IDs)) {
              combined_IDs[[col]] <- NA
            }
            if (!col %in% colnames(current_IDs_df)) {
              current_IDs_df[[col]] <- NA
            }
          }

          # Reorder columns to match
          combined_IDs <- combined_IDs[, all_cols]
          current_IDs_df <- current_IDs_df[, all_cols]
        }

        combined_IDs <- rbind(combined_IDs, current_IDs_df)
      }
    }
  }

  # Create merged original metadata
  merged_original_metadata <- NULL
  if (length(combined_original_metadata) > 0) {
    # Start with first object's metadata
    merged_original_metadata <- combined_original_metadata[[1]]

    if (length(combined_original_metadata) > 1) {
      for (i in 2:length(combined_original_metadata)) {
        current_meta <- combined_original_metadata[[i]]

        # Find common factors (rows)
        common_factors <- intersect(rownames(merged_original_metadata), rownames(current_meta))
        unique_to_current <- setdiff(rownames(current_meta), rownames(merged_original_metadata))

        # Add unique factors from current object
        if (length(unique_to_current) > 0) {
          # Find common samples (columns) to add new factors for
          common_samples <- intersect(colnames(merged_original_metadata), colnames(current_meta))

          if (length(common_samples) > 0) {
            new_rows <- current_meta[unique_to_current, common_samples, drop = FALSE]
            # Add NA columns for samples not in current object
            missing_samples <- setdiff(colnames(merged_original_metadata), colnames(current_meta))
            if (length(missing_samples) > 0) {
              na_cols <- matrix(NA, nrow = length(unique_to_current), ncol = length(missing_samples))
              colnames(na_cols) <- missing_samples
              rownames(na_cols) <- unique_to_current
              new_rows <- cbind(new_rows, na_cols)
              # Reorder columns to match merged_original_metadata
              new_rows <- new_rows[, colnames(merged_original_metadata), drop = FALSE]
            }
            merged_original_metadata <- rbind(merged_original_metadata, new_rows)
          }
        }

        # Check for discrepancies in common factors
        for (factor_name in common_factors) {
          common_samples_check <- intersect(colnames(merged_original_metadata), colnames(current_meta))
          if (length(common_samples_check) > 0) {
            merged_vals <- merged_original_metadata[factor_name, common_samples_check]
            current_vals <- current_meta[factor_name, common_samples_check]

            # Check for discrepancies (excluding NA values)
            non_na_indices <- !is.na(merged_vals) & !is.na(current_vals)
            if (any(non_na_indices)) {
              discrepancies <- merged_vals[non_na_indices] != current_vals[non_na_indices]
              if (any(discrepancies, na.rm = TRUE)) {
                warning(paste0("Discrepancies found in factor '", factor_name,
                               "' between objects. Keeping values from first object (",
                               names(combined_original_metadata)[1], ")."))
              }
            }
          }
        }
      }
    }
  }

  # Create combined romics object
  combined_romics_object <- createRomicsObject(norm_data, metadata = m, main_factor = f)

  # Add layers with proper handling of empty data
  if (nrow(missing_data) > 0) {
    combined_romics_object$missingdata <- missing_data
  }

  if (nrow(og_data) > 0) {
    combined_romics_object$original_data <- og_data
  }

  # Add original metadata layer
  if (!is.null(merged_original_metadata)) {
    combined_romics_object$original_metadata <- merged_original_metadata
  }

  # Add combined IDs layer if any objects had IDs
  if (nrow(combined_IDs) > 0) {
    # Set rownames to match the combined data rownames
    rownames(combined_IDs) <- combined_IDs$combined_ID
    # Remove the helper column
    combined_IDs$combined_ID <- NULL

    combined_romics_object$IDs <- combined_IDs
    message("Combined IDs layer created with ", nrow(combined_IDs), " entries")
    message("IDs have been modified to include source object identifier (@mode)")
  }

  # Add log transformation steps
  if (log[1] == 2) {
    combined_romics_object$steps <- c(combined_romics_object$steps,
                                      paste0("date|", gsub(" ", "_", format(Sys.time(), "%b_%d_%Y_%X")), "|log2transform"),
                                      "fun|log2transform(combined_romics_object)")
  }
  if (log[1] == 10) {
    combined_romics_object$steps <- c(combined_romics_object$steps,
                                      paste0("date|", gsub(" ", "_", format(Sys.time(), "%b_%d_%Y_%X")), "|log10transform"),
                                      "fun|log10transform(combined_romics_object)")
  }

  # Add combination step
  combined_romics_object$steps <- c(combined_romics_object$steps,
                                    paste0("date|", gsub(" ", "_", format(Sys.time(), "%b_%d_%Y_%X")), "|combineRomicsObjects"),
                                    paste0("fun|combineRomicsObjects(", paste(names(romics_objects), collapse = ", "), ")"))

  # Print summary
  message("Successfully combined ", length(romics_objects), " romics objects:")
  for (i in 1:length(romics_objects)) {
    message("  - ", names(romics_objects)[i], ": ", nrow(romics_objects[[i]]$data), " features")
  }
  message("Total features in combined object: ", nrow(combined_romics_object$data))

  # Report on column handling
  if (!is.null(common_original_cols) && length(all_original_cols) > 1) {
    all_orig_unique <- unique(unlist(all_original_cols))
    dropped_orig <- setdiff(all_orig_unique, common_original_cols)
    if (length(dropped_orig) > 0) {
      message("Dropped ", length(dropped_orig), " uncommon columns from original_data: ",
              paste(dropped_orig, collapse = ", "))
    }
  }

  if (!is.null(common_missing_cols) && length(all_missing_cols) > 1) {
    all_missing_unique <- unique(unlist(all_missing_cols))
    dropped_missing <- setdiff(all_missing_unique, common_missing_cols)
    if (length(dropped_missing) > 0) {
      message("Dropped ", length(dropped_missing), " uncommon columns from missingdata: ",
              paste(dropped_missing, collapse = ", "))
    }
  }

  return(combined_romics_object)
}

#' romicsCombineSpatialObjects()
#' @description Combines multiple spatial romics_objects into a single object by arranging them in a grid layout.
#' This function is designed for spatial omics data where samples have x,y coordinates. Each romics_object
#' is positioned in a grid with specified spacing, and coordinates are adjusted to prevent overlap.
#' Sample names are prefixed with object names to ensure uniqueness even when objects have identical sample names.
#' @param romics_objects A named list of romics_objects to be combined. All objects must have 'x' and 'y'
#'        coordinates in their metadata and should contain the same features (rows in data layer).
#'        Names will be used as prefixes for sample names.
#' @param xMinShift Numeric. Minimum spacing between objects along the x-axis. Default: 100
#' @param yMinShift Numeric. Minimum spacing between objects along the y-axis. Default: 100
#' @param xElements Integer. Number of objects to place horizontally before moving to next row. Default: 2
#' @param yElements Integer. Number of objects to place vertically. Default: 2
#' @param center_align Logical. If TRUE (default), objects are centered within their grid positions.
#'        If FALSE, objects are positioned at grid corner coordinates.
#' @param add_object_factor Logical. If TRUE (default), adds a factor indicating which
#'        original object each sample came from.
#' @param preserve_factors Logical. If TRUE (default), attempts to preserve all factors
#'        from input objects. Conflicts are resolved by prefixing with object names.
#' @details This function combines spatial romics objects by:
#'          1. Validating that all objects have compatible structure
#'          2. Calculating grid positions for each object
#'          3. Adjusting x,y coordinates to position objects in the grid
#'          4. Prefixing sample names with object names to ensure uniqueness
#'          5. Combining data, metadata, and other layers
#'          6. Optionally preserving factors from original objects
#'
#'          Grid positioning:
#'          - Objects are placed left-to-right, then bottom-to-top
#'          - Position (0,0) is bottom-left of the grid
#'          - Spacing ensures minimum distance between object boundaries
#'          - Center alignment positions objects at grid cell centers
#'
#'          Sample name handling:
#'          - All sample names are prefixed with their source object name
#'          - Example: "voxel1" from object "heart" becomes "heart_voxel1"
#'          - This prevents conflicts when objects have identical sample names
#'
#'          Coordinate transformation:
#'          - Each object's coordinates are shifted to fit in assigned grid position
#'          - Original relative positions within each object are preserved
#'          - Minimum bounding rectangles are calculated for proper spacing
#' @return Returns a combined romics_object containing all samples from input objects
#'         with adjusted spatial coordinates and prefixed sample names. The object will contain
#'         a new factor indicating source object (if add_object_factor = TRUE).
#' @author Geremy Clair
#' @export
romicsCombineSpatialObjects <- function(romics_objects, xMinShift = 100, yMinShift = 100,
                                        xElements = 2, yElements = 2, center_align = TRUE,
                                        add_object_factor = TRUE, preserve_factors = TRUE) {

  # Input validation
  if (missing(romics_objects) || !is.list(romics_objects) || length(romics_objects) == 0) {
    stop("romics_objects must be a non-empty list of romics_objects")
  }

  n_objects <- length(romics_objects)

  # Check if objects are named
  if (is.null(names(romics_objects))) {
    names(romics_objects) <- paste0("Object_", 1:n_objects)
    message("romics_objects list was not named. Using default names: Object_1, Object_2, etc.")
  } else if (any(names(romics_objects) == "" | is.na(names(romics_objects)))) {
    stop("All elements in romics_objects list must have non-empty names")
  }

  object_names <- names(romics_objects)

  # Validate grid dimensions
  xElements <- as.integer(xElements)
  yElements <- as.integer(yElements)
  max_objects <- xElements * yElements

  if (n_objects > max_objects) {
    stop(paste("Grid size (", xElements, "x", yElements, " = ", max_objects,
               ") is too small for ", n_objects, " objects", sep = ""))
  }

  # Validate romics objects and get common features
  all_features <- lapply(romics_objects, function(obj) {
    if (!is.romicsObject(obj)) {
      stop("All elements must be valid romics_objects")
    }
    if (!"x" %in% rownames(obj$metadata) || !"y" %in% rownames(obj$metadata)) {
      stop("All objects must have 'x' and 'y' coordinates in metadata")
    }
    return(rownames(obj$data))
  })

  common_features <- Reduce(intersect, all_features)
  if (length(common_features) == 0) {
    stop("No common features found across all objects")
  }

  # Calculate object dimensions and grid positions
  object_info <- list()

  for (i in seq_along(romics_objects)) {
    obj <- romics_objects[[i]]
    obj_name <- object_names[i]
    metadata_t <- data.frame(t(obj$metadata), stringsAsFactors = FALSE)

    px_x <- as.numeric(metadata_t$x)
    px_y <- as.numeric(metadata_t$y)
    valid_coords <- !is.na(px_x) & !is.na(px_y)

    if (sum(valid_coords) == 0) {
      stop(paste("Object", obj_name, "has no samples with valid coordinates"))
    }

    px_x_valid <- px_x[valid_coords]
    px_y_valid <- px_y[valid_coords]

    # Calculate bounding box
    x_min <- min(px_x_valid)
    x_max <- max(px_x_valid)
    y_min <- min(px_y_valid)
    y_max <- max(px_y_valid)

    width <- x_max - x_min
    height <- y_max - y_min

    # Calculate grid position
    grid_x <- (i - 1) %% xElements
    grid_y <- floor((i - 1) / xElements)

    object_info[[i]] <- list(
      obj = obj,
      name = obj_name,
      original_bounds = list(x_min = x_min, x_max = x_max, y_min = y_min, y_max = y_max),
      dimensions = list(width = width, height = height),
      grid_pos = list(x = grid_x, y = grid_y),
      valid_coords = valid_coords
    )
  }

  # Calculate grid cell dimensions
  max_width <- max(sapply(object_info, function(info) info$dimensions$width))
  max_height <- max(sapply(object_info, function(info) info$dimensions$height))

  cell_width <- max_width + xMinShift
  cell_height <- max_height + yMinShift

  # Process and combine all objects
  combined_data_list <- list()
  combined_metadata_list <- list()
  combined_missingdata_list <- list()
  combined_embeddings_list <- list()

  for (i in seq_along(object_info)) {
    info <- object_info[[i]]
    obj <- info$obj
    obj_name <- info$name

    # Calculate target position
    target_x <- info$grid_pos$x * cell_width
    target_y <- info$grid_pos$y * cell_height

    if (center_align) {
      target_x <- target_x + (cell_width - info$dimensions$width) / 2
      target_y <- target_y + (cell_height - info$dimensions$height) / 2
    }

    # Calculate coordinate shifts
    x_shift <- target_x - info$original_bounds$x_min
    y_shift <- target_y - info$original_bounds$y_min

    # Process metadata
    metadata_t <- data.frame(t(obj$metadata), stringsAsFactors = FALSE)
    metadata_t <- metadata_t[info$valid_coords, ]

    # Adjust coordinates
    metadata_t$x <- as.numeric(metadata_t$x) + x_shift
    metadata_t$y <- as.numeric(metadata_t$y) + y_shift

    # Add object identifier
    if (add_object_factor) {
      metadata_t$source_object <- obj_name
    }

    # Create new sample names with object prefix
    original_sample_names <- rownames(metadata_t)
    new_sample_names <- paste0(obj_name, "_", original_sample_names)
    rownames(metadata_t) <- new_sample_names

    # Process data layer
    obj_data <- obj$data[common_features, info$valid_coords, drop = FALSE]
    colnames(obj_data) <- new_sample_names

    # Process missing data layer
    obj_missingdata <- NULL
    if (!is.null(obj$missingdata)) {
      obj_missingdata <- obj$missingdata[common_features, info$valid_coords, drop = FALSE]
      colnames(obj_missingdata) <- new_sample_names
    }

    # Process embeddings
    obj_embeddings <- NULL
    if (!is.null(obj$embeddings)) {
      obj_embeddings <- obj$embeddings[, info$valid_coords, drop = FALSE]
      colnames(obj_embeddings) <- new_sample_names
    }

    # Store processed components
    combined_data_list[[i]] <- obj_data
    combined_metadata_list[[i]] <- metadata_t
    combined_missingdata_list[[i]] <- obj_missingdata
    combined_embeddings_list[[i]] <- obj_embeddings
  }

  # Combine all components
  combined_data <- do.call(cbind, combined_data_list)
  combined_metadata_df <- do.call(rbind, combined_metadata_list)
  combined_metadata <- data.frame(t(combined_metadata_df), stringsAsFactors = FALSE)

  # Combine missing data
  combined_missingdata <- NULL
  if (any(!sapply(combined_missingdata_list, is.null))) {
    for (i in seq_along(combined_missingdata_list)) {
      if (is.null(combined_missingdata_list[[i]])) {
        combined_missingdata_list[[i]] <- matrix(FALSE,
                                                 nrow = length(common_features),
                                                 ncol = ncol(combined_data_list[[i]]))
        rownames(combined_missingdata_list[[i]]) <- common_features
        colnames(combined_missingdata_list[[i]]) <- colnames(combined_data_list[[i]])
      }
    }
    combined_missingdata <- do.call(cbind, combined_missingdata_list)
  }

  # Combine embeddings
  combined_embeddings <- NULL
  if (any(!sapply(combined_embeddings_list, is.null))) {
    valid_embeddings <- combined_embeddings_list[!sapply(combined_embeddings_list, is.null)]
    if (length(valid_embeddings) > 0) {
      combined_embeddings <- do.call(cbind, valid_embeddings)
      message("Embeddings combined successfully")
    }
  }

  # Create combined romics object
  combined_object <- romics_objects[[1]]  # Use first object as template
  combined_object$data <- combined_data
  combined_object$metadata <- combined_metadata
  combined_object$missingdata <- combined_missingdata
  combined_object$embeddings <- combined_embeddings

  # Set main factor
  combined_object$main_factor <- if (add_object_factor) "source_object" else "x"

  # Clear statistics as they're no longer valid
  combined_object$statistics <- NULL

  # Update steps
  arguments <- as.list(match.call())
  combined_object <- romicsUpdateSteps(combined_object, arguments)

  # Report summary
  total_samples <- ncol(combined_data)
  message(paste("Successfully combined", n_objects, "objects into", total_samples, "total samples"))
  message(paste("Using", length(common_features), "common features"))
  message(paste("Grid arrangement:", xElements, "x", yElements, "with cell size",
                round(cell_width, 2), "x", round(cell_height, 2)))

  return(combined_object)
}
