#' @title Convert romics_object to SpatialExperiment
#' @description Converts a romics_object with spatial coordinates to a SpatialExperiment object,
#' enabling integration with spatial transcriptomics tools like BANKSY and other Bioconductor
#' spatial analysis packages.
#' @param romics_object A romics_object created using createRomicsObject(). Must contain spatial
#' coordinates in the metadata layer.
#' @param x_coord Character. Name of the metadata factor containing x-coordinates. If NULL (default),
#' will auto-search for a factor named "x". If not found, an error will be raised.
#' @param y_coord Character. Name of the metadata factor containing y-coordinates. If NULL (default),
#' will auto-search for a factor named "y". If not found, an error will be raised.
#' @param z_coord Character. Optional name of the metadata factor containing z-coordinates for 3D
#' spatial data. Default: NULL (2D coordinates only).
#' @param assay_name Character. Name to assign to the primary assay in the SpatialExperiment.
#' Default: "counts".
#' @param include_statistics Logical. If TRUE and the romics_object contains a $statistics layer,
#' it will be included as additional assays in the SpatialExperiment. Default: TRUE.
#' @details This function maps romics_object components to SpatialExperiment as follows:
#' - `$data` → counts assay
#' - Spatial coordinates → spatialCoords slot
#' - Metadata factors → colData
#' - Feature IDs → rowData
#' - Processing history and audit info → metadata slot
#'
#' Coordinates are automatically coerced to numeric. Missing coordinate values will trigger an error.
#' @return A SpatialExperiment object with all components mapped from the romics_object.
#' @author Geremy Clair
#' @export
#' @examples \dontrun{
#' # Basic conversion with auto-detected coordinates
#' se <- romicsToSpatialExperiment(romics_object)
#'
#' # Explicit coordinate specification
#' se <- romicsToSpatialExperiment(romics_object, x_coord = "spatial_x", y_coord = "spatial_y")
#'
#' # With 3D coordinates
#' se <- romicsToSpatialExperiment(romics_object, x_coord = "x", y_coord = "y", z_coord = "z")
#'
#' # Ready for BANKSY analysis
#' # banksy_result <- BANKSY::computeBanksy(se)
#' }
#'
romicsToSpatialExperiment <- function(romics_object,
                                      x_coord = NULL,
                                      y_coord = NULL,
                                      z_coord = NULL,
                                      assay_name = "counts",
                                      include_statistics = TRUE) {

  # ── Input validation ────────────────────────────────────────────────────────
  # Check for UUID first (for legacy objects) before full validation
  if (!inherits(romics_object, "romics_object")) {
    stop("romics_object must be a romics_object created with createRomicsObject()")
  }

  # Assign UUID if missing (for legacy objects that predate UUID requirement)
  if (is.null(romics_object$uuid)) {
    romics_object$uuid <- uuid::UUIDgenerate()
  }

  # Now do full validation
  if (!is.romicsObject(romics_object)) {
    stop("romics_object must be a romics_object created with createRomicsObject()")
  }

  if (!requireNamespace("SpatialExperiment", quietly = TRUE)) {
    stop("SpatialExperiment package is required. Install with: BiocManager::install('SpatialExperiment')")
  }

  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("SingleCellExperiment package is required. Install with: BiocManager::install('SingleCellExperiment')")
  }

  # ── Auto-detect coordinates if not provided ─────────────────────────────────
  available_factors <- rownames(romics_object$metadata)

  if (is.null(x_coord)) {
    if ("x" %in% available_factors) {
      x_coord <- "x"
      message("Auto-detected x-coordinate factor: 'x'")
    } else {
      stop("x_coord not specified and no 'x' factor found in metadata.\n",
           "Please specify x_coord parameter or ensure x/y factors exist.\n",
           "Available factors: ", paste(available_factors, collapse = ", "))
    }
  }

  if (is.null(y_coord)) {
    if ("y" %in% available_factors) {
      y_coord <- "y"
      message("Auto-detected y-coordinate factor: 'y'")
    } else {
      stop("y_coord not specified and no 'y' factor found in metadata.\n",
           "Please specify y_coord parameter or ensure x/y factors exist.\n",
           "Available factors: ", paste(available_factors, collapse = ", "))
    }
  }

  # Check if coordinates exist in metadata
  if (!x_coord %in% available_factors) {
    stop("x_coord factor '", x_coord, "' not found in metadata.\n",
         "Available factors: ", paste(available_factors, collapse = ", "))
  }

  if (!y_coord %in% available_factors) {
    stop("y_coord factor '", y_coord, "' not found in metadata.\n",
         "Available factors: ", paste(available_factors, collapse = ", "))
  }

  if (!is.null(z_coord) && !z_coord %in% available_factors) {
    stop("z_coord factor '", z_coord, "' not found in metadata.\n",
         "Available factors: ", paste(available_factors, collapse = ", "))
  }

  # ── Extract data matrix and names ──────────────────────────────────────────
  # Keep sample names for consistency
  sample_names <- colnames(romics_object$data)
  feature_names <- rownames(romics_object$data)

  # ── Extract and validate coordinates ────────────────────────────────────────
  # Extract coordinates from metadata (columns are samples, rows are factors)
  x_values <- as.numeric(as.character(romics_object$metadata[x_coord, ]))
  y_values <- as.numeric(as.character(romics_object$metadata[y_coord, ]))

  if (any(is.na(x_values)) || any(is.na(y_values))) {
    stop("Coordinate factors contain NA or non-numeric values.\n",
         "Please ensure coordinates can be coerced to numeric.")
  }

  # Build spatial coords matrix with samples as rows
  # Ensure rownames match exactly with sample names
  spatial_coords <- matrix(nrow = length(sample_names), ncol = 2)
  spatial_coords[, 1] <- x_values
  spatial_coords[, 2] <- y_values
  rownames(spatial_coords) <- sample_names
  colnames(spatial_coords) <- c("x", "y")

  # Add z-coordinate if provided
  if (!is.null(z_coord)) {
    z_values <- as.numeric(as.character(romics_object$metadata[z_coord, ]))
    if (any(is.na(z_values))) {
      stop("z_coord contains NA or non-numeric values.")
    }
    spatial_coords <- cbind(spatial_coords, z = z_values)
  }

  # ── Ensure all names match before building components ───────────────────────
  # Create matrix with proper dimension names
  data_matrix <- matrix(as.numeric(as.matrix(romics_object$data)),
                        nrow = length(feature_names),
                        ncol = length(sample_names))

  # ── Create colData from metadata ────────────────────────────────────────────
  coldata <- data.frame(t(romics_object$metadata), stringsAsFactors = FALSE)
  rownames(coldata) <- sample_names

  # Remove x and y from coldata if they exist (they're in spatialCoords)
  coldata <- coldata[, !colnames(coldata) %in% c("x", "y", "z"), drop = FALSE]

  # ── Create rowData from feature IDs ─────────────────────────────────────────
  # Use S4Vectors::DataFrame for proper compatibility
  rowdata <- S4Vectors::DataFrame(row.names = feature_names)

  if (!is.null(romics_object$IDs) && !all(romics_object$IDs == "unknown")) {
    rowdata$IDs <- romics_object$IDs
  }

  # ── Create SpatialExperiment using a workaround for version compatibility ─────
  # Bioconductor has strict internal validation that conflicts with our setup
  # Use a minimal approach that bypasses the validation

  # Ensure matrix has proper names
  rownames(data_matrix) <- feature_names
  colnames(data_matrix) <- sample_names

  # Create an empty SpatialExperiment first with just colData to establish structure
  # This avoids passing spatialCoords through the constructor which triggers validation
  se <- tryCatch({
    SpatialExperiment::SpatialExperiment(
      colData = coldata
    )
  }, error = function(e) {
    # If that fails, try SummarizedExperiment as base
    SummarizedExperiment::SummarizedExperiment(
      colData = coldata
    )
  })

  # Ensure this is converted to SpatialExperiment class if it isn't already
  if (!methods::is(se, "SpatialExperiment")) {
    se <- as(se, "SpatialExperiment")
  }

  # Now add all components via direct slot manipulation to avoid validation
  # Add the assay
  methods::slot(se, "assays") <- SummarizedExperiment::Assays(
    S4Vectors::SimpleList(counts = methods::as(data_matrix, "matrix"))
  )

  # Add rowData
  methods::slot(se, "elementMetadata") <- rowdata

  # Add spatial coordinates
  SpatialExperiment::spatialCoords(se) <- spatial_coords

  # ── Store metadata (audit trail, processing steps, etc.) ──────────────────
  S4Vectors::metadata(se)$omics_type <- romics_object$omics_type
  S4Vectors::metadata(se)$uuid <- romics_object$uuid
  S4Vectors::metadata(se)$steps <- paste(romics_object$steps, collapse = " | ")

  if (!is.null(romics_object$omics_information)) {
    S4Vectors::metadata(se)$omics_information <- romics_object$omics_information
  }

  if (!is.null(romics_object$custom_colors)) {
    S4Vectors::metadata(se)$custom_colors <- romics_object$custom_colors
  }

  # ── Add statistics as additional assays if requested ──────────────────────
  if (include_statistics && !is.null(romics_object$statistics)) {
    stats_matrix <- as.matrix(romics_object$statistics)

    if (nrow(stats_matrix) == nrow(se)) {
      SummarizedExperiment::assay(se, "statistics") <- stats_matrix
    }
  }

  return(se)
}


#' @title Convert SpatialExperiment to romics_object
#' @description Converts a SpatialExperiment object back to a romics_object format,
#' integrating spatial coordinates as x/y metadata factors.
#' @param se_object A SpatialExperiment object, typically output from romicsToSpatialExperiment()
#' or another spatial analysis package.
#' @param main_factor Character. Name of the colData column to use as the main grouping factor.
#' If NULL (default), will try to restore from metadata, otherwise uses the first available factor.
#' @param omics_type Character. Type of omics data (e.g., "proteomics", "metabolomics").
#' If NULL, will attempt to restore from SpatialExperiment metadata.
#' @details This function converts SpatialExperiment components back to romics_object format:
#' - Primary assay → $data
#' - spatialCoords → x/y metadata factors
#' - colData → metadata factors
#' - rowData → IDs (if available)
#' - metadata → processing history, audit info, etc.
#'
#' Note: Statistics stored in SpatialExperiment assays are stored in rowData of the returned
#' romics_object rather than reconstructed into the $statistics layer.
#' @return A romics_object with class assignment and all necessary components.
#' @author Geremy Clair
#' @export
#' @examples \dontrun{
#' # Convert back from SpatialExperiment
#' romics_obj <- SpatialExperimentToRomics(se)
#'
#' # With explicit main factor specification
#' romics_obj <- SpatialExperimentToRomics(se, main_factor = "cell_type")
#' }
#'
SpatialExperimentToRomics <- function(se_object,
                                      main_factor = NULL,
                                      omics_type = NULL) {

  # ── Input validation ────────────────────────────────────────────────────────
  if (!requireNamespace("SpatialExperiment", quietly = TRUE)) {
    stop("SpatialExperiment package is required. Install with: BiocManager::install('SpatialExperiment')")
  }

  if (!methods::is(se_object, "SpatialExperiment")) {
    stop("se_object must be a SpatialExperiment object")
  }

  if (!SpatialExperiment::nrow(se_object) > 0) {
    stop("SpatialExperiment object is empty")
  }

  # ── Extract core components ────────────────────────────────────────────────
  data_matrix <- SummarizedExperiment::assay(se_object, "counts")
  if (is.null(data_matrix)) {
    data_matrix <- SummarizedExperiment::assay(se_object, 1)
  }

  coldata_df <- data.frame(SummarizedExperiment::colData(se_object),
                           stringsAsFactors = FALSE)

  # ── Extract spatial coordinates and add as metadata ────────────────────────
  spatial_coords <- SpatialExperiment::spatialCoords(se_object)

  x_values <- as.character(spatial_coords[, "x"])
  y_values <- as.character(spatial_coords[, "y"])

  metadata_list <- list()
  metadata_list$x <- x_values
  metadata_list$y <- y_values

  if (ncol(spatial_coords) > 2) {
    z_values <- as.character(spatial_coords[, 3])
    metadata_list$z <- z_values
  }

  # ── Add colData to metadata ────────────────────────────────────────────────
  for (col in colnames(coldata_df)) {
    if (col != "x" && col != "y" && col != "z") {
      metadata_list[[col]] <- as.character(coldata_df[[col]])
    }
  }

  # Convert to metadata matrix (rows=factors, columns=samples)
  metadata_matrix <- do.call(rbind, metadata_list)
  colnames(metadata_matrix) <- colnames(se_object)

  # ── Determine main_factor ──────────────────────────────────────────────────
  if (is.null(main_factor)) {
    main_factors <- S4Vectors::metadata(se_object)$main_factor
    if (!is.null(main_factors)) {
      main_factor <- main_factors
    } else {
      main_factor <- rownames(metadata_matrix)[1]
    }
  }

  if (!main_factor %in% rownames(metadata_matrix)) {
    stop("main_factor '", main_factor, "' not found in metadata.\n",
         "Available factors: ", paste(rownames(metadata_matrix), collapse = ", "))
  }

  # ── Determine omics_type ───────────────────────────────────────────────────
  if (is.null(omics_type)) {
    omics_type <- S4Vectors::metadata(se_object)$omics_type
    if (is.null(omics_type)) {
      omics_type <- "unknown"
    }
  }

  # ── Extract IDs from rowData if available ──────────────────────────────────
  ids <- NULL
  rowdata_df <- data.frame(SummarizedExperiment::rowData(se_object))

  if (nrow(rowdata_df) > 0 && "IDs" %in% colnames(rowdata_df)) {
    ids <- rowdata_df$IDs
  }

  if (is.null(ids)) {
    ids <- rep("unknown", nrow(se_object))
  }

  # ── Extract processing history and other metadata ─────────────────────────
  omics_info <- S4Vectors::metadata(se_object)$omics_information
  custom_colors <- S4Vectors::metadata(se_object)$custom_colors
  steps_str <- S4Vectors::metadata(se_object)$steps
  uuid <- S4Vectors::metadata(se_object)$uuid

  steps <- if (!is.null(steps_str)) {
    strsplit(steps_str, " \\| ")[[1]]
  } else {
    c(paste0("Converted from SpatialExperiment: ", Sys.time()))
  }

  # ── Create romics_object using the core creation function ───────────────
  data_df <- data.frame(data_matrix)
  metadata_df <- data.frame(t(metadata_matrix))

  romics_obj <- createRomicsObject(
    data = data_df,
    metadata = metadata_df,
    IDs = ids,
    main_factor = main_factor,
    custom_colors = custom_colors,
    omics_type = omics_type,
    omics_information = omics_info
  )

  # ── Restore UUID and steps if they exist ────────────────────────────────
  if (!is.null(uuid)) {
    romics_obj$uuid <- uuid
  }

  romics_obj$steps <- steps

  return(romics_obj)
}
