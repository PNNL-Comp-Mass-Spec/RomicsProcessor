#' @title Convert romics_object to SpatialExperiment
#' @description Converts a romics_object with spatial coordinates to a SpatialExperiment object
#' with proper rownames and colnames for native BANKSY and other Bioconductor spatial analysis tools.
#' @param romics_object A romics_object created using createRomicsObject(). Must contain spatial
#' coordinates in the metadata layer.
#' @param x_coord Character. Name of the metadata factor containing x-coordinates. If NULL (default),
#' will auto-search for a factor named "x". If not found, an error will be raised.
#' @param y_coord Character. Name of the metadata factor containing y-coordinates. If NULL (default),
#' will auto-search for a factor named "y". If not found, an error will be raised.
#' @param z_coord Character. Optional name of the metadata factor containing z-coordinates for 3D
#' spatial data. Default: NULL (2D coordinates only).
#' @param assay_name Character. Name to assign to the primary assay in the SpatialExperiment.
#' Default: "counts" (Note: currently assay_name parameter is reserved for future use).
#' @param include_statistics Logical. If TRUE and the romics_object contains a $statistics layer,
#' it will be included as additional assays in the SpatialExperiment. Default: TRUE.
#' @param clusters_to_import Character vector. Names of metadata factors to import as cluster assignments in colData.
#' For example, c("leiden_clusters", "kmeans_clusters"). These will be converted to factors. Default: NULL (no clusters imported).
#' @details This function maps romics_object components to SpatialExperiment with proper dimension names:
#' - `$data` → counts assay with rownames (features) and colnames (samples)
#' - Spatial coordinates → spatialCoords slot with rownames (samples)
#' - Metadata factors → colData with rownames (samples)
#' - Feature IDs → rowData with rownames (features)
#' - Processing history and audit info → metadata slot
#' - PCA/UMAP embeddings (if present) → reducedDims as "romicsPCA" and "romicsUMAP"
#' - Cluster assignments (if specified) → colData as factors
#' - Statistics layer (if present) → additional assay named "statistics"
#' - in_tissue column → automatically added to colData (assumes all points in tissue if not present)
#'
#' All dimension names are validated before construction to ensure compatibility with BANKSY's
#' native functions like Banksy::computeBanksy(), Banksy::runBanksyPCA(), and Banksy::clusterBanksy().
#'
#' Coordinates are automatically coerced to numeric. Missing coordinate values will trigger an error.
#' Embeddings from the romics_object are renamed with "romics" prefix to distinguish them from
#' any embeddings computed directly on the SpatialExperiment.
#' @return A SpatialExperiment object with properly named dimensions, ready for BANKSY analysis.
#' @author Geremy Clair
#' @export
#' @examples \dontrun{
#' # Basic conversion with auto-detected coordinates
#' spe <- romicsToSpatialExperiment(romics_object)
#'
#' # With PCA/UMAP and cluster import
#' # (assumes romics_object has PCA, UMAP computed and "leiden_clusters" in metadata)
#' spe <- romicsToSpatialExperiment(
#'   romics_object,
#'   clusters_to_import = c("leiden_clusters", "main_factor")
#' )
#'
#' # Verify imported embeddings
#' SingleCellExperiment::reducedDimNames(spe)  # Shows "romicsPCA", "romicsUMAP"
#' colnames(SpatialExperiment::colData(spe))   # Shows imported cluster factors
#'
#' # Explicit coordinate specification
#' spe <- romicsToSpatialExperiment(romics_object, x_coord = "spatial_x", y_coord = "spatial_y")
#'
#' # With 3D coordinates
#' spe <- romicsToSpatialExperiment(romics_object, x_coord = "x", y_coord = "y", z_coord = "z")
#'
#' # Ready for BANKSY analysis (using imported embeddings if available)
#' spe <- Banksy::computeBanksy(spe, assay_name = "counts", compute_aggrAssay = TRUE)
#' spe <- Banksy::clusterBanksy(spe, use_pca = TRUE, k = 25, resolution = 1.0)
#' }
#'
romicsToSpatialExperiment <- function(romics_object,
                                      x_coord = NULL,
                                      y_coord = NULL,
                                      z_coord = NULL,
                                      assay_name = "counts",
                                      include_statistics = TRUE,
                                      clusters_to_import = NULL) {

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

  # Add in_tissue column if it doesn't exist (assume all points are in tissue)
  if (!"in_tissue" %in% colnames(coldata)) {
    coldata$in_tissue <- 1L
  }

  # ── Import specified cluster assignments ────────────────────────────────────
  if (!is.null(clusters_to_import)) {
    if (!is.character(clusters_to_import)) {
      stop("clusters_to_import must be a character vector of metadata factor names")
    }

    # Validate that requested cluster factors exist
    available_factors <- rownames(romics_object$metadata)
    missing_factors <- setdiff(clusters_to_import, available_factors)
    if (length(missing_factors) > 0) {
      warning("The following cluster factors were not found in metadata and will be skipped: ",
              paste(missing_factors, collapse = ", "))
      clusters_to_import <- intersect(clusters_to_import, available_factors)
    }

    # Import each cluster factor
    for (cluster_name in clusters_to_import) {
      if (!(cluster_name %in% colnames(coldata))) {
        # Extract from metadata and convert to factor
        cluster_values <- as.character(romics_object$metadata[cluster_name, ])
        coldata[[cluster_name]] <- as.factor(cluster_values)
        message("Imported cluster factor '", cluster_name, "' to colData")
      }
    }
  }

  # ── Create rowData from feature IDs ─────────────────────────────────────────
  # Create DataFrame with feature names as a column for BANKSY compatibility
  # Empty rowData with just rownames causes validation errors in Banksy
  rowdata <- S4Vectors::DataFrame(
    feature_name = feature_names,
    row.names = feature_names
  )

  if (!is.null(romics_object$IDs) && !all(romics_object$IDs == "unknown")) {
    rowdata$IDs <- romics_object$IDs
  }

  # ── Create SpatialExperiment with proper names for BANKSY compatibility ─────
  # Ensure all matrices have proper rownames and colnames BEFORE construction
  # This is critical for BANKSY and other spatial analysis tools

  rownames(data_matrix) <- feature_names
  colnames(data_matrix) <- sample_names

  # Ensure colData has proper rownames (samples as rows)
  stopifnot(nrow(coldata) == length(sample_names),
            rownames(coldata) == sample_names)

  # Ensure spatial coordinates have proper rownames (samples as rows)
  stopifnot(nrow(spatial_coords) == length(sample_names),
            rownames(spatial_coords) == sample_names)

  # Ensure rowData has proper rownames (features as rows)
  stopifnot(nrow(rowdata) == length(feature_names),
            rownames(rowdata) == feature_names)

  # Create SpatialExperiment directly with all components
  # All components must have matching, properly-named dimensions
  se <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = data_matrix),
    rowData = rowdata,
    colData = coldata,
    spatialCoords = spatial_coords
  )

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

  # ── Import PCA and UMAP embeddings if they exist ───────────────────────────
  if (!is.null(romics_object$embeddings)) {
    embedding_matrix <- romics_object$embeddings

    # Extract PCA components (from romics_object)
    pca_cols <- grep("^pca_component_", rownames(embedding_matrix), value = TRUE)
    if (length(pca_cols) > 0) {
      pca_matrix <- t(embedding_matrix[pca_cols, , drop = FALSE])
      colnames(pca_matrix) <- sub("pca_component_", "romicsPC", pca_cols)
      SingleCellExperiment::reducedDim(se, "romicsPCA") <- pca_matrix
      message("Imported ", ncol(pca_matrix), " PCA components to SpatialExperiment (stored as 'romicsPCA')")
    }

    # Extract UMAP components (from romics_object)
    umap_cols <- grep("^umap_component_", rownames(embedding_matrix), value = TRUE)
    if (length(umap_cols) > 0) {
      umap_matrix <- t(embedding_matrix[umap_cols, , drop = FALSE])
      colnames(umap_matrix) <- sub("umap_component_", "romicsUMAP", umap_cols)
      SingleCellExperiment::reducedDim(se, "romicsUMAP") <- umap_matrix
      message("Imported ", ncol(umap_matrix), " UMAP components to SpatialExperiment (stored as 'romicsUMAP')")
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
#' - reducedDims (romicsPCA, romicsUMAP) → $embeddings layer
#' - Other reducedDims → $embeddings layer with their original names
#'
#' Note: Statistics stored in SpatialExperiment assays are stored in rowData of the returned
#' romics_object rather than reconstructed into the $statistics layer.
#' Embeddings from reducedDims are restored to the $embeddings layer with their original prefixes.
#' @return A romics_object with class assignment and all necessary components.
#' @author Geremy Clair
#' @export
#' @examples \dontrun{
#' # Round-trip conversion preserves embeddings
#' # Start with romics_object with PCA and UMAP
#' spe <- romicsToSpatialExperiment(romics_object)
#'
#' # PCA and UMAP are stored as romicsPCA and romicsUMAP
#' SingleCellExperiment::reducedDimNames(spe)  # "romicsPCA", "romicsUMAP"
#'
#' # Convert back - embeddings are restored
#' romics_obj_restored <- SpatialExperimentToRomics(spe)
#'
#' # Verify embeddings were restored
#' head(romics_obj_restored$embeddings)  # Shows pca_component_*, umap_component_*
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

  # ── Import PCA and UMAP embeddings from reducedDims if they exist ───────────
  reduced_dim_names <- SingleCellExperiment::reducedDimNames(se_object)

  if (length(reduced_dim_names) > 0) {
    embedding_list <- list()

    # Import romicsPCA (convert back to pca_component_ format)
    if ("romicsPCA" %in% reduced_dim_names) {
      pca_matrix <- SingleCellExperiment::reducedDim(se_object, "romicsPCA")
      pca_t <- t(pca_matrix)
      # Rename columns back to pca_component_ format
      colnames(pca_t) <- sub("^romicsPC", "pca_component_", colnames(pca_t))
      embedding_list <- c(embedding_list, list(pca_t))
      message("Restored romicsPCA embeddings to $embeddings as pca_component_*")
    }

    # Import romicsUMAP (convert back to umap_component_ format)
    if ("romicsUMAP" %in% reduced_dim_names) {
      umap_matrix <- SingleCellExperiment::reducedDim(se_object, "romicsUMAP")
      umap_t <- t(umap_matrix)
      # Rename columns back to umap_component_ format
      colnames(umap_t) <- sub("^romicsUMAP", "umap_component_", colnames(umap_t))
      embedding_list <- c(embedding_list, list(umap_t))
      message("Restored romicsUMAP embeddings to $embeddings as umap_component_*")
    }

    # Import any other reducedDims (preserve their names)
    other_dims <- setdiff(reduced_dim_names, c("romicsPCA", "romicsUMAP"))
    if (length(other_dims) > 0) {
      for (dim_name in other_dims) {
        dim_matrix <- SingleCellExperiment::reducedDim(se_object, dim_name)
        dim_t <- t(dim_matrix)
        rownames(dim_t) <- paste0(dim_name, "_component_", 1:nrow(dim_t))
        embedding_list <- c(embedding_list, list(dim_t))
        message("Restored '", dim_name, "' embeddings to $embeddings")
      }
    }

    # Combine all embeddings into one matrix
    if (length(embedding_list) > 0) {
      romics_obj$embeddings <- do.call(rbind, embedding_list)
    }
  }

  return(romics_obj)
}


#' @title Convert romics_object to Seurat object
#' @description Converts a romics_object to a Seurat object with metadata and optional embeddings,
#' compatible with Seurat's analysis workflows and visualization functions.
#' @param romics_object A romics_object created using createRomicsObject().
#' @param assay_name Character. Name for the Seurat assay. Default: "RNA"
#' @param clusters_to_import Character vector. Names of metadata factors to import as active ident and metadata.
#' Default: NULL (imports main_factor only)
#' @param include_embeddings Logical. If TRUE and embeddings exist (PCA/UMAP), import them to Seurat reductions.
#' Default: TRUE
#' @details This function creates a Seurat object with:
#' - Data matrix as the assay
#' - Metadata factors as object metadata
#' - Optional PCA/UMAP as reductions (named "pca" and "umap")
#' - Processing history preserved in metadata
#' @return A Seurat object compatible with standard Seurat workflows
#' @author Geremy Clair
#' @export
#' @examples \dontrun{
#' # Basic conversion
#' seurat_obj <- romicsToSeurat(romics_object)
#'
#' # With embeddings
#' seurat_obj <- romicsToSeurat(romics_object, include_embeddings = TRUE)
#'
#' # Ready for Seurat analysis
#' Seurat::DimPlot(seurat_obj, reduction = "umap")
#' }
#'
romicsToSeurat <- function(romics_object,
                           assay_name = "RNA",
                           clusters_to_import = NULL,
                           include_embeddings = TRUE) {

  if (!is.romicsObject(romics_object) || missing(romics_object)) {
    stop("romics_object must be a romics_object created with createRomicsObject()")
  }

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required. Install with: install.packages('Seurat')")
  }

  # Extract components
  data_matrix <- romics_object$data
  metadata_df <- data.frame(t(romics_object$metadata), check.names = FALSE)

  # Create Seurat object
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = data_matrix,
    assay = assay_name,
    meta.data = metadata_df
  )

  # Set active ident to main_factor
  Seurat::Idents(seurat_obj) <- romics_object$main_factor

  # Add additional cluster metadata if specified
  if (!is.null(clusters_to_import)) {
    if (!is.character(clusters_to_import)) {
      stop("clusters_to_import must be a character vector of metadata factor names")
    }

    available_factors <- rownames(romics_object$metadata)
    valid_factors <- intersect(clusters_to_import, available_factors)
    if (length(valid_factors) < length(clusters_to_import)) {
      missing <- setdiff(clusters_to_import, available_factors)
      warning("Some cluster factors not found in metadata and skipped: ", paste(missing, collapse = ", "))
    }
  }

  # Import embeddings if requested
  if (include_embeddings && !is.null(romics_object$embeddings)) {
    embedding_matrix <- romics_object$embeddings

    # Import PCA
    pca_cols <- grep("^pca_component_", rownames(embedding_matrix), value = TRUE)
    if (length(pca_cols) > 0) {
      pca_matrix <- t(embedding_matrix[pca_cols, , drop = FALSE])
      colnames(pca_matrix) <- sub("pca_component_", "PC", pca_cols)
      seurat_obj[["pca"]] <- Seurat::CreateDimReducObject(embeddings = pca_matrix, key = "PC_")
      message("Imported ", ncol(pca_matrix), " PCA components to Seurat")
    }

    # Import UMAP
    umap_cols <- grep("^umap_component_", rownames(embedding_matrix), value = TRUE)
    if (length(umap_cols) > 0) {
      umap_matrix <- t(embedding_matrix[umap_cols, , drop = FALSE])
      colnames(umap_matrix) <- sub("umap_component_", "UMAP", umap_cols)
      seurat_obj[["umap"]] <- Seurat::CreateDimReducObject(embeddings = umap_matrix, key = "UMAP_")
      message("Imported ", ncol(umap_matrix), " UMAP components to Seurat")
    }
  }

  # Store metadata
  seurat_obj@misc$omics_type <- romics_object$omics_type
  seurat_obj@misc$uuid <- romics_object$uuid
  seurat_obj@misc$steps <- romics_object$steps

  message("Successfully converted romics_object to Seurat object")
  return(seurat_obj)
}


#' @title Convert Seurat object to romics_object
#' @description Converts a Seurat object back to a romics_object format,
#' preserving metadata and optional embeddings.
#' @param seurat_object A Seurat object
#' @param main_factor Character. Name of the metadata column to use as main factor.
#' If NULL, uses the active identity. Default: NULL
#' @param omics_type Character. Type of omics data (e.g., "proteomics", "metabolomics").
#' If NULL, attempts to retrieve from Seurat misc. Default: NULL
#' @details This function converts Seurat object components back to romics_object format:
#' - Assay data → $data
#' - Metadata → metadata factors
#' - Reductions (pca, umap) → $embeddings layer
#' - misc metadata → metadata/audit info
#' @return A romics_object with class assignment
#' @author Geremy Clair
#' @export
#' @examples \dontrun{
#' # Convert Seurat object back to romics
#' romics_obj <- SeuratToRomics(seurat_object)
#'
#' # With specific main factor
#' romics_obj <- SeuratToRomics(seurat_object, main_factor = "cell_type")
#' }
#'
SeuratToRomics <- function(seurat_object,
                           main_factor = NULL,
                           omics_type = NULL) {

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required. Install with: install.packages('Seurat')")
  }

  if (!methods::is(seurat_object, "Seurat")) {
    stop("seurat_object must be a Seurat object")
  }

  # Extract data from default assay
  # Use 'layer' argument for SeuratObject 5.0.0+, fall back to 'slot' for older versions
  tryCatch({
    data_matrix <- Seurat::GetAssayData(seurat_object, layer = "counts")
  }, error = function(e) {
    data_matrix <- Seurat::GetAssayData(seurat_object, slot = "counts")
  })

  if (is.null(data_matrix) || nrow(data_matrix) == 0) {
    tryCatch({
      data_matrix <- Seurat::GetAssayData(seurat_object, layer = "data")
    }, error = function(e) {
      data_matrix <- Seurat::GetAssayData(seurat_object, slot = "data")
    })
  }

  # Convert sparse matrix to dense if needed
  if (methods::is(data_matrix, "sparseMatrix")) {
    message("Converting sparse matrix to dense (", format(object.size(data_matrix), units = "auto"), ")")
    data_matrix <- as.matrix(data_matrix)
  }

  # Transpose if needed (Seurat has features x cells, romics needs cells x features)
  message("Generating the metadata matrix")
  # Check by comparing dimensions to metadata
  metadata_df <- seurat_object@meta.data
  n_cells <- nrow(metadata_df)
  cell_names <- rownames(metadata_df)

  # Transpose if needed (Seurat has features x cells, romics needs cells x features)
  # Check by comparing dimensions to metadata
  metadata_df <- seurat_object@meta.data
  n_cells <- nrow(metadata_df)
  cell_names <- rownames(metadata_df)

  # Transpose data if needed (Seurat: features x cells -> romics: cells x features)
  if (nrow(data_matrix) == n_cells) {
    message("Data is already in cells x features format")
  } else if (ncol(data_matrix) == n_cells) {
    message("Generating the data matrix")
    data_matrix <- t(data_matrix)
  } else {
    stop("Data matrix dimensions (", nrow(data_matrix), " x ", ncol(data_matrix),
         ") don't match number of cells (", n_cells, ")")
  }

  # Data is now cells x features
  cell_ids <- rownames(data_matrix)
  feature_names <- colnames(data_matrix)

  # Store the original Seurat metadata for later use
  seurat_metadata_stored <- data.frame(seurat_object@meta.data)

  # Determine main_factor from Seurat metadata
  if (is.null(main_factor)) {
    # Look for common identity/cluster column names in order of preference
    common_names <- c("seurat_clusters", "cluster", "ident", "cell_type", "celltype",
                      "CellType", "cell_annotation", "annotation")

    if(sum(common_names %in% colnames(seurat_metadata_stored))==0){
      message("No cell defining or cluster factor was found in the Seurat, the main_factor will remain undefined. The list of available factor is the following:")
      message(paste(colnames(seurat_metadata_stored),collapse = " ; "))
      seurat_metadata_stored$main_factor<-"undefined"
      main_factor<-"main_factor"
    }else{
      for (name in common_names) {
        if (name %in% colnames(seurat_metadata_stored)) {
          main_factor <- name
          message("Using Seurat metadata column '", name, "' as main factor")
          break
        }
      }}}

    if (!main_factor %in% colnames(seurat_metadata_stored)) {
    stop("main_factor '", main_factor, "' not found in Seurat metadata.\n",
         "Available metadata columns: ", paste(colnames(seurat_metadata_stored), collapse = ", "))
  }

  # Determine omics_type
  if (is.null(omics_type)) {
    omics_type <- seurat_object@misc$omics_type
    if (is.null(omics_type)) {
      omics_type <- "unknown"
    }
  }

  # Convert data to data frame: cells x (ID + features)
  data_df <- as.data.frame(data_matrix)

  #data_df <- cbind(ID = cell_ids, data_df)
  #and add the names of the features
  data_df <- as.data.frame(t(data_df))
  data_df <- cbind(IDs=rownames(data_df),data_df)

  #extract the X/Y/Z coordinates exists they are appended to the metadata
  if(!is.null(GetTissueCoordinates(seurat_object))){
  message("Adding spatial coordinates to the metadata")
  coordinates<- GetTissueCoordinates(seurat_object)
  coordinates<-coordinates[,colnames(coordinates) %in% c("x","y","z","X","Y","Z")]
  colnames(coordinates)<-tolower(colnames(coordinates))
  #now let's merge the coordinates with the metadata
  seurat_metadata_stored<-cbind(seurat_metadata_stored,coordinates)
  }

  message("Finalizing the romics_object generation")
  #Lets orient the object appropriately
  seurat_metadata_stored<-t(data.frame(seurat_metadata_stored))
  #and add the first column
  seurat_metadata_stored<-cbind(IDs=rownames(seurat_metadata_stored), seurat_metadata_stored)

  data<-data.frame(data_df)
  metadata<-data.frame(seurat_metadata_stored)

  # Create romics_object using createRomicsObject
  romics_obj <- createRomicsObject(
    data = data,
    metadata = metadata,
    main_factor = main_factor,
    omics_type = omics_type
  )


  # Import reductions as embeddings
  reduction_names <- Seurat::Reductions(seurat_object)

  if (length(reduction_names) > 0) {
    embedding_list <- list()

    # Convert pca reduction
    if ("pca" %in% reduction_names) {
      pca_matrix <- Seurat::Embeddings(seurat_object, reduction = "pca")
      pca_t <- t(pca_matrix)
      colnames(pca_t) <- paste0("pca_component_", 1:ncol(pca_t))
      embedding_list[[length(embedding_list) + 1]] <- pca_t
      message("Restored PCA reduction to embeddings")
    }

    # Convert umap reduction
    if ("umap" %in% reduction_names) {
      umap_matrix <- Seurat::Embeddings(seurat_object, reduction = "umap")
      umap_t <- t(umap_matrix)
      colnames(umap_t) <- paste0("umap_component_", 1:ncol(umap_t))
      embedding_list[[length(embedding_list) + 1]] <- umap_t
      message("Restored UMAP reduction to embeddings")
    }

    # Convert other reductions
    other_reductions <- setdiff(reduction_names, c("pca", "umap"))
    if (length(other_reductions) > 0) {
      for (red_name in other_reductions) {
        red_matrix <- Seurat::Embeddings(seurat_object, reduction = red_name)
        red_t <- t(red_matrix)
        rownames(red_t) <- paste0(red_name, "_component_", 1:nrow(red_t))
        embedding_list[[length(embedding_list) + 1]] <- red_t
        message("Restored '", red_name, "' reduction to embeddings")
      }
    }

    if (length(embedding_list) > 0) {
      romics_obj$embeddings <- do.call(rbind, embedding_list)
    }
  }

  message("Successfully converted Seurat object to romics_object")
  return(romics_obj)
}


#' @title Convert romics_object to pmartR object
#' @description Converts a romics_object to a pmartR object (requires the pmartR package).
#' This enables use of pmartR's powerful omics data analysis and visualization functions.
#' @param romics_object A romics_object created using createRomicsObject()
#' @param data_type Character. The type of pmartR object to create: "proData" (proteomics),
#' "pepData" (peptide), "lipidData" (lipidomics), or "metabData" (metabolomics). Default: "proData"
#' @details This function creates a pmartR object with proper dimensions:
#' - Feature IDs from romics_object rownames
#' - Expression data from romics_object$data
#' - Sample metadata from romics_object$metadata
#' - Data scale (log vs abundance) determined from processing history
#' @return A pmartR object of the specified type
#' @author Geremy Clair
#' @export
#' @examples \dontrun{
#' # Convert to proteomics pmartR object
#' pmartR_obj <- romicsToPmartR(romics_object, data_type = "proData")
#'
#' # Convert to metabolomics pmartR object
#' pmartR_obj <- romicsToPmartR(romics_object, data_type = "metabData")
#'
#' # Ready for pmartR analysis
#' pmartR_obj <- pmartR::edata_transform(pmartR_obj, data_scale = "log2")
#' }
#'
romicsToPmartR <- function(romics_object, data_type = "proData") {

  if (!is.romicsObject(romics_object) || missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }

  if (!requireNamespace("pmartR", quietly = TRUE)) {
    stop("The pmartR package is required. Please install it with: install.packages('pmartR')")
  }

  if (!data_type %in% c("lipidData", "proData", "pepData", "metabData")) {
    stop("data_type must be one of: 'lipidData', 'proData', 'pepData', 'metabData'")
  }

  # Create expression data matrix with Mass_Tag_ID
  e_data <- data.frame(
    Mass_Tag_ID = 1:nrow(romics_object$data),
    romics_object$data,
    check.names = FALSE
  )
  rownames(e_data) <- NULL

  # Create sample metadata with SampleID
  f_data <- data.frame(
    SampleID = colnames(romics_object$metadata),
    t(romics_object$metadata),
    check.names = FALSE
  )
  rownames(f_data) <- NULL

  # Create feature metadata with IDs
  e_meta <- data.frame(
    Mass_Tag_ID = e_data$Mass_Tag_ID,
    IDs = rownames(romics_object$data),
    check.names = FALSE
  )
  rownames(e_meta) <- NULL

  # Determine data scale from processing history
  data_scale <- "abundance"
  if (sum(grepl("log10transform|log2transform", paste(romics_object$steps, collapse = " "))) > 0) {
    data_scale <- "log"
  }

  # Create pmartR object based on type
  pmartR_obj <- switch(data_type,
    "pepData" = pmartR::as.pepData(
      e_data = e_data,
      f_data = f_data,
      e_meta = e_meta,
      edata_cname = "Mass_Tag_ID",
      fdata_cname = "SampleID",
      emeta_cname = "Mass_Tag_ID",
      data_scale = data_scale,
    ),
    "proData" = pmartR::as.proData(
      e_data = e_data,
      f_data = f_data,
      e_meta = e_meta,
      edata_cname = "Mass_Tag_ID",
      fdata_cname = "SampleID",
      emeta_cname = "Mass_Tag_ID",
      data_scale = data_scale
    ),
    "lipidData" = pmartR::as.lipidData(
      e_data = e_data,
      f_data = f_data,
      e_meta = e_meta,
      edata_cname = "Mass_Tag_ID",
      fdata_cname = "SampleID",
      emeta_cname = "Mass_Tag_ID",
      data_scale = data_scale
    ),
    "metabData" = pmartR::as.metabData(
      e_data = e_data,
      f_data = f_data,
      e_meta = e_meta,
      edata_cname = "Mass_Tag_ID",
      fdata_cname = "SampleID",
      emeta_cname = "Mass_Tag_ID",
      data_scale = data_scale
    )
  )

  message("Successfully created ", class(pmartR_obj)[1], " pmartR object")
  return(pmartR_obj)
}
