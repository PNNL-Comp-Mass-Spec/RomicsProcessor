#' plotMap()
#' @description This function enables to create a spatial map of the features contained in an romics_object. It can either represent segments contained in a specific factor or features data. To ensure that the maps conserve appropriate proportionality all maps generated using this function are squares.
#' @param romics_object A romics_object created using romicsCreateObject(), the metadata has to contain at least a 'x' and 'y' coordinate in the factor layer in order to plot the map
#' @param factor The rowname of a metadata factor to be used for the plotting of regions associated with a given factor, if a factor is selected the feature map won't be shown.
#' @param factor_colors A character vector containing a list of colors to be used to represent the different levels of the factor. Default romics_colors will be used if not specified.
#' @param ROI A numeric vector containing 3 parameters, c(i,ii,iii): (i) 'x axis minimum' indicating the starting of the map on the x axis, (ii) 'y axis minimum' indicates the starting of the map on the y axis, and (iii) the lateral distance indicating the width and height of the map, to ensure that the maps conserve appropriate proportionality all maps generated using this function are squares.
#' @param feature A character vector containing the feature name. The feature name has to be the rowname of the 'romics_object$data' layer.
#' @param point_size A numeric vector indicating the size of the points generated for the map (1 by default)
#' @param point_alpha A numeric vector indicating the alpha of the points generated for the map (has to be comprised between 0 and 1, will be 1 = opaque by default).
#' @param scale_intensities Boolean. Whether to scale feature intensities. Default: TRUE
#' @param scale_colors A character vector that indicates the colors that should be used for intensity maps by default the viridis colors will be used
#' @param scale_limits A numerical vector that indicates the limits to be used for the color_scale, the function will automatically use the "squish" option that will make the points below the minimum or above the maximum appear with the extreme value of the scale.
#' @param na_value A color that indicates how NA values should be displayed when a color scale is used.
#' @param plotly Boolean. If TRUE, generates an interactive plotly plot instead of a static ggplot. Default: FALSE
#' @return This function generate spatial maps of spatial omics data as ggplot2 or plotly interactive plots.
#' @author Geremy Clair
#' @export
plotMap <- function(romics_object,
                     factor = NULL,
                     factor_colors = ROP_colors,
                     ROI = NULL,
                     feature = NULL,
                     point_size = 1,
                     point_alpha = 1,
                     scale_intensities = TRUE,
                     scale_colors = viridis::viridis(3),
                     scale_limits = NULL,
                     na_value = "gray20",
                     plotly = FALSE){

  # Input validation
  if(!is.romicsObject(romics_object)){
    stop("The object was not created with the function 'romicsCreateObject()'.")
  }

  # Parameter defaults and validation
  if(missing(factor)){factor <- NULL}
  if(missing(ROI)){ROI <- NULL}
  if(missing(feature)){feature <- NULL}
  if(missing(factor_colors)){factor_colors <- ROP_colors}
  if(missing(point_size)){point_size <- 1}
  if(missing(point_alpha)){point_alpha <- 1}
  if(missing(scale_intensities)){scale_intensities <- TRUE}
  if(missing(scale_colors)){scale_colors <- viridis::viridis(3)}
  if(missing(na_value)){na_value <- "gray20"}
  if(missing(plotly)){plotly <- FALSE}

  # Validate plotly parameter
  if(!is.logical(plotly) || length(plotly) != 1) {
    stop("plotly must be TRUE or FALSE")
  }

  # Import metadata
  m <- data.frame(t(romics_object$metadata))

  # Ensure that metadata coordinates are present and numerical
  if("x" %in% colnames(m)){
    m$x <- as.numeric(m$x)
  } else {
    stop("The coordinate 'x' is not a factor of the metadata.")
  }

  if("y" %in% colnames(m)){
    m$y <- as.numeric(m$y)
  } else {
    stop("The coordinate 'y' is not a factor of the metadata.")
  }

  # Calculate plot boundaries and filter data if ROI is specified
  if(!is.null(ROI)){
    if(length(ROI) != 3 || !is.numeric(ROI)) {
      stop("ROI must be a numeric vector of length 3: c(x_min, y_min, distance)")
    }
    xlim_bounds <- c(ROI[1], ROI[3] + ROI[1])
    ylim_bounds <- c(ROI[2], ROI[3] + ROI[2])
    # Filter data to ROI - keep points within the ROI boundaries
    roi_filter <- m$x >= ROI[1] & m$x <= (ROI[3] + ROI[1]) &
      m$y >= ROI[2] & m$y <= (ROI[3] + ROI[2])
    m_filtered <- m[roi_filter, ]
  } else {
    # Auto-calculate square boundaries
    xmin <- min(m$x, na.rm = TRUE)
    ymin <- min(m$y, na.rm = TRUE)
    xmax <- max(m$x, na.rm = TRUE)
    ymax <- max(m$y, na.rm = TRUE)
    x_amp <- xmax - xmin
    y_amp <- ymax - ymin
    amp <- max(c(x_amp, y_amp))
    xlim_bounds <- c(xmin, xmin + amp)
    ylim_bounds <- c(ymin, ymin + amp)
    # No filtering needed
    m_filtered <- m
    roi_filter <- rep(TRUE, nrow(m))
  }

  # Handle different plotting scenarios
  if(is.null(factor) && is.null(feature)){
    # Default white points
    p <- ggplot2::ggplot(m_filtered, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_point(colour = "white", size = point_size, alpha = point_alpha, stroke = NA) +
      theme_map() +
      ggplot2::xlim(xlim_bounds) +
      ggplot2::ylim(ylim_bounds)

  } else if(!is.null(factor)){
    # Factor-based coloring

    # Check if the selected factor is different from the main factor
    current_main_factor <- romics_object$main_factor

    if(factor != current_main_factor) {
      message(paste0("The selected factor '", factor, "' is not the main factor ('", current_main_factor, "')."))
      message("Changing romics_object factor to ensure consistent colors...")
      romics_object <- romicsChangeFactor(romics_object, main_factor = factor)
      # Re-import metadata after factor change
      m <- data.frame(t(romics_object$metadata))
      # Re-ensure coordinates are numerical
      m$x <- as.numeric(m$x)
      m$y <- as.numeric(m$y)
      # Re-apply ROI filtering if needed
      if(!is.null(ROI)){
        roi_filter <- m$x >= ROI[1] & m$x <= (ROI[3] + ROI[1]) &
          m$y >= ROI[2] & m$y <= (ROI[3] + ROI[2])
        m_filtered <- m[roi_filter, ]
      } else {
        m_filtered <- m
        roi_filter <- rep(TRUE, nrow(m))
      }
    }

    if(suppressWarnings({sum(is.na(as.numeric(m_filtered[, colnames(m_filtered) == factor]))) > 0})){
      message(paste0("The factor ", factor, " is non-numerical, it will be considered categorical"))

      # Get factor values for the filtered data
      factor_values <- m_filtered[, colnames(m_filtered) == factor]
      unique_factor_levels <- unique(factor_values)

      # Use provided factor_colors, or fall back to romics_object colors
      if(!identical(factor_colors, ROP_colors)) {
        # User provided custom colors
        color_mapping <- factor_colors[1:length(unique_factor_levels)]
        names(color_mapping) <- unique_factor_levels
      } else {
        # Use colors from the romics_object
        colors_romics <- as.character(t(romics_object$metadata["colors_romics", ]))

        # Filter colors to match ROI if applied
        if(!is.null(ROI)){
          colors_filtered <- colors_romics[roi_filter]
        } else {
          colors_filtered <- colors_romics
        }

        # Create a proper color mapping
        color_mapping <- c()
        for(level in unique_factor_levels) {
          # Get the color for the first occurrence of this level
          level_color <- colors_filtered[factor_values == level][1]
          color_mapping[level] <- level_color
        }
      }

      p <- ggplot2::ggplot(m_filtered, ggplot2::aes(x = x, y = y, colour = factor(factor_values, levels = unique_factor_levels))) +
        ggplot2::geom_point(size = point_size, alpha = point_alpha, stroke = NA) +
        theme_map() +
        ggplot2::xlim(xlim_bounds) +
        ggplot2::ylim(ylim_bounds) +
        ggplot2::scale_color_manual(values = color_mapping, drop = FALSE) +
        ggplot2::labs(colour = factor)

    } else {
      message(paste0("The factor ", factor, " is numerical, a gradient will be generated from its values"))
      p <- ggplot2::ggplot(m_filtered, ggplot2::aes(x = x, y = y, colour = as.numeric(m_filtered[, colnames(m_filtered) == factor]))) +
        ggplot2::geom_point(size = point_size, alpha = point_alpha, stroke = NA) +
        theme_map() +
        ggplot2::xlim(xlim_bounds) +
        ggplot2::ylim(ylim_bounds) +
        ggplot2::scale_colour_gradientn(colours = scale_colors, na.value = na_value) +
        ggplot2::labs(colour = factor)
    }

  } else if(is.null(factor) && !is.null(feature)){
    # Feature-based coloring
    # Validate feature exists
    if(!feature %in% rownames(romics_object$data)) {
      stop(paste("Feature", feature, "not found in romics_object data."))
    }

    # Extract feature data for ALL points first
    feature_data <- romics_object$data[rownames(romics_object$data) == feature, , drop = FALSE]
    feature_values_all <- as.numeric(feature_data[1, ])

    # Then filter to match the ROI-filtered metadata
    feature_values <- feature_values_all[roi_filter]

    # Debug: print some info
    message(paste("Feature values range:", round(min(feature_values, na.rm = TRUE), 3), "to", round(max(feature_values, na.rm = TRUE), 3)))
    message(paste("Number of points in ROI:", nrow(m_filtered)))
    message(paste("Number of feature values:", length(feature_values)))

    # Scale feature if requested
    if(scale_intensities == TRUE){
      feature_values <- as.numeric(scale(feature_values))
      message(paste("After scaling, range:", round(min(feature_values, na.rm = TRUE), 3), "to", round(max(feature_values, na.rm = TRUE), 3)))
    }

    # Ensure same length as filtered coordinate data
    if(length(feature_values) != nrow(m_filtered)) {
      stop(paste("Mismatch: feature has", length(feature_values), "values but filtered metadata has", nrow(m_filtered), "rows"))
    }

    # Add scaled values to filtered dataframe
    m_filtered$feature_values <- feature_values

    # Create plot
    p <- ggplot2::ggplot(m_filtered, ggplot2::aes(x = x, y = y, colour = feature_values)) +
      ggplot2::geom_point(size = point_size, alpha = point_alpha, stroke = NA) +
      theme_map() +
      ggplot2::xlim(xlim_bounds) +
      ggplot2::ylim(ylim_bounds) +
      ggplot2::ggtitle(label = feature) +
      ggplot2::labs(colour = if(scale_intensities) paste(feature, "(scaled)") else feature)

    # Add color scale
    if(is.null(scale_limits)){
      p <- p + ggplot2::scale_colour_gradientn(colours = scale_colors, na.value = na_value)
    } else {
      p <- p + ggplot2::scale_colour_gradientn(colours = scale_colors,
                                               limits = scale_limits,
                                               oob = scales::squish,
                                               na.value = na_value)
    }
  }

  # Return appropriate format
  if(plotly) {
    if(requireNamespace("plotly", quietly = TRUE)) {
      return(plotly::ggplotly(p))
    } else {
      warning("plotly package not available. Returning ggplot instead.")
      return(p)
    }
  } else {
    return(p)
  }
}

#' multiFeatureMap()
#' @description Creates a spatial map displaying multiple features from a romics_object simultaneously, with each feature represented by different colors and intensities mapped to transparency.
#' @param romics_object A romics_object created using romicsCreateObject(), the metadata must contain 'x' and 'y' coordinates.
#' @param features A character vector containing the feature names to plot. Must match row names in the data.
#' @param colors A character vector of colors to represent different features. If NULL, uses default color palette. Default: NULL
#' @param scale_intensities Logical. Whether to scale feature values (recommended for multiple features). Default: TRUE
#' @param ROI A numeric vector of length 3: c(x_min, y_min, distance) defining the region of interest. If NULL, uses full data range. Default: NULL
#' @param point_size Numeric. Size of the points. Default: 1
#' @param alpha_range Numeric vector of length 2. Min and max transparency values c(min, max). Values between 0-1. Default: c(0.1, 0.8)
#' @param alpha_limits Numeric vector of length 2. Intensity values corresponding to min/max transparency. If NULL, uses data range. Default: NULL
#' @param facet Logical. Whether to create separate panels for each feature. Default: FALSE
#' @param na_value Color for NA values. Default: "gray90"
#' @param plotly Logical. If TRUE, creates an interactive plotly plot. Default: FALSE
#' @param title Character. Plot title. If NULL, uses auto-generated title. Default: NULL
#' @param white_background Logical. If TRUE, uses white background with black grid. If FALSE, uses black background with white grid. Default: TRUE
#' @param show_point_locations Logical. Whether to show background points indicating all point locations. Default: FALSE
#' @details This function overlays multiple features on a single spatial map, using color to distinguish features and transparency to represent intensity values. Feature names with spaces and special characters are properly handled.
#' @return Returns either a ggplot2 or plotly interactive plot showing the spatial distribution of multiple features.
#' @author Geremy Clair
#' @export
multiFeatureMap <- function(romics_object,
                            features = NULL,
                            colors = NULL,
                            scale_intensities = TRUE,
                            ROI = NULL,
                            point_size = 1,
                            alpha_range = c(0.1, 0.8),
                            alpha_limits = NULL,
                            facet = FALSE,
                            na_value = "gray90",
                            plotly = FALSE,
                            title = NULL,
                            white_background = TRUE,
                            show_point_locations = FALSE){

  # Check required packages
  required_packages <- c("ggplot2", "reshape2")
  if(plotly) {
    required_packages <- c(required_packages, "plotly")
  }

  for(pkg in required_packages) {
    if(!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required. Please install it with: install.packages('", pkg, "')", sep = ""))
    }
  }

  # Input validation (same as before...)
  if(!is.romicsObject(romics_object)){
    stop("The object was not created with the function 'romicsCreateObject()'.")
  }

  if(is.null(features) || missing(features)) {
    stop("At least one feature must be provided in the 'features' parameter")
  }

  if(!is.character(features)) {
    stop("features must be a character vector")
  }

  # Validate features exist
  missing_features <- features[!features %in% rownames(romics_object$data)]
  if(length(missing_features) > 0) {
    stop(paste("Feature(s) not found in data:", paste(missing_features, collapse = ", ")))
  }

  # Additional parameter validation
  if(!is.logical(white_background) || length(white_background) != 1) {
    stop("white_background must be TRUE or FALSE")
  }

  if(!is.logical(show_point_locations) || length(show_point_locations) != 1) {
    stop("show_point_locations must be TRUE or FALSE")
  }

  # Set default colors
  if(is.null(colors)) {
    if(exists("ROP_colors")) {
      colors <- ROP_colors[1:length(features)]
    } else {
      colors <- scales::hue_pal()(length(features))
    }
  } else if(length(colors) < length(features)) {
    warning("Fewer colors provided than features. Extending with default colors.")
    additional_colors_needed <- length(features) - length(colors)
    additional_colors <- scales::hue_pal()(additional_colors_needed)
    colors <- c(colors, additional_colors)
  }

  message(paste("Creating multi-feature map for", length(features), "features"))

  # Import metadata
  m <- data.frame(t(romics_object$metadata))

  # Ensure coordinates are present and numerical
  if(!"x" %in% colnames(m)){
    stop("The coordinate 'x' is not present in the metadata.")
  }
  if(!"y" %in% colnames(m)){
    stop("The coordinate 'y' is not present in the metadata.")
  }

  m$x <- as.numeric(m$x)
  m$y <- as.numeric(m$y)

  # Extract coordinates
  coordinates <- m[, c("x", "y")]

  # FIXED: Proper handling of feature names with spaces
  # Create safe mapping for feature names
  safe_feature_names <- make.names(features, unique = TRUE)
  feature_mapping <- setNames(features, safe_feature_names)

  # Extract feature data directly as matrix to preserve names
  feature_matrix <- romics_object$data[features, , drop = FALSE]

  # Scale features if requested
  if(scale_intensities) {
    feature_matrix <- t(scale(t(feature_matrix)))
    message("Feature intensities have been scaled (mean=0, sd=1)")
  }

  # Transpose and convert to data frame with safe names
  feature_data <- data.frame(t(feature_matrix))
  colnames(feature_data) <- safe_feature_names

  # Combine coordinates and feature data
  plot_data <- cbind(coordinates, feature_data)

  # Create long format using safe names
  plot_data_long <- reshape2::melt(plot_data, id.vars = c("x", "y"),
                                   variable.name = "feature_safe",
                                   value.name = "intensity")

  # Map back to original feature names
  plot_data_long$feature <- feature_mapping[as.character(plot_data_long$feature_safe)]
  plot_data_long$feature_safe <- NULL

  # Ensure proper data types
  plot_data_long$x <- as.numeric(plot_data_long$x)
  plot_data_long$y <- as.numeric(plot_data_long$y)
  plot_data_long$intensity <- as.numeric(plot_data_long$intensity)
  plot_data_long$feature <- as.factor(plot_data_long$feature)

  # Set alpha limits
  if(is.null(alpha_limits)) {
    if(scale_intensities) {
      alpha_limits <- c(-2, 2)
    } else {
      alpha_limits <- range(plot_data_long$intensity, na.rm = TRUE)
    }
  }

  # Set title
  if(is.null(title)) {
    if(length(features) <= 3) {
      title <- paste("Multi-feature Map:", paste(features, collapse = ", "))
    } else {
      title <- paste("Multi-feature Map (", length(features), "features)")
    }
  }

  # UPDATED: Create base plot with customizable theme
  p <- ggplot2::ggplot(plot_data_long, ggplot2::aes(x = x, y = y))

  # Apply theme based on background choice
  if(white_background) {
    p <- p + theme_map() +  # Assuming theme_map() has white background
      ggplot2::theme(
        panel.background = ggplot2::element_rect(fill = "white"),
        panel.grid.major = ggplot2::element_line(color = "black", size = 0.2),
        panel.grid.minor = ggplot2::element_line(color = "black", size = 0.1)
      )
    background_color <- "gray90"
  } else {
    p <- p + theme_map() +
      ggplot2::theme(
        panel.background = ggplot2::element_rect(fill = "black"),
        panel.grid.major = ggplot2::element_line(color = "white", size = 0.2),
        panel.grid.minor = ggplot2::element_line(color = "white", size = 0.1)
      )
    background_color <- "white"
  }

  p <- p + ggplot2::ggtitle(title) +
    ggplot2::labs(x = "X Coordinate", y = "Y Coordinate",
                  color = "Feature", alpha = "Intensity")

  # UPDATED: Add background points only if requested
  if(show_point_locations) {
    p <- p + ggplot2::geom_point(color = background_color, size = point_size,
                                 alpha = 0.3, stroke = 0)  # No outline
  }

  # UPDATED: Add feature points with no outline (stroke = 0)
  p <- p +
    ggplot2::geom_point(ggplot2::aes(color = feature, alpha = intensity),
                        size = point_size, stroke = 0) +  # No outline
    ggplot2::scale_alpha_continuous(range = alpha_range, limits = alpha_limits,
                                    na.value = 0, oob = scales::squish) +
    ggplot2::scale_color_manual(values = setNames(colors, features))

  # Set plot boundaries
  if(!is.null(ROI)){
    p <- p + ggplot2::xlim(ROI[1], ROI[1] + ROI[3]) +
      ggplot2::ylim(ROI[2], ROI[2] + ROI[3])
  } else {
    # Auto-calculate square boundaries
    xmin <- min(coordinates$x, na.rm = TRUE)
    ymin <- min(coordinates$y, na.rm = TRUE)
    xmax <- max(coordinates$x, na.rm = TRUE)
    ymax <- max(coordinates$y, na.rm = TRUE)
    x_amp <- xmax - xmin
    y_amp <- ymax - ymin
    amp <- max(c(x_amp, y_amp))
    p <- p + ggplot2::xlim(xmin, xmin + amp) +
      ggplot2::ylim(ymin, ymin + amp)
  }

  # Add faceting if requested
  if(facet) {
    p <- p + ggplot2::facet_wrap(~feature)
  }

  # Print summary
  message("Plot summary:")
  message(paste("  Features plotted:", paste(features, collapse = ", ")))
  message(paste("  Background:", ifelse(white_background, "white", "black")))
  message(paste("  Show locations:", show_point_locations))

  # Return appropriate format
  if(plotly) {
    return(plotly::ggplotly(p))
  } else {
    return(p)
  }
}

#' romicsCreateFactorFromROIs()
#' @description Creates metadata factors from 2D Regions of Interest (ROI) definitions.
#' This function takes a list of ROI coordinates and creates a new factor in the romics_object metadata
#' indicating which samples fall within each ROI. ROIs are defined as rectangular regions using
#' minimum coordinates and width parameters.
#' @param romics_object A romics_object created using romicsCreateObject() that contains 'x' and 'y'
#'        coordinates in its metadata layer
#' @param list_ROI A named list where each element is a numeric vector of length 3: c(X, Y, W).
#'        X = minimum x-coordinate, Y = minimum y-coordinate, W = width of the square ROI.
#'        The maximum coordinates are calculated as X+W and Y+W.
#'        Example: list(ROI1 = c(10, 20, 50), ROI2 = c(100, 150, 30))
#' @param factor_name A character string specifying the name for the new factor to be created
#'        in the metadata layer. Default: "ROIs"
#' @param separator A character string used to separate multiple ROI names when a sample
#'        falls within multiple ROIs. Default: ";"
#' @param outside_label A character string to label samples that don't fall within any ROI.
#'        If NULL, these samples will be left as NA. Default: "Outside"
#' @details This function creates rectangular ROIs and assigns samples to them based on their
#'          x and y coordinates stored in the metadata. Each ROI is defined by:
#'          - X: minimum x-coordinate
#'          - Y: minimum y-coordinate
#'          - W: width (and height) of the square region
#'          - Maximum coordinates: (X+W, Y+W)
#'
#'          Samples can belong to multiple ROIs if their coordinates fall within overlapping regions.
#'          In such cases, ROI names are concatenated using the specified separator.
#'
#'          The function requires 'x' and 'y' factors to exist in the romics_object metadata.
#'          These coordinates are typically present in spatial omics datasets like imaging mass spectrometry,
#'          spatial transcriptomics, or multiplexed immunofluorescence data.
#' @return Returns the romics_object with a new factor added to the metadata layer containing
#'         ROI assignments for each sample. Samples outside all ROIs will be labeled according
#'         to outside_label parameter.
#' @author Geremy Clair
#' @export
romicsCreateFactorFromROIs <- function(romics_object, list_ROI, factor_name = "ROIs",
                                       separator = ";", outside_label = "Outside") {

  arguments <- as.list(match.call())

  # Input validation
  if (!is.romicsObject(romics_object) || missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }

  if (missing(list_ROI) || !is.list(list_ROI) || length(list_ROI) == 0) {
    stop("list_ROI must be a non-empty list of ROI definitions")
  }

  if (!is.character(factor_name) || length(factor_name) != 1 || factor_name == "") {
    stop("factor_name must be a non-empty character string")
  }

  if (!is.character(separator) || length(separator) != 1) {
    stop("separator must be a character string of length 1")
  }

  # Check if ROI list is properly named
  if (is.null(names(list_ROI)) || any(names(list_ROI) == "") || any(is.na(names(list_ROI)))) {
    stop("All elements in list_ROI must have non-empty names")
  }

  # Validate ROI definitions
  for (i in seq_along(list_ROI)) {
    roi <- list_ROI[[i]]
    roi_name <- names(list_ROI)[i]

    if (!is.numeric(roi) || length(roi) != 3) {
      stop(paste("ROI", roi_name, "must be a numeric vector of length 3: c(X, Y, W)"))
    }

    if (any(is.na(roi)) || any(!is.finite(roi))) {
      stop(paste("ROI", roi_name, "contains NA or infinite values"))
    }

    if (roi[3] <= 0) {
      stop(paste("ROI", roi_name, "width (W) must be positive"))
    }
  }

  # Check for required coordinates in metadata
  if (!"x" %in% rownames(romics_object$metadata)) {
    stop("'x' coordinate factor not found in romics_object metadata")
  }

  if (!"y" %in% rownames(romics_object$metadata)) {
    stop("'y' coordinate factor not found in romics_object metadata")
  }

  # Extract metadata and coordinates
  metadata <- romics_object$metadata
  metadata <- data.frame(t(metadata), stringsAsFactors = FALSE)

  # Extract coordinates
  px_x <- as.numeric(metadata$x)
  px_y <- as.numeric(metadata$y)

  # Check for missing coordinates
  if (any(is.na(px_x)) || any(is.na(px_y))) {
    warning("Some samples have missing x or y coordinates and will be excluded from ROI assignment")
  }

  # Initialize ROI assignment matrix
  roi_assignments <- data.frame(matrix(ncol = length(list_ROI), nrow = nrow(metadata)))
  colnames(roi_assignments) <- names(list_ROI)
  rownames(roi_assignments) <- rownames(metadata)

  # Process each ROI
  n_samples_in_rois <- 0
  for (i in seq_along(list_ROI)) {
    roi_name <- names(list_ROI)[i]
    roi <- list_ROI[[i]]

    roi_x_min <- roi[1]
    roi_y_min <- roi[2]
    roi_width <- roi[3]
    roi_x_max <- roi_x_min + roi_width
    roi_y_max <- roi_y_min + roi_width

    # Find samples within this ROI
    within_roi <- px_x >= roi_x_min & px_x <= roi_x_max &
      px_y >= roi_y_min & px_y <= roi_y_max &
      !is.na(px_x) & !is.na(px_y)

    roi_assignments[within_roi, i] <- roi_name
    n_samples_in_rois <- n_samples_in_rois + sum(within_roi, na.rm = TRUE)

    message(paste("ROI", roi_name, ": found", sum(within_roi, na.rm = TRUE), "samples"))
  }

  # Create final ROI factor by combining assignments
  metadata[[factor_name]] <- apply(roi_assignments, 1, function(row) {
    non_empty <- row[nzchar(row) & !is.na(row)]
    if (length(non_empty) == 0) {
      return(if(is.null(outside_label)) NA_character_ else outside_label)
    } else {
      return(paste(non_empty, collapse = separator))
    }
  })

  # Count samples in different categories
  n_outside <- sum(is.na(metadata[[factor_name]]) |
                     (!is.null(outside_label) & metadata[[factor_name]] == outside_label), na.rm = TRUE)
  n_multiple <- sum(grepl(separator, metadata[[factor_name]], fixed = TRUE), na.rm = TRUE)

  # Report summary
  message(paste("ROI assignment summary:"))
  message(paste("- Samples in ROIs:", n_samples_in_rois))
  message(paste("- Samples outside ROIs:", n_outside))
  if (n_multiple > 0) {
    message(paste("- Samples in multiple ROIs:", n_multiple))
  }

  # Check if factor already exists and warn
  if (factor_name %in% rownames(romics_object$metadata)) {
    warning(paste("Factor", factor_name, "already exists and will be replaced"))
    romics_object$metadata <- romics_object$metadata[rownames(romics_object$metadata) != factor_name, ]
  }

  # Convert back to romics format and add new factor
  romics_object$metadata <- data.frame(t(metadata), stringsAsFactors = FALSE)

  message(paste("Factor", factor_name, "was added to the romics metadata"))

  # Update processing steps
  romics_object <- romicsUpdateSteps(romics_object, arguments)

  return(romics_object)
}

#' romicsCreateROIsFromFactor()
#' @description Creates a list of ROI definitions from factor levels in romics_object metadata.
#' This function generates rectangular ROIs that encompass all samples belonging to each level
#' of a specified factor. Each ROI is defined as the minimum bounding rectangle containing
#' all samples of that factor level.
#' @param romics_object A romics_object created using romicsCreateObject() that contains 'x', 'y'
#'        coordinates and the specified factor in its metadata layer
#' @param factor_name A character string specifying the name of the factor to use for ROI creation.
#'        If NULL (default), uses the main_factor of the romics_object.
#' @param padding A numeric value specifying additional padding to add around each ROI boundary.
#'        This ensures ROIs don't touch the exact sample boundaries. Default: 0
#' @param min_roi_size A numeric value specifying the minimum width/height for any ROI.
#'        ROIs smaller than this will be expanded symmetrically. Default: 1
#' @param check_overlap Logical. If TRUE (default), checks whether generated ROIs are fully distinct
#'        and issues warnings for overlapping ROIs.
#' @param square_rois Logical. If TRUE (default), creates square ROIs where width = height = max(width, height).
#'        If FALSE, creates rectangular ROIs with independent width and height.
#' @details This function creates ROIs by:
#'          1. Extracting x,y coordinates for each factor level
#'          2. Computing the bounding rectangle for each level's samples
#'          3. Adding optional padding and enforcing minimum size
#'          4. Checking for overlaps between ROIs if requested
#'
#'          Each ROI is returned as a vector c(X, Y, W) where:
#'          - X: minimum x-coordinate
#'          - Y: minimum y-coordinate
#'          - W: width (and height for square ROIs)
#'
#'          The function requires 'x' and 'y' factors in the metadata and will fail if
#'          any factor level has no samples or samples with missing coordinates.
#'
#'          Overlap checking compares all ROI pairs and reports any overlapping regions,
#'          which might indicate that the factor levels are not spatially distinct.
#' @return Returns a named list where each element is a numeric vector c(X, Y, W) defining
#'         an ROI for the corresponding factor level. Names correspond to factor levels.
#' @examples
#' \dontrun{
#' # Create ROIs from main factor
#' roi_list <- romicsCreateROIsFromFactor(romics_obj)
#'
#' # Create ROIs from specific factor with padding
#' roi_list <- romicsCreateROIsFromFactor(romics_obj,
#'                                        factor_name = "tissue_regions",
#'                                        padding = 5)
#'
#' # Create rectangular ROIs without overlap checking
#' roi_list <- romicsCreateROIsFromFactor(romics_obj,
#'                                        factor_name = "cell_types",
#'                                        square_rois = FALSE,
#'                                        check_overlap = FALSE)
#'
#' # Use the generated ROIs to create a new factor
#' romics_obj <- romicsCreateFactorFromROIs(romics_obj, list_ROI = roi_list)
#' }
#' @seealso \code{\link{romicsCreateFactorFromROIs}} to create factors from ROI definitions
#' @author Geremy Clair
#' @export
romicsCreateROIsFromFactor <- function(romics_object, factor_name = NULL, padding = 0,
                                       min_roi_size = 1, check_overlap = TRUE, square_rois = TRUE) {

  # Input validation
  if (!is.romicsObject(romics_object) || missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }

  if (!is.numeric(padding) || length(padding) != 1 || padding < 0) {
    stop("padding must be a non-negative numeric value")
  }

  if (!is.numeric(min_roi_size) || length(min_roi_size) != 1 || min_roi_size <= 0) {
    stop("min_roi_size must be a positive numeric value")
  }

  if (!is.logical(check_overlap) || length(check_overlap) != 1) {
    stop("check_overlap must be a logical value (TRUE or FALSE)")
  }

  if (!is.logical(square_rois) || length(square_rois) != 1) {
    stop("square_rois must be a logical value (TRUE or FALSE)")
  }

  # Handle factor selection
  if (is.null(factor_name)) {
    factor_name <- romics_object$main_factor
    message(paste("Using main factor:", factor_name))
  } else {
    if (!is.character(factor_name) || length(factor_name) != 1 || factor_name == "") {
      stop("factor_name must be a non-empty character string")
    }
    if (!(factor_name %in% rownames(romics_object$metadata))) {
      stop(paste("Factor", factor_name, "not found in romics_object metadata"))
    }
  }

  # Check for required coordinates in metadata
  if (!"x" %in% rownames(romics_object$metadata)) {
    stop("'x' coordinate factor not found in romics_object metadata")
  }

  if (!"y" %in% rownames(romics_object$metadata)) {
    stop("'y' coordinate factor not found in romics_object metadata")
  }

  # Extract metadata and coordinates
  metadata <- data.frame(t(romics_object$metadata), stringsAsFactors = FALSE)

  # Extract coordinates and factor values
  px_x <- as.numeric(metadata$x)
  px_y <- as.numeric(metadata$y)
  factor_values <- as.character(metadata[[factor_name]])

  # Check for missing data
  missing_coords <- is.na(px_x) | is.na(px_y)
  missing_factor <- is.na(factor_values) | factor_values == ""

  if (any(missing_coords)) {
    warning(paste(sum(missing_coords), "samples have missing x or y coordinates and will be excluded"))
  }

  if (any(missing_factor)) {
    warning(paste(sum(missing_factor), "samples have missing factor values and will be excluded"))
  }

  # Filter out samples with missing data
  valid_samples <- !missing_coords & !missing_factor
  px_x <- px_x[valid_samples]
  px_y <- px_y[valid_samples]
  factor_values <- factor_values[valid_samples]

  if (length(px_x) == 0) {
    stop("No samples with valid coordinates and factor values found")
  }

  # Get unique factor levels
  factor_levels <- unique(factor_values)
  factor_levels <- factor_levels[factor_levels != "" & !is.na(factor_levels)]

  if (length(factor_levels) == 0) {
    stop("No valid factor levels found")
  }

  message(paste("Creating ROIs for", length(factor_levels), "factor levels:", paste(factor_levels, collapse = ", ")))

  # Create ROIs for each factor level
  roi_list <- list()

  for (level in factor_levels) {
    # Get samples for this level
    level_samples <- factor_values == level

    if (sum(level_samples) == 0) {
      warning(paste("No samples found for factor level:", level))
      next
    }

    level_x <- px_x[level_samples]
    level_y <- px_y[level_samples]

    # Calculate bounding box
    x_min <- min(level_x, na.rm = TRUE)
    x_max <- max(level_x, na.rm = TRUE)
    y_min <- min(level_y, na.rm = TRUE)
    y_max <- max(level_y, na.rm = TRUE)

    # Calculate dimensions
    width <- x_max - x_min
    height <- y_max - y_min

    # Handle single point case or enforce minimum size
    if (width < min_roi_size) {
      center_x <- (x_min + x_max) / 2
      x_min <- center_x - min_roi_size / 2
      width <- min_roi_size
    }

    if (height < min_roi_size) {
      center_y <- (y_min + y_max) / 2
      y_min <- center_y - min_roi_size / 2
      height <- min_roi_size
    }

    # Make square if requested
    if (square_rois) {
      max_dim <- max(width, height)
      # Expand the smaller dimension symmetrically
      if (width < max_dim) {
        center_x <- x_min + width / 2
        x_min <- center_x - max_dim / 2
      }
      if (height < max_dim) {
        center_y <- y_min + height / 2
        y_min <- center_y - max_dim / 2
      }
      width <- height <- max_dim
    }

    # Add padding
    x_min <- x_min - padding
    y_min <- y_min - padding
    width <- width + 2 * padding
    if (!square_rois) {
      height <- height + 2 * padding
    }

    # Store ROI (for square ROIs, width = height, so we only store width)
    roi_list[[level]] <- c(x_min, y_min, width)

    message(paste("ROI", level, ": X=", round(x_min, 2), ", Y=", round(y_min, 2),
                  ", W=", round(width, 2), " (", sum(level_samples), " samples)", sep = ""))
  }

  if (length(roi_list) == 0) {
    stop("No valid ROIs could be created")
  }

  # Check for overlaps if requested
  if (check_overlap && length(roi_list) > 1) {
    message("Checking for ROI overlaps...")

    overlaps_found <- FALSE
    roi_names <- names(roi_list)

    for (i in 1:(length(roi_list) - 1)) {
      for (j in (i + 1):length(roi_list)) {
        roi1 <- roi_list[[i]]
        roi2 <- roi_list[[j]]
        name1 <- roi_names[i]
        name2 <- roi_names[j]

        # Extract coordinates
        x1_min <- roi1[1]
        y1_min <- roi1[2]
        x1_max <- x1_min + roi1[3]
        y1_max <- y1_min + roi1[3]

        x2_min <- roi2[1]
        y2_min <- roi2[2]
        x2_max <- x2_min + roi2[3]
        y2_max <- y2_min + roi2[3]

        # Check for overlap
        x_overlap <- max(0, min(x1_max, x2_max) - max(x1_min, x2_min))
        y_overlap <- max(0, min(y1_max, y2_max) - max(y1_min, y2_min))

        if (x_overlap > 0 && y_overlap > 0) {
          overlap_area <- x_overlap * y_overlap
          roi1_area <- roi1[3] * roi1[3]
          roi2_area <- roi2[3] * roi2[3]

          overlap_pct1 <- (overlap_area / roi1_area) * 100
          overlap_pct2 <- (overlap_area / roi2_area) * 100

          warning(paste("ROI overlap detected between", name1, "and", name2,
                        "- Overlap area:", round(overlap_area, 2),
                        paste0("(", round(overlap_pct1, 1), "% of ", name1, ", "),
                        paste0(round(overlap_pct2, 1), "% of ", name2, ")")))
          overlaps_found <- TRUE
        }
      }
    }

    if (!overlaps_found) {
      message("No ROI overlaps detected - all ROIs are fully distinct")
    } else {
      warning("ROI overlaps detected! Factor levels may not be spatially distinct. Consider increasing factor separation or reducing padding.")
    }
  }

  message(paste("Successfully created", length(roi_list), "ROIs from factor", factor_name))

  return(roi_list)
}
