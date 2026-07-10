# Helper function to collapse data by factor level (average/mean by factor)
collapse_data_by_factor <- function(romics_object, data_matrix, factor = "main") {
  if(is.null(factor)) {
    factor <- "main"
  }
  if(factor == "main") {
    factor <- romics_object$main_factor
  }
  if(is.null(factor)) {
    return(data_matrix)
  }

  # Extract factor
  factor_values <- romicsExtractFactor(romics_object, factor = factor)
  factor_levels <- levels(factor_values)

  # Create collapsed matrix
  collapsed_data <- matrix(NA, nrow = nrow(data_matrix), ncol = length(factor_levels),
                          dimnames = list(rownames(data_matrix), as.character(factor_levels)))

  # Collapse each factor level
  for(i in seq_along(factor_levels)) {
    level <- factor_levels[i]
    level_mask <- factor_values == level
    collapsed_data[, i] <- rowMeans(data_matrix[, level_mask, drop = FALSE], na.rm = TRUE)
  }

  # Remove rows that are entirely NaN (all missing values)
  valid_rows <- !apply(is.nan(collapsed_data), 1, all)
  collapsed_data <- collapsed_data[valid_rows, , drop = FALSE]

  return(collapsed_data)
}

#' romicsFilterFeature()
#' @description Filters Romics_objects based on their feature names or statistics note that ALL the filters will be taken in consideration simultaneously
#' @param romics_object A romics_object created using romicsCreateObject()
#' @param ANOVA_filter Either 'none', 'p' or 'padj'. Indicates if an the ANOVA filter has to be used to filter features (only the features below the filter will be retained)
#' @param p Numerical of length 1 indicating the value of the ANOVA_filter cutoff (anything below this value will be conserved).
#' @param feature_names A character vector that enable to filter based on the names of the features, if feature_names set to 'none' the filter won't be applied.
#' @param statCol A column contained in the statistical layer of the romics_object, the list of columns can be obtained by using the function romicsCalculatedStats().
#' @param statCol_filter Character to indicate how this column should be filtered (e.g. '<=0.05','>0.05','==1', '==TRUE', '>2')
#' @param statCol2 A column contained in the statistical layer of the romics_object, the list of columns can be obtained by using the function romicsCalculatedStats().
#' @param statCol2_filter Character to indicate how this column should be filtered (e.g. '<=0.05','>0.05','==1', '==TRUE', '>2')
#' @param mode Either 'keep' or 'drop' to indicate if the features should be kept or dropped based on the filters.
#' @details This function applies multiple filters simultaneously to retain or remove features based on:
#' - ANOVA p-values or adjusted p-values
#' - Feature name patterns
#' - Statistical column criteria
#' All filters are combined using AND logic (all conditions must be met).
#' @return Returns a romics_object with filtered features
#' @author Geremy Clair
#' @export
romicsFilterFeature <- function(romics_object,
                                ANOVA_filter = "none",
                                p = 0.05,
                                feature_names = "none",
                                statCol = "none",
                                statCol_filter = "<=0.05",
                                statCol2 = "none",
                                statCol2_filter = "<=0.05",
                                mode = "keep"){
  arguments <- as.list(match.call())

  # Input validation
  if(!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }
  if(missing(ANOVA_filter)){ANOVA_filter = "none"}
  if(!ANOVA_filter %in% c("p", "padj", "none")){
    stop("ANOVA_filter should be either 'p', 'padj' or 'none'")
  }
  if(missing(p)){p <- 0.05}
  if(!is.numeric(p) || p > 1 || p < 0){
    stop("p should be numeric and comprised between 0 and 1")
  }

  # ANOVA Filter
  if(ANOVA_filter == "none"){
    ANOVA_filter_result <- rep(TRUE, nrow(romics_object$data))
  } else {
    if(ANOVA_filter == "p"){
      if(sum(grepl("ANOVA.*", colnames(romics_object$statistics)) &
             grepl("._p", colnames(romics_object$statistics)) &
             !grepl("._padj", colnames(romics_object$statistics))) == 0){
        warning("The ANOVA_p has not been calculated, no filtering was applied")
        ANOVA_filter_result <- rep(TRUE, nrow(romics_object$data))
      } else {
        ANOVA <- romics_object$statistics[, grepl("ANOVA.*", colnames(romics_object$statistics)) &
                                            grepl("._p", colnames(romics_object$statistics)) &
                                            !grepl("._padj", colnames(romics_object$statistics))]
        ANOVA_filter_result <- !is.na(ANOVA) & ANOVA < p
      }
    } else {
      if(sum(grepl("ANOVA.*", colnames(romics_object$statistics)) &
             grepl("._padj", colnames(romics_object$statistics))) == 0){
        warning("The ANOVA_padj has not been calculated, no filtering was applied")
        ANOVA_filter_result <- rep(TRUE, nrow(romics_object$data))
      } else {
        ANOVA <- romics_object$statistics[, grepl("ANOVA.*", colnames(romics_object$statistics)) &
                                            grepl("._padj", colnames(romics_object$statistics))]
        ANOVA_filter_result <- !is.na(ANOVA) & ANOVA < p
      }
    }
  }

  # Feature names filter validation
  if(missing(feature_names)){feature_names <- "none"}
  if(!is.character(feature_names) || length(feature_names) < 1){
    stop("'feature_names' has to be a character vector of length >= 1")
  }

  # Create feature names filter
  if(feature_names == "none" || feature_names[1] == "none"){
    feature_name_filter <- rep(TRUE, nrow(romics_object$data))
  } else {
    feature_name_filter <- rep(FALSE, nrow(romics_object$data))
    for(i in 1:length(feature_names)){
      feature_name_filter <- feature_name_filter | grepl(feature_names[i], rownames(romics_object$data))
    }
  }

  # Statistical columns filter validation
  if(missing(statCol)){statCol <- "none"}
  if(!is.character(statCol) || length(statCol) != 1){
    stop("statCol should be a character vector of length 1")
  }
  if(missing(statCol2)){statCol2 <- "none"}
  if(!is.character(statCol2) || length(statCol2) != 1){
    stop("statCol2 should be a character vector of length 1")
  }

  # Filter based on statCol
  if(statCol == "none"){
    statCol_filter_result <- rep(TRUE, nrow(romics_object$data))
  } else {
    # Check if the statCol exists
    if(!statCol %in% colnames(romics_object$statistics)){
      cat(paste0("'", statCol, "' is not a column of the statistics layer, below are the usable columns:\n"))
      print(colnames(romics_object$statistics))
      stop("Invalid statCol specified")
    }
    # Check if the statCol_filter exists
    if(missing(statCol_filter)){
      stop("The 'statCol_filter' was missing, the stat column was not filtered")
    } else {
      text <- paste0("romics_object$statistics$`", statCol, "`", statCol_filter)
    }
    # Create the filter
    statCol_filter_result <- eval(parse(text = text))
  }

  # Filter based on statCol2
  if(statCol2 == "none"){
    statCol2_filter_result <- rep(TRUE, nrow(romics_object$data))
  } else {
    # Check if the statCol2 exists
    if(!statCol2 %in% colnames(romics_object$statistics)){
      cat(paste0("'", statCol2, "' is not a column of the statistics layer, below are the usable columns:\n"))
      print(colnames(romics_object$statistics))
      stop("Invalid statCol2 specified")
    }
    # Check if the statCol2_filter exists
    if(missing(statCol2_filter)){
      stop("The 'statCol2_filter' was missing, the stat column was not filtered")
    } else {
      text <- paste0("romics_object$statistics$`", statCol2, "`", statCol2_filter)
    }
    # Create the filter
    statCol2_filter_result <- eval(parse(text = text))
  }

  # Global filter (all conditions must be TRUE)
  filter <- feature_name_filter & ANOVA_filter_result & statCol_filter_result & statCol2_filter_result

  # Check mode
  if(missing(mode)){mode <- "keep"}
  if(!mode %in% c("keep", "drop")){
    stop("'mode' has to be either 'keep' or 'drop'")
  }

  # Reverse the filter if mode == drop
  if(mode == "drop"){
    filter <- !filter
  }

  # Count features before and after filtering
  n_features_before <- nrow(romics_object$data)
  n_features_after <- sum(filter)

  # Apply filter to all relevant layers
  romics_object$data <- romics_object$data[filter, , drop = FALSE]
  romics_object$missingdata <- romics_object$missingdata[filter, , drop = FALSE]

  # Filter statistics layer if it exists
  if(!is.null(romics_object$statistics) && nrow(romics_object$statistics) == n_features_before) {
    romics_object$statistics <- romics_object$statistics[filter, , drop = FALSE]
  }

  # Update steps
  romics_object <- romicsUpdateSteps(romics_object, arguments)

  # Informative message
  message(paste("Feature filtering complete:", n_features_after, "features retained out of", n_features_before))

  return(romics_object)
}

# Helper function to check if a color is valid
is.color <- function(x) {
  tryCatch({
    col2rgb(x)
    TRUE
  }, error = function(e) FALSE)
}

#' romicsHeatmap()
#' @description Plots a scaled heatmap of the data layer from a romics_object using heatmap.2 from the gplots package. Data can be filtered based on the statistics layer. Optionally records cluster assignments to the statistics layer. Can also collapse samples by factor level before plotting.
#' @param romics_object A romics_object created using romicsCreateObject().
#' @param color_palette Character vector of colors. Default: viridis(20)
#' @param color_boundaries Numerical vector of length 2 for min and max of color scale. Default: c(-2,2)
#' @param sample_hclust Boolean. Whether to cluster samples. Default: TRUE
#' @param sample_hclust_method_dist Distance method for samples. Default: "euclidean"
#' @param sample_hclust_method_hclust Agglomeration method for samples. Default: "ward.D"
#' @param feature_hclust Boolean. Whether to cluster features. Default: TRUE
#' @param feature_hclust_method_dist Distance method for features. Default: "euclidean"
#' @param feature_hclust_method_hclust Agglomeration method for features. Default: "ward.D"
#' @param feature_hclust_number Number of clusters for feature coloring. Default: 1
#' @param show_feature_clusters Boolean. Whether to show cluster assignments as a sidebar. Default: FALSE. Requires feature_hclust_number > 1.
#' @param record_feature_clusters Boolean. Whether to save cluster assignments to the statistics layer as 'hclust_clusters' column. Default: FALSE. Requires feature_hclust_number > 1.
#' @param feature_cluster_colors Character vector of colors for cluster sidebar. If NULL, uses a default color palette. Length should match number of unique clusters.
#' @param show_cluster_legend Boolean. Whether to show a legend for cluster colors when show_feature_clusters = TRUE. Default: TRUE.
#' @param cluster_legend_position Position for cluster legend: "topright", "topleft", "bottomright", "bottomleft", "top", "bottom", "left", "right". Default: "topright"
#' @param cluster_legend_cex Size of legend text. Default: 0.8
#' @param collapse_by_factor Character. Name of the factor to collapse samples by (average within each factor level). Default: "none" (no collapsing). If specified, samples are grouped by factor level and averaged before plotting.
#' @param collapse_method Character. Method for collapsing: "mean" or "median". Default: "mean"
#' @param statCol Statistical column for filtering. Default: "none"
#' @param statCol_filter Filter expression for statCol. Default: "<=0.05"
#' @param statCol2 Second statistical column for filtering. Default: "none"
#' @param statCol2_filter Filter expression for statCol2. Default: "<=0.05"
#' @param feature_list Character vector of feature names to include in the heatmap. Default: "all"
#' @param notecol Color for cell note text. Default: "black"
#' @param density.info Density plot on color key: 'histogram', 'density', or 'none'. Default: "none"
#' @param trace Trace lines: "column", "row", "both", or "none". Default: "none"
#' @param labRow Show row labels. Can be TRUE, FALSE, or character vector. Default: FALSE
#' @param cexCol Size of column labels. Default: 1
#' @param margins Margins for labels c(bottom, right). Default: c(15, 5)
#' @param key.title Title for color key. Default: "Scaled Heatmap"
#' @param key.xlab X-axis label for color key. Default: "Z-scores"
#' @param ... Additional parameters for heatmap.2
#' @details When record_feature_clusters = TRUE, the function will save cluster assignments to romics_object$statistics$hclust_clusters
#' and return the updated romics_object. The heatmap is still displayed.
#'
#' When record_feature_clusters = FALSE (default), the function only displays the heatmap and returns NULL invisibly.
#'
#' Note: Cluster assignments are only recorded for features that pass the filtering criteria and appear in the heatmap.
#' Features not in the heatmap will have NA in the hclust_clusters column.
#'
#' When collapse_by_factor is specified, samples are grouped by the specified factor level and averaged
#' (or median-collapsed) before plotting. This reduces the number of columns in the heatmap and allows
#' for easier visualization of group-level patterns. The original romics_object is not modified.
#'
#' @return If record_feature_clusters = TRUE, returns the updated romics_object. Otherwise, returns NULL invisibly. The heatmap is displayed in both cases.
#' @author Geremy Clair and Harsh Bhotika
#' @export
romicsHeatmap <- function(romics_object,
                          color_palette = NULL,
                          color_boundaries = c(-2, 2),
                          sample_hclust = TRUE,
                          sample_hclust_method_dist = "euclidean",
                          sample_hclust_method_hclust = "ward.D",
                          feature_hclust = TRUE,
                          feature_hclust_method_dist = "euclidean",
                          feature_hclust_method_hclust = "ward.D",
                          feature_hclust_number = 1,
                          show_feature_clusters = FALSE,
                          record_feature_clusters = FALSE,
                          feature_cluster_colors = NULL,
                          show_cluster_legend = TRUE,
                          cluster_legend_position = "topright",
                          cluster_legend_cex = 0.8,
                          collapse_by_factor = "none",
                          collapse_method = "mean",
                          statCol = "none",
                          statCol_filter = "<=0.05",
                          statCol2 = "none",
                          statCol2_filter = "<=0.05",
                          feature_list = "all",
                          notecol = "black",
                          density.info = "none",
                          trace = "none",
                          labRow = FALSE,
                          cexCol = 1,
                          margins = c(15, 5),
                          key.title = "Scaled Heatmap",
                          key.xlab = "Z-scores",
                          ...) {
  arguments <- as.list(match.call())
  required_packages <- c("gplots", "viridis", "dendextend")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required. Please install it with: install.packages('", pkg, "')", sep = ""))
    }
  }
  if (!is.romicsObject(romics_object) || missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }
  if (is.null(color_palette)) {
    color_palette <- viridis::viridis(20)
  }
  if (!is.numeric(color_boundaries) || length(color_boundaries) != 2) {
    stop("color_boundaries must be numerical vector of length 2")
  }
  if(missing(feature_list)){feature_list <- "all"}
  if(!is.character(feature_list) && feature_list != "all"){
    warning("'feature_list' should be a character vector or 'all'. 'all' was used by default.")
    feature_list <- "all"
  }
  if(missing(show_feature_clusters)){show_feature_clusters <- FALSE}
  if(!is.logical(show_feature_clusters)){
    warning("show_feature_clusters must be TRUE or FALSE. Setting to FALSE.")
    show_feature_clusters <- FALSE
  }
  if(missing(record_feature_clusters)){record_feature_clusters <- FALSE}
  if(!is.logical(record_feature_clusters)){
    warning("record_feature_clusters must be TRUE or FALSE. Setting to FALSE.")
    record_feature_clusters <- FALSE
  }
  if(missing(show_cluster_legend)){show_cluster_legend <- TRUE}
  if(!is.logical(show_cluster_legend)){
    warning("show_cluster_legend must be TRUE or FALSE. Setting to TRUE.")
    show_cluster_legend <- TRUE
  }
  if(missing(cluster_legend_position)){cluster_legend_position <- "topright"}
  valid_positions <- c("topright", "topleft", "bottomright", "bottomleft", "top", "bottom", "left", "right")
  if(!cluster_legend_position %in% valid_positions){
    warning(paste("Invalid cluster_legend_position. Must be one of:", paste(valid_positions, collapse = ", "), ". Using 'topright'."))
    cluster_legend_position <- "topright"
  }
  if(missing(cluster_legend_cex)){cluster_legend_cex <- 0.8}
  if(!is.numeric(cluster_legend_cex) || cluster_legend_cex <= 0){
    warning("cluster_legend_cex must be a positive number. Using 0.8.")
    cluster_legend_cex <- 0.8
  }
  # ---- COLLAPSE BY FACTOR VALIDATION ----
  if(missing(collapse_by_factor)){collapse_by_factor <- "none"}
  # Handle boolean collapse_by_factor (convert TRUE to main factor)
  if(is.logical(collapse_by_factor)) {
    if(isTRUE(collapse_by_factor)) {
      collapse_by_factor <- romics_object$main_factor
      if(is.null(collapse_by_factor)) collapse_by_factor <- "none"
    } else {
      collapse_by_factor <- "none"
    }
  }
  if(!is.character(collapse_by_factor)){
    warning("collapse_by_factor must be a character string, logical, or NULL. Using 'none'.")
    collapse_by_factor <- "none"
  }
  if(missing(collapse_method)){collapse_method <- "mean"}
  if(!collapse_method %in% c("mean", "median")){
    warning("collapse_method must be 'mean' or 'median'. Using 'mean'.")
    collapse_method <- "mean"
  }
  collapse_enabled <- collapse_by_factor != "none"
  if(collapse_enabled) {
    # Validate factor exists
    factor_names <- try(romicsFactorNames(romics_object), silent = TRUE)
    if(inherits(factor_names, "try-error")) {
      stop("Could not extract factor names from romics_object")
    }
    if(!collapse_by_factor %in% factor_names) {
      stop("Factor '", collapse_by_factor, "' not found in metadata. Available factors: ",
           paste(factor_names, collapse = ", "))
    }
    message("Collapsing samples by factor: ", collapse_by_factor, " (method: ", collapse_method, ")")
  }
  if((show_feature_clusters || record_feature_clusters) && feature_hclust_number <= 1) {
    warning("show_feature_clusters or record_feature_clusters requires feature_hclust_number > 1. These options will be ignored.")
    show_feature_clusters <- FALSE
    record_feature_clusters <- FALSE
  }
  valid_dist_methods <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
  valid_hclust_methods <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
  if (!sample_hclust_method_dist %in% valid_dist_methods) {
    stop("Invalid sample_hclust_method_dist. Must be one of: ", paste(valid_dist_methods, collapse = ", "))
  }
  if (!feature_hclust_method_dist %in% valid_dist_methods) {
    stop("Invalid feature_hclust_method_dist. Must be one of: ", paste(valid_dist_methods, collapse = ", "))
  }
  if (!sample_hclust_method_hclust %in% valid_hclust_methods) {
    stop("Invalid sample_hclust_method_hclust. Must be one of: ", paste(valid_hclust_methods, collapse = ", "))
  }
  if (!feature_hclust_method_hclust %in% valid_hclust_methods) {
    stop("Invalid feature_hclust_method_hclust. Must be one of: ", paste(valid_hclust_methods, collapse = ", "))
  }
  if (!is.numeric(feature_hclust_number) || feature_hclust_number <= 0) {
    stop("feature_hclust_number should be a positive integer")
  }
  # ---- COLLAPSE SAMPLES IF REQUESTED ----
  data_to_plot <- romics_object$data
  if(collapse_enabled) {
    message("Filtering data...")
    filtered_object <- romicsFilterFeature(romics_object,
                                           ANOVA_filter = "none",
                                           p = 0.05,
                                           feature_names = "none",
                                           statCol = statCol,
                                           statCol_filter = statCol_filter,
                                           statCol2 = statCol2,
                                           statCol2_filter = statCol2_filter)
    # Extract factor
    factor_data <- romicsExtractFactor(filtered_object, factor = collapse_by_factor)
    factor_levels <- levels(factor_data)

    message("Collapsing ", ncol(filtered_object$data), " samples to ",
            length(factor_levels), " factor levels using '", collapse_method, "' method...")

    # Create collapsed data matrix
    collapsed_data <- matrix(NA,
                             nrow = nrow(filtered_object$data),
                             ncol = length(factor_levels),
                             dimnames = list(rownames(filtered_object$data),
                                             as.character(factor_levels)))

    # Collapse each factor level
    for (i in seq_along(factor_levels)) {
      level <- factor_levels[i]
      level_mask <- factor_data == level

      if (collapse_method == "mean") {
        collapsed_data[, i] <- rowMeans(filtered_object$data[, level_mask, drop = FALSE], na.rm = TRUE)
      } else if (collapse_method == "median") {
        collapsed_data[, i] <- apply(filtered_object$data[, level_mask, drop = FALSE], 1,
                                     median, na.rm = TRUE)
      }
    }

    data_to_plot <- collapsed_data
    message("Collapsed data dimensions: ", nrow(data_to_plot), " features x ",
            ncol(data_to_plot), " factor levels")
  } else {
    message("Filtering data...")
    filtered_object <- romicsFilterFeature(romics_object,
                                           ANOVA_filter = "none",
                                           p = 0.05,
                                           feature_names = "none",
                                           statCol = statCol,
                                           statCol_filter = statCol_filter,
                                           statCol2 = statCol2,
                                           statCol2_filter = statCol2_filter)
    data_to_plot <- filtered_object$data
  }
  # ---- FEATURE FILTERING ----
  if(!identical(feature_list, "all")) {
    available_features <- feature_list[feature_list %in% rownames(data_to_plot)]
    missing_features <- feature_list[!feature_list %in% rownames(data_to_plot)]
    if(length(missing_features) > 0) {
      warning(paste("The following features from feature_list were not found in the data:",
                    paste(missing_features, collapse = ", ")))
    }
    if(length(available_features) == 0) {
      stop("None of the features in feature_list were found in the data after statistical filtering.")
    }
    data_to_plot <- data_to_plot[available_features, , drop = FALSE]
    message(paste("Filtered to", length(available_features), "features from feature_list"))
  }
  if (nrow(data_to_plot) == 0) {
    stop("No features remain after filtering. Please adjust your filter criteria or feature_list.")
  }
  message(paste("Plotting heatmap with", nrow(data_to_plot), "features and", ncol(data_to_plot), "samples"))
  feature_hclust_number <- as.integer(feature_hclust_number)
  # ---- SCALING ----
  message("Scaling data...")
  scaled_data <- t(scale(t(data_to_plot), center = TRUE, scale = TRUE))
  scaled_data_for_clustering <- scaled_data
  scaled_data_for_clustering[is.na(scaled_data_for_clustering)] <- -20

  # ---- SAMPLE CLUSTERING ----
  col_dendrogram <- FALSE
  if (sample_hclust && ncol(data_to_plot) > 1) {
    message("Performing sample clustering...")

    # If collapsed by factor, create colors based on factor levels
    if(collapse_enabled) {
      # Get unique factor levels and assign colors
      factor_data <- romicsExtractFactor(filtered_object, factor = collapse_by_factor)
      factor_levels <- levels(factor_data)

      # Create a color palette for factor levels
      n_levels <- length(factor_levels)
      if (n_levels <= 8) {
        level_colors <- RColorBrewer::brewer.pal(max(3, n_levels), "Set2")
      } else {
        level_colors <- grDevices::rainbow(n_levels)
      }

      # Map levels to colors
      level_to_color <- setNames(level_colors[1:n_levels], factor_levels)
      colors_dend <- level_to_color[as.character(factor_levels)]
    } else {
      # Original behavior: use metadata if available
      if (!is.null(romics_object$metadata) && nrow(romics_object$metadata) > 0) {
        colors_dend <- as.character(romics_object$metadata[1, ])
      } else {
        colors_dend <- rep("gray", ncol(data_to_plot))
      }
    }

    distance_samples <- dist(t(data_to_plot), method = sample_hclust_method_dist)
    hc_samples <- hclust(distance_samples, method = sample_hclust_method_hclust)
    sample_dd <- as.dendrogram(hc_samples)

    if(collapse_enabled) {
      # Map colors directly to collapsed levels
      colors_dend <- level_to_color[colnames(data_to_plot)]
    } else {
      colors_dend <- colors_dend[hc_samples$order]
    }
    col_dendrogram <- dendextend::color_branches(sample_dd, col = as.vector(colors_dend))
  }

  # ---- FEATURE CLUSTERING ----
  row_dendrogram <- FALSE
  heatmap_clusters <- NULL
  cluster_colors_for_legend <- NULL
  if (feature_hclust && nrow(data_to_plot) > 1) {
    message("Performing feature clustering...")
    dist_feature <- dist(scaled_data_for_clustering, method = feature_hclust_method_dist)
    feature_hc <- hclust(dist_feature, method = feature_hclust_method_hclust)
    feature_dd <- as.dendrogram(feature_hc)
    if (feature_hclust_number > 1) {
      feature_dd <- dendextend::color_branches(feature_dd, k = feature_hclust_number)
      feature_dd <- dendextend::color_labels(feature_dd, k = feature_hclust_number)
      heatmap_clusters <- cutree(feature_hc, k = feature_hclust_number)
      names(heatmap_clusters) <- rownames(scaled_data)
      message(paste("Features assigned to", feature_hclust_number, "clusters"))
      cluster_counts <- table(heatmap_clusters)
      for(i in 1:length(cluster_counts)) {
        message(paste("  Cluster", i, ":", cluster_counts[i], "features"))
      }
    }
    row_dendrogram <- feature_dd
  }

  # ---- CLUSTER SIDEBAR ----
  row_side_colors <- NULL
  if(!is.null(heatmap_clusters)) {
    message("Preparing cluster sidebar...")
    unique_clusters <- sort(unique(heatmap_clusters))
    n_clusters <- length(unique_clusters)
    if(is.null(feature_cluster_colors)) {
      if(requireNamespace("RColorBrewer", quietly = TRUE)) {
        if(n_clusters <= 12) {
          feature_cluster_colors <- RColorBrewer::brewer.pal(min(max(3, n_clusters), 12), "Set3")[1:n_clusters]
        } else {
          feature_cluster_colors <- rainbow(n_clusters)
        }
      } else {
        feature_cluster_colors <- rainbow(n_clusters)
      }
    } else {
      if(length(feature_cluster_colors) < n_clusters) {
        warning(paste("Not enough colors provided for", n_clusters, "clusters. Using default palette."))
        feature_cluster_colors <- rainbow(n_clusters)
      }
    }

    # Create color mapping for clusters
    cluster_color_map <- setNames(feature_cluster_colors[1:n_clusters], unique_clusters)
    row_side_colors <- as.character(cluster_color_map[as.character(heatmap_clusters)])

    # Only set up legend if show_feature_clusters is TRUE
    if(show_feature_clusters) {
      cluster_colors_for_legend <- feature_cluster_colors[1:n_clusters]
      names(cluster_colors_for_legend) <- paste("Cluster", unique_clusters)
    }
  }

  # ---- RECORD CLUSTERS (ONLY IF NOT COLLAPSED) ----
  if(record_feature_clusters && !is.null(heatmap_clusters) && !collapse_enabled) {
    message("Recording cluster assignments to statistics layer...")
    reference_features <- rownames(romics_object$data)
    cluster_assignments <- rep(NA, length(reference_features))
    names(cluster_assignments) <- reference_features
    cluster_assignments[names(heatmap_clusters)] <- heatmap_clusters
    if("hclust_clusters" %in% colnames(romics_object$statistics)){
      romics_object$statistics$hclust_clusters <- NULL
      warning("The hclust_clusters column was previously present in statistics. It has been removed and replaced with the new clustering results.")
    }
    if(nrow(romics_object$statistics) != nrow(romics_object$data)) {
      message(paste("Statistics layer has", nrow(romics_object$statistics),
                    "rows but data layer has", nrow(romics_object$data), "rows."))
      message("Aligning statistics layer to match data layer...")
      common_features <- intersect(rownames(romics_object$statistics), reference_features)
      romics_object$statistics <- romics_object$statistics[common_features, , drop = FALSE]
      romics_object$statistics <- romics_object$statistics[reference_features[reference_features %in% common_features], , drop = FALSE]
      message(paste("Statistics layer now has", nrow(romics_object$statistics), "rows"))
    }
    romics_object$statistics$hclust_clusters <- cluster_assignments[rownames(romics_object$statistics)]
    n_clustered <- sum(!is.na(romics_object$statistics$hclust_clusters))
    n_total <- nrow(romics_object$statistics)
    message(paste("Cluster assignments recorded:", n_clustered, "features clustered,",
                  n_total - n_clustered, "features not in heatmap (set to NA)"))
    romics_object <- romicsUpdateSteps(romics_object, arguments)
  } else if(record_feature_clusters && collapse_enabled) {
    warning("record_feature_clusters is ignored when collapse_by_factor is enabled.")
  }

  # ---- HEATMAP RENDERING ----
  palette <- colorRampPalette(color_palette)(n = 299)
  col_breaks <- seq(color_boundaries[1], color_boundaries[2], length = 300)

  row_labels <- labRow
  if (is.logical(labRow)) {
    if (labRow) {
      row_labels <- rownames(scaled_data)
    } else {
      row_labels <- FALSE
    }
  } else if (is.character(labRow)) {
    if (length(labRow) != nrow(scaled_data)) {
      warning("Length of labRow doesn't match number of rows. Using rownames instead.")
      row_labels <- rownames(scaled_data)
    }
  }

  message("Generating heatmap...")
  extra_args <- list(...)
  legacy_params <- c("stat_sidebar", "stat_sidebar_columns", "stat_sidebar_tests",
                     "stat_sidebar_padj_type", "stat_sidebar_threshold",
                     "stat_sidebar_colors", "stat_sidebar_max_tests",
                     "stat_sidebar_method", "RowSideColor", "ANOVA_filter", "p")
  extra_args <- extra_args[!names(extra_args) %in% legacy_params]

  heatmap_args <- list(
    x = scaled_data,
    notecol = notecol,
    density.info = density.info,
    trace = trace,
    col = palette,
    breaks = col_breaks,
    Rowv = row_dendrogram,
    Colv = col_dendrogram,
    labRow = row_labels,
    margins = margins,
    cexCol = cexCol,
    keysize = 1,
    key.title = key.title,
    key.xlab = key.xlab
  )

  if(!is.null(row_side_colors) && length(row_side_colors) > 0) {
    # Validate colors before adding to heatmap
    valid_colors <- all(sapply(row_side_colors, function(x) {
      tryCatch(is.color(x), error = function(e) FALSE)
    }))
    if(valid_colors) {
      heatmap_args$RowSideColors <- row_side_colors
    } else {
      warning("Invalid colors in row_side_colors. Cluster sidebar will not be displayed.")
    }
  }

  heatmap_args <- c(heatmap_args, extra_args)
  do.call(gplots::heatmap.2, heatmap_args)

  if(show_feature_clusters && show_cluster_legend && !is.null(cluster_colors_for_legend)) {
    legend(cluster_legend_position,
           legend = names(cluster_colors_for_legend),
           fill = cluster_colors_for_legend,
           border = "black",
           bty = "n",
           cex = cluster_legend_cex,
           title = "Feature Clusters")
  }

  message("Heatmap complete!")

  if(record_feature_clusters && !collapse_enabled) {
    message("Returning updated romics_object with cluster assignments")
    return(romics_object)
  } else {
    invisible(NULL)
  }
}


#' romicsComplexHeatmap()
#' @description Creates a heatmap with statistical annotations from a romics object
#' @param romics_object A romics object containing $data, $statistics, and $metadata layers
#' @param comparison_prefix Character string of the comparison prefix (e.g., "Alveoli_pleura_vs_Bronchi_bronchiole"). If NULL, auto-detects from statistics columns
#' @param factor_name Factor name from metadata for sample grouping. If NULL, uses main factor. Use romicsFactorNames() to see available factors
#' @param feature_list Character vector of features to include in heatmap. If provided, creates subset using romicsFeatureSubset(). If NULL, uses all features
#' @param p_threshold Numeric threshold for p-values (default 0.05)
#' @param use_adjusted_p TRUE or FALSE, if TRUE uses adjusted p-values, if FALSE uses raw p-values
#' @param scale_data TRUE or FALSE, if TRUE scales the data by rows (z-score)
#' @param cluster_rows TRUE or FALSE, if TRUE performs hierarchical clustering of rows
#' @param cluster_cols TRUE or FALSE, if TRUE performs hierarchical clustering of columns
#' @param show_row_names TRUE or FALSE, if TRUE shows feature names on rows
#' @param show_column_names TRUE or FALSE, if TRUE shows sample names on columns
#' @param heatmap_colors Character vector of colors for the heatmap. Default: viridis(100)
#' @param color_breaks Numeric vector of breaks for color mapping. If NULL, uses automatic breaks
#' @param allow_missing TRUE or FALSE, if TRUE allows missing values in the heatmap
#' @param missing_color Color for missing values in the heatmap (default: "grey")
#' @param significance_color Color for significant results in all test annotations (default: "red")
#' @details Creates a heatmap with side annotations automatically detected from statistics columns (T-test, Wilcoxon, GLM binomial, ANOVA) and sample groups from metadata. If feature_list is provided, automatically subsets the romics object first.
#' @return ComplexHeatmap object
#' @author Geremy Clair
#' @export
romicsComplexHeatmap <- function(romics_object,
                                 comparison_prefix = NULL,
                                 factor_name = NULL,
                                 feature_list = NULL,
                                 p_threshold = 0.05,
                                 use_adjusted_p = TRUE,
                                 scale_data = TRUE,
                                 cluster_rows = TRUE,
                                 cluster_cols = TRUE,
                                 show_row_names = FALSE,
                                 show_column_names = TRUE,
                                 heatmap_colors = NULL,
                                 color_breaks = NULL,
                                 allow_missing = TRUE,
                                 missing_color = "grey",
                                 significance_color = "red",
                                 collapse_by_factor = NULL) {

  # Check for required packages
  required_packages <- c("ComplexHeatmap", "circlize", "RColorBrewer")
  missing_packages <- c()
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing_packages <- c(missing_packages, pkg)
    }
  }
  if (length(missing_packages) > 0) {
    stop("Required package(s) not installed: ", paste(missing_packages, collapse = ", "),
         "\nPlease install with: BiocManager::install(c('", paste(missing_packages, collapse = "', '"), "'))")
  }

  # Load required libraries
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)

  # Set default palette to viridis if not specified
  if (is.null(heatmap_colors)) {
    heatmap_colors <- viridis::viridis(100)
  }

  # Create feature subset if feature_list is provided
  if (!is.null(feature_list)) {
    if (!is.character(feature_list)) {
      stop("feature_list must be a character vector")
    }
    message("Creating subset with ", length(feature_list), " features...")
    romics_object <- romicsFeatureSubset(romics_object, feature_list = feature_list)
    message("Subset created with ", nrow(romics_object$data), " features")
    if (nrow(romics_object$data) == 0) {
      stop("No features found after subsetting. Check that feature_list matches rownames in romics_object$data")
    }
  }

  # Extract data and statistics
  data_matrix <- as.matrix(romics_object$data)
  stats_df <- romics_object$statistics
  stat_cols <- colnames(stats_df)

  # Auto-detect all comparison prefixes if not provided
  if (is.null(comparison_prefix)) {
    # Define test patterns to look for
    test_patterns <- c("Ttest_p", "Ttest_padj",
                       "Wilcoxon_p", "Wilcoxon_padj",
                       "glmBinomialTest_p", "glmBinomialTest_padj",
                       "ANOVA_p", "ANOVA_padj")

    # Find all columns ending with test patterns
    possible_prefixes <- c()
    for (pattern in test_patterns) {
      matching_cols <- stat_cols[grepl(paste0("_", pattern, "$"), stat_cols)]
      if (length(matching_cols) > 0) {
        # Extract prefix by removing the test pattern
        prefixes <- gsub(paste0("_", pattern, "$"), "", matching_cols)
        possible_prefixes <- c(possible_prefixes, prefixes)
      }
    }

    # Get unique prefixes and store all for later use
    all_comparisons <- unique(possible_prefixes)

    if (length(all_comparisons) == 0) {
      stop("Could not auto-detect comparison prefix. Please specify comparison_prefix parameter.")
    }

    message("Multiple comparisons found: ", paste(all_comparisons, collapse = ", "))
  } else {
    # User specified comparison_prefix - can be single value or vector
    all_comparisons <- as.character(comparison_prefix)
  }

  message("Using comparisons: ", paste(all_comparisons, collapse = ", "))

  # Dynamically find all available test columns for each comparison
  # Pattern: {comparison_prefix}_{test_name}_p or {comparison_prefix}_{test_name}_padj
  # Look for both _p and _padj, prioritize based on use_adjusted_p
  p_suffix_primary <- if (use_adjusted_p) "_padj" else "_p"
  p_suffix_secondary <- if (use_adjusted_p) "_p" else "_padj"

  all_available_tests <- list()
  for (comp_prefix in all_comparisons) {
    # Find all columns matching: comp_prefix_*_p or comp_prefix_*_padj
    pattern <- paste0("^", gsub("([.|()\\^{}+*$?]|\\[|\\])", "\\\\\\1", comp_prefix), "_.*(_p|_padj)$")
    matching_cols <- stat_cols[grepl(pattern, stat_cols)]

    available_tests <- list()
    test_name_map <- list()  # Map test_name to list of available suffixes

    if (length(matching_cols) > 0) {
      for (p_col in matching_cols) {
        # Extract test name from column: comp_prefix_testname_p -> testname
        test_name <- gsub(paste0("^", gsub("([.|()\\^{}+*$?]|\\[|\\])", "\\\\\\1", comp_prefix), "_"), "", p_col)
        suffix <- if (grepl("_padj$", test_name)) "_padj" else "_p"
        test_name <- gsub("(_p|_padj)$", "", test_name)

        # Store which suffix(es) are available for this test
        if (!(test_name %in% names(test_name_map))) {
          test_name_map[[test_name]] <- list()
        }
        test_name_map[[test_name]][[suffix]] <- p_col
      }

      # Now prioritize: use primary suffix if available, fall back to secondary
      for (test_name in names(test_name_map)) {
        if (p_suffix_primary %in% names(test_name_map[[test_name]])) {
          available_tests[[test_name]] <- test_name_map[[test_name]][[p_suffix_primary]]
        } else if (p_suffix_secondary %in% names(test_name_map[[test_name]])) {
          available_tests[[test_name]] <- test_name_map[[test_name]][[p_suffix_secondary]]
        }
      }
    }

    if (length(available_tests) > 0) {
      all_available_tests[[comp_prefix]] <- available_tests
    }
  }

  if (length(all_available_tests) == 0) {
    stop("No statistical test columns found for any comparisons",
         "\nLooking for columns matching: prefix_testname_p or prefix_testname_padj",
         "\nExample: ALV_vs_others_Ttest_p, ALV_vs_others_bobtest_p")
  }

  # Report all detected tests
  for (comp_prefix in names(all_available_tests)) {
    message("  ", comp_prefix, ": ", paste(names(all_available_tests[[comp_prefix]]), collapse = ", "))
  }

  # Handle missing values
  n_na <- sum(is.na(data_matrix))
  n_inf <- sum(is.infinite(data_matrix))
  n_nan <- sum(is.nan(data_matrix))
  total_missing <- n_na + n_inf + n_nan

  if (total_missing > 0) {
    message("Found problematic values in data:")
    message("  NA values: ", n_na)
    message("  Infinite values: ", n_inf)
    message("  NaN values: ", n_nan)
    message("  Total missing: ", total_missing, " out of ", length(data_matrix),
            " values (", round(100 * total_missing / length(data_matrix), 2), "%)")
  }

  # Store original missing pattern for display
  original_na_pattern <- is.na(data_matrix) | is.infinite(data_matrix) | is.nan(data_matrix)

  # Handle infinite/NaN values (convert to NA for proper handling)
  data_matrix[is.infinite(data_matrix) | is.nan(data_matrix)] <- NA

  # Create a matrix for clustering (imputed) and one for display (with NAs)
  display_matrix <- data_matrix
  clustering_matrix <- data_matrix

  if (allow_missing) {
    message("Preserving missing values for display...")
    if (cluster_rows || cluster_cols) {
      # Replace NA/Inf/NaN with row means in CLUSTERING matrix only
      for (i in 1:nrow(clustering_matrix)) {
        row_vals <- clustering_matrix[i, ]
        finite_vals <- row_vals[is.finite(row_vals)]
        if (length(finite_vals) > 0) {
          row_mean <- mean(finite_vals, na.rm = TRUE)
          clustering_matrix[i, !is.finite(row_vals)] <- row_mean
        } else {
          overall_mean <- mean(clustering_matrix[is.finite(clustering_matrix)], na.rm = TRUE)
          if (is.finite(overall_mean)) {
            clustering_matrix[i, ] <- overall_mean
          } else {
            clustering_matrix[i, ] <- 0
          }
        }
      }
      message("Applied row-wise imputation for clustering (display matrix preserves NAs)")
    }
    data_matrix_for_heatmap <- display_matrix
    data_matrix_for_clustering <- clustering_matrix
  } else {
    complete_rows <- complete.cases(data_matrix)
    if (sum(complete_rows) == 0) {
      stop("No complete rows found. Set allow_missing=TRUE to handle missing values.")
    }
    data_matrix <- data_matrix[complete_rows, , drop = FALSE]
    stats_df <- stats_df[complete_rows, , drop = FALSE]
    data_matrix_for_heatmap <- data_matrix
    data_matrix_for_clustering <- data_matrix
    message("Removed ", sum(!complete_rows), " rows with missing values")
  }

  # Scale data if requested
  if (scale_data) {
    row_sds <- apply(data_matrix_for_clustering, 1, sd, na.rm = TRUE)
    zero_var_rows <- row_sds == 0 | is.na(row_sds)

    if (any(zero_var_rows)) {
      message("Found ", sum(zero_var_rows), " rows with zero variance. Removing for scaling...")
      data_matrix_for_clustering <- data_matrix_for_clustering[!zero_var_rows, , drop = FALSE]
      data_matrix_for_heatmap <- data_matrix_for_heatmap[!zero_var_rows, , drop = FALSE]
      stats_df <- stats_df[!zero_var_rows, , drop = FALSE]
    }

    if (nrow(data_matrix_for_clustering) == 0) {
      stop("No rows with variance remaining for scaling.")
    }

    scaled_clustering <- t(scale(t(data_matrix_for_clustering)))
    scaled_clustering[is.nan(scaled_clustering)] <- 0

    scaled_display <- t(scale(t(data_matrix_for_heatmap)))

    data_matrix_for_heatmap <- scaled_display
    data_matrix_for_clustering <- scaled_clustering
    heatmap_name <- "Z-score"
  } else {
    heatmap_name <- "Abundance"
  }

  # Use the display matrix for the heatmap
  data_matrix <- data_matrix_for_heatmap

  # Collapse by factor if requested
  if (!is.null(collapse_by_factor) && collapse_by_factor != FALSE) {
    collapse_factor <- if(isTRUE(collapse_by_factor)) factor_name else collapse_by_factor
    # Create temporary romics object with filtered data for collapsing (consistent with current data state)
    temp_romics <- romics_object
    temp_romics$data <- as.matrix(romics_object$data[rownames(data_matrix), , drop = FALSE])
    temp_romics$missingdata <- romics_object$missingdata[rownames(data_matrix), colnames(data_matrix), drop = FALSE]
    temp_romics$metadata <- romics_object$metadata[, colnames(data_matrix), drop = FALSE]
    data_matrix <- collapse_data_by_factor(temp_romics, data_matrix, collapse_factor)
    data_matrix_for_clustering <- collapse_data_by_factor(temp_romics, data_matrix_for_clustering, collapse_factor)
  }

  # Ensure stats_df rows match data_matrix rows (in case some were removed)
  if (!all(rownames(data_matrix) %in% rownames(stats_df))) {
    message("Warning: Some rows in data_matrix not found in stats_df. Using intersection.")
  }
  common_rows <- intersect(rownames(data_matrix), rownames(stats_df))
  stats_df <- stats_df[common_rows, , drop = FALSE]
  data_matrix <- data_matrix[common_rows, , drop = FALSE]
  data_matrix_for_clustering <- data_matrix_for_clustering[common_rows, , drop = FALSE]

  # Check if we still have data
  if (nrow(data_matrix) == 0) {
    stop("No valid data remaining.")
  }

  # Extract sample groups from metadata
  if (is.null(factor_name)) {
    sample_groups <- romicsExtractFactor(romics_object, factor = "main")
  } else {
    sample_groups <- romicsExtractFactor(romics_object, factor = factor_name)
  }

  # After collapsing, update sample_groups and col_annotation to match new data_matrix columns
  # When collapsed, columns are factor levels, so create a named factor from them
  was_collapsed <- FALSE
  if (!all(colnames(data_matrix) %in% names(sample_groups))) {
    # Data was collapsed - columns are now factor levels
    was_collapsed <- TRUE
    sample_groups <- factor(colnames(data_matrix), levels = unique(colnames(data_matrix)))
    names(sample_groups) <- colnames(data_matrix)
  } else {
    # Not collapsed - subset as before
    sample_groups <- sample_groups[colnames(data_matrix)]
  }
  group_levels <- levels(sample_groups)

  # Recreate column annotation if needed
  col_annotation <- NULL
  if (length(group_levels) > 1) {
    if ("colors_romics" %in% rownames(romics_object$metadata)) {
      tryCatch({
        # Get the original factor used for collapsing
        original_factor <- if(was_collapsed) collapse_factor else factor_name
        if(is.null(original_factor)) {
          original_factor <- romics_object$main_factor
        }

        # Get mapping of original samples to factor levels
        original_factor_values <- romicsExtractFactor(romics_object, factor = original_factor)

        # Get all original sample colors
        original_metadata_cols <- colnames(romics_object$metadata)
        all_sample_colors <- romics_object$metadata["colors_romics", original_metadata_cols, drop = TRUE]

        # For each group level (which is now a column in data_matrix after collapsing)
        # find which original samples belong to it and get a representative color
        group_colors <- c()
        for (group in group_levels) {
          # Find original samples in this group
          group_samples <- names(original_factor_values)[original_factor_values == group]
          group_samples_in_metadata <- intersect(group_samples, names(all_sample_colors))
          if (length(group_samples_in_metadata) > 0) {
            group_color <- as.character(all_sample_colors[group_samples_in_metadata[1]])
            group_colors[group] <- group_color
          }
        }

        if (length(group_colors) == length(group_levels)) {
          message("Recreated column annotation with colors: ", paste(paste0(names(group_colors), " = ", group_colors), collapse = ", "))
        } else {
          stop("Missing colors for some groups: ", paste(group_levels[!(group_levels %in% names(group_colors))], collapse = ", "))
        }
      }, error = function(e) {
        message("Error extracting colors from metadata: ", e$message)
        message("Using default colors instead")
        group_colors <<- brewer.pal(min(max(length(group_levels), 3), 8), "Set2")[1:length(group_levels)]
        names(group_colors) <<- group_levels
      })
    } else {
      message("colors_romics row not found in metadata, using default colors")
      group_colors <- brewer.pal(min(max(length(group_levels), 3), 8), "Set2")[1:length(group_levels)]
      names(group_colors) <- group_levels
    }

    col_annotation <- HeatmapAnnotation(
      Group = sample_groups,
      col = list(Group = group_colors),
      annotation_name_gp = gpar(fontsize = 10)
    )
  }

  # Create significance annotations for all comparisons and tests
  row_annotation <- NULL

  annotation_components <- list()
  color_components <- list()

  for (comp_prefix in names(all_available_tests)) {
    available_tests <- all_available_tests[[comp_prefix]]

    for (test_name in names(available_tests)) {
      p_col <- available_tests[[test_name]]

      if (p_col %in% colnames(stats_df)) {
        sig_values <- factor(
          ifelse(stats_df[[p_col]] < p_threshold, "Significant", "Not Sig."),
          levels = c("Significant", "Not Sig.")
        )
        # Create unique name for each comparison-test combination
        annotation_name <- paste0(comp_prefix, " (", test_name, ")")
        annotation_components[[annotation_name]] <- sig_values
        color_components[[annotation_name]] <- c("Significant" = significance_color, "Not Sig." = "lightgrey")
      }
    }
  }

  if (length(annotation_components) > 0) {
    args_list <- annotation_components
    args_list$col <- color_components
    args_list$annotation_name_gp <- gpar(fontsize = 8)
    row_annotation <- do.call(rowAnnotation, args_list)
  }


  # Create color scheme for heatmap
  if (is.null(color_breaks)) {
    finite_vals <- data_matrix[is.finite(data_matrix)]
    if (scale_data) {
      breaks <- seq(-2, 2, length.out = length(heatmap_colors))
    } else {
      min_val <- min(finite_vals, na.rm = TRUE)
      max_val <- max(finite_vals, na.rm = TRUE)
      breaks <- seq(min_val, max_val, length.out = length(heatmap_colors))
    }
  } else {
    breaks <- color_breaks
  }

  # Handle different numbers of colors - if user provided explicit breaks, adjust if needed
  if (length(breaks) != length(heatmap_colors)) {
    if (length(heatmap_colors) == 2 && length(breaks) == 3) {
      breaks <- c(breaks[1], breaks[3])
    } else {
      # Interpolate breaks to match number of colors
      breaks <- seq(min(breaks), max(breaks), length.out = length(heatmap_colors))
    }
  }

  col_fun <- colorRamp2(breaks, heatmap_colors)

  # Create clustering objects if needed
  if (cluster_rows || cluster_cols) {
    if (cluster_rows) {
      row_clust <- hclust(dist(data_matrix_for_clustering, method = "euclidean"), method = "complete")
    } else {
      row_clust <- FALSE
    }

    if (cluster_cols) {
      col_clust <- hclust(dist(t(data_matrix_for_clustering), method = "euclidean"), method = "complete")
    } else {
      col_clust <- FALSE
    }
  } else {
    row_clust <- FALSE
    col_clust <- FALSE
  }

  # Create the main heatmap with error handling
  ht <- tryCatch({
    Heatmap(
      data_matrix,
      name = heatmap_name,
      col = col_fun,
      na_col = missing_color,
      # Clustering - use precomputed clustering objects
      cluster_rows = row_clust,
      cluster_columns = col_clust,
      # Annotations
      right_annotation = row_annotation,
      top_annotation = col_annotation,
      # Names and labels
      show_row_names = show_row_names,
      show_column_names = show_column_names,
      row_names_gp = gpar(fontsize = 8),
      column_names_gp = gpar(fontsize = 10),
      # Heatmap appearance
      border = TRUE,
      heatmap_legend_param = list(
        title_gp = gpar(fontsize = 12),
        labels_gp = gpar(fontsize = 10)
      )
    )
  }, error = function(e) {
    message("Clustering failed. Creating heatmap without clustering...")
    Heatmap(
      data_matrix,
      name = heatmap_name,
      col = col_fun,
      na_col = missing_color,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      right_annotation = row_annotation,
      top_annotation = col_annotation,
      show_row_names = show_row_names,
      show_column_names = show_column_names,
      row_names_gp = gpar(fontsize = 8),
      column_names_gp = gpar(fontsize = 10),
      border = TRUE,
      heatmap_legend_param = list(
        title_gp = gpar(fontsize = 12),
        labels_gp = gpar(fontsize = 10)
      )
    )
  })

  # Draw the heatmap
  ht_drawn <- draw(ht,
                   heatmap_legend_side = "bottom",
                   annotation_legend_side = "bottom",
                   legend_grouping = "original")

  # Print summary statistics
  n_features <- nrow(data_matrix)
  message("\nHeatmap Summary:")
  message("Total features plotted: ", n_features)

  for (test_name in names(available_tests)) {
    if (test_name %in% names(annotation_components)) {
      n_sig <- sum(annotation_components[[test_name]] == "Significant", na.rm = TRUE)
      message(test_name, " significant (", p_type, " < ", p_threshold, "): ", n_sig)
    }
  }

  return(ht_drawn)
}



#' heatmapFeatures()
#' @description Creates a heatmap visualization of selected features from a romics_object with hierarchical clustering and optional collapsing by factor groups.
#' @param romics_object A romics_object created with createRomicsObject()
#' @param factor Character string indicating which factor to use for organizing columns (default: "main" uses the main factor)
#' @param scale_feature Boolean indicating whether to scale feature values across samples (default: TRUE)
#' @param feature_list Character vector of feature names to include in the heatmap. If NULL, all features are used.
#' @param filter_by_stat_column Character string specifying a statistics column name to filter features (e.g., "mean_Group1")
#' @param stat_column_filter Character string specifying the filter condition as an expression (e.g., "> 100" or "< 0.05")
#' @param gradient_colors Character vector of colors for custom gradient, or NULL to use viridis_option
#' @param viridis_option Character string specifying viridis color palette: "viridis", "plasma", "inferno", "magma", "cividis", "rocket", "mako", or "turbo" (default: "viridis")
#' @param show_completeness Boolean indicating whether to display data completeness percentages as a heatmap annotation (default: TRUE)
#' @param show_dendrogram Boolean indicating whether to display hierarchical clustering dendrograms (default: TRUE)
#' @param show_feature_names Boolean indicating whether to display feature names on the heatmap (default: FALSE)
#' @param show_clusters Boolean indicating whether to color-code hierarchical clusters on the heatmap (default: FALSE)
#' @param show_cluster_legend Boolean indicating whether to display a legend for cluster colors (default: FALSE, only used when show_clusters=TRUE)
#' @param n_clusters Numeric value specifying the number of clusters to cut the dendrogram into (required when show_clusters=TRUE, must be >= 2)
#' @param clustering_method Character string specifying the hierarchical clustering method passed to hclust() (default: "ward.D")
#' @param collapse_by_factor Character string or TRUE indicating whether to collapse/average features by a different factor before creating the heatmap (default: NULL). If TRUE, uses the same factor as specified by 'factor' parameter.
#' @details This function creates a ComplexHeatmap visualization of features grouped by the specified factor. Features can be filtered by providing a statistics column and filter expression. The heatmap can optionally show data completeness as a side heatmap, hierarchical clustering dendrograms, and cluster assignments. Using collapse_by_factor creates a secondary level of averaging, useful for viewing group-level patterns.
#' @return A ComplexHeatmap heatmap object
#' @author Geremy Clair
#' @export
heatmapFeatures <- function(romics_object,
                            factor = "main",
                            scale_feature = TRUE,
                            feature_list = NULL,
                            filter_by_stat_column = NULL,
                            stat_column_filter = NULL,
                            gradient_colors = NULL,
                            viridis_option = "viridis",
                            show_completeness = TRUE,
                            show_dendrogram = TRUE,
                            show_feature_names = FALSE,
                            show_clusters = FALSE,
                            show_cluster_legend = FALSE,
                            n_clusters = NULL,
                            clustering_method = "ward.D",
                            collapse_by_factor = NULL) {
  # ---- VALIDATION ----
  if (!RomicsProcessor::is.romicsObject(romics_object)) {
    stop("Invalid romics_object")
  }
  if (missing(factor)) {
    factor <- "main"
  }
  if (!factor %in% c("main", RomicsProcessor::romicsFactorNames(romics_object))) {
    stop("The selected factor is not in the list of factors")
  }
  if (!is.logical(scale_feature)) {
    stop("scale_feature must be TRUE or FALSE")
  }
  if (!is.logical(show_completeness)) {
    stop("show_completeness must be TRUE or FALSE")
  }
  if (!is.logical(show_feature_names)) {
    stop("show_feature_names must be TRUE or FALSE")
  }
  if (!is.logical(show_dendrogram)) {
    stop("show_dendrogram must be TRUE or FALSE")
  }
  if (!is.logical(show_clusters)) {
    stop("show_clusters must be TRUE or FALSE")
  }
  if (!is.logical(show_cluster_legend)) {
    stop("show_cluster_legend must be TRUE or FALSE")
  }

  if (show_clusters && is.null(n_clusters)) {
    stop("n_clusters must be specified when show_clusters is TRUE")
  }

  if (!is.null(n_clusters) && (!is.numeric(n_clusters) || n_clusters < 2)) {
    stop("n_clusters must be a numeric value >= 2")
  }

  # Validate viridis_option
  valid_viridis_options <- c("viridis", "plasma", "inferno", "magma", "cividis", "rocket", "mako", "turbo")
  if (!viridis_option %in% valid_viridis_options) {
    stop("viridis_option must be one of: ", paste(valid_viridis_options, collapse = ", "))
  }

  if (is.null(feature_list)) {
    feature_list <- rownames(romics_object$data)
  }
  # ---- PREPARE OBJECT ----
  if (factor != "main") {
    romics_object <- RomicsProcessor::romicsChangeFactor(romics_object, main_factor = factor)
  }
  # Remove statistics to recalculate
  romics_object <- romics_object[!names(romics_object) == "statistics"]
  class(romics_object) <- "romics_object"
  # Calculate statistics
  romics_object <- RomicsProcessor::romicsMean(romics_object)
  romics_object <- RomicsProcessor::romicsPercentComplete(romics_object)
  df <- romics_object$statistics
  # ---- APPLY STAT-BASED FILTERING ----
  if (!is.null(filter_by_stat_column) && !is.null(stat_column_filter)) {
    # Find the stat column
    if (!filter_by_stat_column %in% colnames(df)) {
      stop("Stat column '", filter_by_stat_column, "' not found in statistics. ",
           "Available columns: ", paste(colnames(df)[1:min(5, ncol(df))], collapse = ", "), "...")
    }
    # Parse the filter expression
    stat_values <- df[[filter_by_stat_column]]
    # Evaluate the filter expression
    filter_result <- tryCatch({
      # Create a temporary variable for evaluation
      x <- stat_values
      eval(parse(text = paste("x", stat_column_filter)))
    }, error = function(e) {
      stop("Error parsing filter expression '", stat_column_filter, "': ", e$message)
    })
    # Apply filter
    if (!is.logical(filter_result)) {
      stop("Filter expression must return logical values")
    }
    cat("Filtering: ", sum(filter_result, na.rm = TRUE), " features passed filter out of ",
        length(filter_result), "\n")
    df <- df[filter_result, , drop = FALSE]
    if (nrow(df) == 0) {
      stop("No features passed the filtering criteria")
    }
  }
  # Filter by feature list
  feature_list <- feature_list[feature_list %in% rownames(df)]
  df <- df[rownames(df) %in% feature_list, , drop = FALSE]
  if (nrow(df) == 0) {
    stop("No features remaining after filtering")
  }
  df$Feature <- rownames(df)
  rownames(df) <- NULL

  # ---- HANDLE COLLAPSE BY FACTOR ----
  if (!is.null(collapse_by_factor) && collapse_by_factor != FALSE) {
    # If collapse_by_factor is provided, we need to further aggregate the group means
    # This creates a second level of averaging
    collapse_factor <- if(isTRUE(collapse_by_factor)) factor else collapse_by_factor
    collapsed_temp <- collapse_data_by_factor(romics_object,
                                             romics_object$data[rownames(romics_object$data) %in% df$Feature, ],
                                             collapse_factor)
    # Calculate means from collapsed data
    collapsed_means <- as.data.frame(collapsed_temp)
    collapsed_means$Feature <- rownames(collapsed_temp)
    rownames(collapsed_means) <- NULL

    # Scale if requested
    if (scale_feature) {
      for (col in colnames(collapsed_means)) {
        if (col != "Feature") {
          collapsed_means[[col]] <- scale(as.numeric(collapsed_means[[col]]))
        }
      }
    }

    # Calculate completeness for collapsed data
    collapsed_completeness <- as.data.frame(collapsed_temp)
    for (col in colnames(collapsed_completeness)) {
      original_samples <- which(romicsExtractFactor(romics_object, factor = collapse_by_factor) == col)
      completeness_vals <- rowMeans(!is.na(romics_object$data[rownames(romics_object$data) %in% df$Feature, original_samples, drop = FALSE])) * 100
      collapsed_completeness[[col]] <- completeness_vals
    }
    collapsed_completeness$Feature <- rownames(collapsed_temp)
    rownames(collapsed_completeness) <- NULL

    # Build heatmap data from collapsed data
    heatmap_data <- data.frame()
    groups <- colnames(collapsed_means)[colnames(collapsed_means) != "Feature"]
    for (i in seq_along(groups)) {
      temp_df <- data.frame(
        Feature = collapsed_means$Feature,
        Group = groups[i],
        Mean = collapsed_means[[groups[i]]],
        Percentage_completeness = collapsed_completeness[[groups[i]]],
        stringsAsFactors = FALSE
      )
      heatmap_data <- rbind(heatmap_data, temp_df)
    }
    mat_for_clustering <- as.matrix(collapsed_means[, colnames(collapsed_means) != "Feature"])
    rownames(mat_for_clustering) <- collapsed_means$Feature

  } else {
    # Original behavior: use factor-level means
    # ---- EXTRACT MEAN VALUES ----
    df_mean <- df[grepl("_mean|Feature", colnames(df))]
  if (scale_feature) {
    mean_cols <- grepl("_mean", colnames(df_mean))
    df_mean[, mean_cols] <- t(scale(t(df_mean[, mean_cols])))
  }
  colnames(df_mean) <- gsub("_mean", "", colnames(df_mean))
  # ---- EXTRACT COMPLETENESS VALUES ----
  df_completeness <- df[grepl("_percentage_completeness|Feature", colnames(df))]
  colnames(df_completeness) <- gsub("_percentage_completeness", "", colnames(df_completeness))
  # ---- BUILD HEATMAP DATA ----
  heatmap_data <- data.frame()
  groups <- colnames(df_mean)[colnames(df_mean) != "Feature"]
  for (i in seq_along(groups)) {
    temp_df <- data.frame(
      Feature = df$Feature,
      Group = groups[i],
      Mean = df_mean[[groups[i]]],
      Percentage_completeness = df_completeness[[groups[i]]],
      stringsAsFactors = FALSE
    )
    heatmap_data <- rbind(heatmap_data, temp_df)
  }
  heatmap_data <- as.data.frame(heatmap_data)
  # Mark completely missing data for later coloring
  heatmap_data$is_missing <- heatmap_data$Percentage_completeness == 0
  }  # Close the collapse_by_factor conditional

  # ---- CLUSTERING ON IMPUTED DATA ----
  # Create imputed matrix for clustering only
  if (is.null(collapse_by_factor)) {
    mat_for_clustering_original <- as.matrix(df_mean[, colnames(df_mean) != "Feature"])
    rownames(mat_for_clustering_original) <- df$Feature
    mat_for_clustering <- mat_for_clustering_original
  }
  # mat_for_clustering is already set in the collapse_by_factor block if used
  # Impute using 2 stdev below minimum
  mat_imputed <- mat_for_clustering
  for (col in 1:ncol(mat_imputed)) {
    col_data <- mat_imputed[, col]
    valid_data <- col_data[!is.na(col_data)]
    if (length(valid_data) > 0) {
      min_val <- min(valid_data, na.rm = TRUE)
      sd_val <- sd(valid_data, na.rm = TRUE)
      if (!is.na(sd_val) && sd_val > 0) {
        impute_val <- min_val - (2 * sd_val)
      } else {
        impute_val <- min_val * 0.5
      }
      # Replace NA and 0 values with imputation value
      mat_imputed[is.na(mat_imputed[, col]), col] <- impute_val
    }
  }
  # Perform clustering on imputed data
  clust <- tryCatch({
    stats::hclust(stats::dist(mat_imputed), method = clustering_method)
  }, error = function(e) {
    cat("Warning: Clustering failed on imputed data:", e$message, "\n")
    cat("Attempting clustering on original data...\n")
    tryCatch({
      stats::hclust(stats::dist(mat_for_clustering), method = clustering_method)
    }, error = function(e2) {
      cat("Warning: Clustering failed on original data too\n")
      NULL
    })
  })
  if (!is.null(clust)) {
    feature_order <- clust$labels[clust$order]
    dendrogram_obj <- clust
  } else {
    feature_order <- unique(heatmap_data$Feature)
    dendrogram_obj <- NULL
  }
  heatmap_data$Feature <- factor(heatmap_data$Feature, levels = feature_order)
  heatmap_data <- heatmap_data[order(heatmap_data$Feature), ]

  # ---- CREATE CLUSTER ASSIGNMENTS ----
  cluster_assignments <- NULL
  p_clusters <- NULL
  cluster_colors <- NULL
  clustering_results <- NULL

  if (show_clusters && !is.null(dendrogram_obj) && !is.null(n_clusters)) {
    # Cut dendrogram into clusters
    clusters <- stats::cutree(dendrogram_obj, k = n_clusters)

    # Create a data frame with feature-cluster mapping, ordered by dendrogram
    cluster_df <- data.frame(
      Feature = names(clusters),
      Cluster = clusters
    )
    cluster_df <- cluster_df[match(feature_order, cluster_df$Feature), ]

    # Create clustering results object
    cluster_list <- split(cluster_df$Feature, cluster_df$Cluster)
    clustering_results <- list(
      clusters = cluster_list,
      cluster_assignments = cluster_df,
      dendrogram = dendrogram_obj,
      n_clusters = n_clusters,
      method = clustering_method
    )

    cat("\nCluster assignments:\n")
    cat("Clusters:\n")
    for (i in seq_along(cluster_list)) {
      cat("  Cluster", i, ":", length(cluster_list[[i]]), "features\n")
    }

    # Create cluster sidebar
    cluster_colors <- grDevices::hcl.colors(n_clusters, palette = "Dark3")

    cluster_plot_df <- data.frame(
      Feature = factor(feature_order, levels = feature_order),
      Cluster = cluster_df$Cluster,
      y = seq_along(feature_order)
    )

    p_clusters <- ggplot2::ggplot(cluster_plot_df, ggplot2::aes(x = 1, y = y, fill = factor(Cluster))) +
      ggplot2::geom_tile(width = 1, height = 1) +
      ggplot2::scale_fill_manual(
        values = cluster_colors,
        name = "Cluster",
        guide = ggplot2::guide_legend(ncol = 1)
      ) +
      ggplot2::scale_y_continuous(
        breaks = seq_along(feature_order),
        limits = c(0.5, length(feature_order) + 0.5),
        expand = c(0, 0)
      ) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::theme_void() +
      ggplot2::theme(
        legend.position = "none",
        plot.margin = ggplot2::margin(0, 5, 5, 0)
      )
  }

  # Store title info
  plot_title <- paste(length(feature_order), "Features -", length(groups), "Groups")

  # ---- CREATE HEATMAP PLOT (WITHOUT TITLE AND Y-AXIS LABELS) ----
  if (show_completeness) {
    # With alpha transparency for completeness
    p_heatmap <- ggplot2::ggplot(heatmap_data, ggplot2::aes(
      x = Group,
      y = Feature,
      fill = Mean,
      alpha = Percentage_completeness,
      text = paste0(
        "<b>Feature:</b> ", Feature, "<br>",
        "<b>Group:</b> ", Group, "<br>",
        "<b>Mean:</b> ", round(Mean, 3), "<br>",
        "<b>Completeness:</b> ", round(Percentage_completeness, 1), "%"
      )
    )) +
      ggplot2::geom_tile(color = NA, size = 0) +
      ggplot2::scale_alpha_continuous(
        name = "Completeness (%)",
        range = c(0.3, 1),
        limits = c(0, 100)
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
        axis.text.y = ggplot2::element_blank(),
        axis.title = ggplot2::element_text(size = 12, face = "bold"),
        panel.grid = ggplot2::element_blank(),
        legend.position = "none",
        plot.margin = ggplot2::margin(0, 5, 5, 5)
      ) +
      ggplot2::labs(
        x = "Group",
        y = ""
      )
  } else {
    # Without alpha - use black for missing data
    p_heatmap <- ggplot2::ggplot(heatmap_data, ggplot2::aes(
      x = Group,
      y = Feature,
      fill = Mean,
      color = is_missing,
      text = paste0(
        "<b>Feature:</b> ", Feature, "<br>",
        "<b>Group:</b> ", Group, "<br>",
        "<b>Mean:</b> ", round(Mean, 3), "<br>",
        "<b>Completeness:</b> ", round(Percentage_completeness, 1), "%"
      )
    )) +
      ggplot2::geom_tile(size = 0) +
      ggplot2::scale_color_manual(
        values = c("FALSE" = NA, "TRUE" = "black"),
        guide = "none"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
        axis.text.y = ggplot2::element_blank(),
        axis.title = ggplot2::element_text(size = 12, face = "bold"),
        panel.grid = ggplot2::element_blank(),
        legend.position = "none",
        plot.margin = ggplot2::margin(0, 5, 5, 5)
      ) +
      ggplot2::labs(
        x = "Group",
        y = ""
      )
  }
  # ---- ADD COLOR SCALE ----
  if (is.null(gradient_colors)) {
    # Use viridis color scale with specified option
    if (scale_feature) {
      p_heatmap <- p_heatmap + ggplot2::scale_fill_viridis_c(
        name = "Scaled Mean",
        option = viridis_option,
        na.value = "grey90"
      )
    } else {
      p_heatmap <- p_heatmap + ggplot2::scale_fill_viridis_c(
        name = "Mean Abundance",
        option = viridis_option,
        na.value = "grey90"
      )
    }
  } else {
    # Use custom gradient colors
    if (length(gradient_colors) >= 3) {
      p_heatmap <- p_heatmap + ggplot2::scale_fill_gradient2(
        low = gradient_colors[1],
        mid = gradient_colors[2],
        high = gradient_colors[3],
        name = ifelse(scale_feature, "Scaled Mean", "Mean Abundance"),
        na.value = "grey90"
      )
    } else if (length(gradient_colors) >= 2) {
      p_heatmap <- p_heatmap + ggplot2::scale_fill_gradient(
        low = gradient_colors[1],
        high = gradient_colors[2],
        name = ifelse(scale_feature, "Scaled Mean", "Mean Abundance"),
        na.value = "grey90"
      )
    }
  }

  # ---- CREATE FEATURE NAMES PLOT ----
  p_feature_names <- NULL
  if (show_feature_names) {
    # Create a simple text plot with feature names on the right
    feature_names_df <- data.frame(
      Feature = factor(feature_order, levels = feature_order),
      x = 1,
      y = seq_along(feature_order)
    )

    p_feature_names <- ggplot2::ggplot(feature_names_df, ggplot2::aes(x = x, y = y, label = Feature)) +
      ggplot2::geom_text(hjust = 0, vjust = 0.5, size = 3.5) +
      ggplot2::scale_x_continuous(limits = c(0.8, 3), expand = c(0, 0)) +
      ggplot2::scale_y_continuous(
        breaks = seq_along(feature_order),
        labels = feature_order,
        limits = c(0.5, length(feature_order) + 0.5),
        expand = c(0, 0)
      ) +
      ggplot2::theme_void() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(0, 5, 5, 0)
      )
  }

  # ---- CREATE LEGEND PLOT ----
  # Extract legend from heatmap
  p_heatmap_with_legend <- p_heatmap +
    ggplot2::theme(
      legend.position = "left",
      legend.box = "vertical"
    )

  legend <- cowplot::get_legend(p_heatmap_with_legend)

  # Remove legend from heatmap
  p_heatmap <- p_heatmap +
    ggplot2::theme(legend.position = "none")

  # ---- CREATE CLUSTER LEGEND PLOT ----
  p_cluster_legend <- NULL
  if (show_clusters && show_cluster_legend && !is.null(p_clusters)) {
    cluster_legend_df <- data.frame(
      Cluster = 1:n_clusters,
      y = seq_along(1:n_clusters)
    )

    p_cluster_legend <- ggplot2::ggplot(cluster_legend_df, ggplot2::aes(x = 1, y = y, fill = factor(Cluster))) +
      ggplot2::geom_tile(width = 1, height = 0.8) +
      ggplot2::geom_text(ggplot2::aes(label = Cluster), vjust = 0.5, hjust = 0.5, size = 3, fontface = "bold", color = "white") +
      ggplot2::scale_fill_manual(
        values = cluster_colors,
        guide = "none"
      ) +
      ggplot2::scale_y_continuous(
        breaks = seq_along(1:n_clusters),
        limits = c(0.5, n_clusters + 0.5),
        expand = c(0, 0)
      ) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.margin = ggplot2::margin(0, 0, 5, 5)
      )
  }

  # ---- CREATE DENDROGRAM PLOT ----
  p_final <- p_heatmap

  if (!is.null(dendrogram_obj) && show_dendrogram) {
    tryCatch({
      if (requireNamespace("ggdendro", quietly = TRUE)) {
        # Get dendrogram data
        dendro_data <- ggdendro::dendro_data(dendrogram_obj, type = "rectangle")

        # Get number of features for scaling
        n_features <- length(feature_order)

        # Create dendrogram plot with proper scaling and orientation
        p_dendro <- ggplot2::ggplot(ggdendro::segment(dendro_data)) +
          ggplot2::geom_segment(ggplot2::aes(x = y, y = x, xend = yend, yend = xend),
                                size = 0.5, color = "gray40") +
          ggplot2::scale_y_continuous(
            breaks = seq_along(feature_order),
            labels = feature_order,
            limits = c(0.5, n_features + 0.5),
            expand = c(0, 0)
          ) +
          ggplot2::scale_x_reverse(expand = c(0, 0)) +
          ggplot2::theme_void() +
          ggplot2::theme(
            axis.text.y = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(0, 5, 5, 0),
            panel.spacing = ggplot2::unit(0, "lines")
          )

        # Build layout: legend | dendrogram | clusters | heatmap | rownames | cluster_legend
        plot_elements <- list(legend, p_dendro)
        rel_widths <- c(0.15, 0.2)

        if (show_clusters && !is.null(p_clusters)) {
          plot_elements <- c(plot_elements, list(p_clusters))
          rel_widths <- c(rel_widths, 0.06)
        }

        plot_elements <- c(plot_elements, list(p_heatmap))
        rel_widths <- c(rel_widths, 1)

        if (show_feature_names) {
          plot_elements <- c(plot_elements, list(p_feature_names))
          rel_widths <- c(rel_widths, 0.25)
        }

        if (show_clusters && show_cluster_legend && !is.null(p_cluster_legend)) {
          plot_elements <- c(plot_elements, list(p_cluster_legend))
          rel_widths <- c(rel_widths, 0.08)
        }

        p_combined <- do.call(cowplot::plot_grid, c(
          plot_elements,
          list(nrow = 1, rel_widths = rel_widths, align = "h", axis = "b")
        ))

        # Add title on top
        p_final <- cowplot::plot_grid(
          cowplot::ggdraw() +
            cowplot::draw_text(
              plot_title,
              x = 0.5,
              y = 0.5,
              hjust = 0.5,
              vjust = 0.5,
              size = 14,
              fontface = "bold"
            ),
          p_combined,
          nrow = 2,
          rel_heights = c(0.08, 1)
        )
      } else {
        # Dendrogram package not available
        plot_elements <- list(legend)
        rel_widths <- c(0.15)

        if (show_clusters && !is.null(p_clusters)) {
          plot_elements <- c(plot_elements, list(p_clusters))
          rel_widths <- c(rel_widths, 0.06)
        }

        plot_elements <- c(plot_elements, list(p_heatmap))
        rel_widths <- c(rel_widths, 1)

        if (show_feature_names) {
          plot_elements <- c(plot_elements, list(p_feature_names))
          rel_widths <- c(rel_widths, 0.25)
        }

        if (show_clusters && show_cluster_legend && !is.null(p_cluster_legend)) {
          plot_elements <- c(plot_elements, list(p_cluster_legend))
          rel_widths <- c(rel_widths, 0.08)
        }

        p_combined <- do.call(cowplot::plot_grid, c(
          plot_elements,
          list(nrow = 1, rel_widths = rel_widths, align = "h", axis = "b")
        ))

        p_final <- cowplot::plot_grid(
          cowplot::ggdraw() +
            cowplot::draw_text(
              plot_title,
              x = 0.5,
              y = 0.5,
              hjust = 0.5,
              vjust = 0.5,
              size = 14,
              fontface = "bold"
            ),
          p_combined,
          nrow = 2,
          rel_heights = c(0.08, 1)
        )
      }
    }, error = function(e) {
      cat("Warning: Failed to create dendrogram:", e$message, "\n")
      plot_elements <- list(legend)
      rel_widths <- c(0.15)

      if (show_clusters && !is.null(p_clusters)) {
        plot_elements <- c(plot_elements, list(p_clusters))
        rel_widths <- c(rel_widths, 0.06)
      }

      plot_elements <- c(plot_elements, list(p_heatmap))
      rel_widths <- c(rel_widths, 1)

      if (show_feature_names) {
        plot_elements <- c(plot_elements, list(p_feature_names))
        rel_widths <- c(rel_widths, 0.25)
      }

      if (show_clusters && show_cluster_legend && !is.null(p_cluster_legend)) {
        plot_elements <- c(plot_elements, list(p_cluster_legend))
        rel_widths <- c(rel_widths, 0.08)
      }

      p_combined <- do.call(cowplot::plot_grid, c(
        plot_elements,
        list(nrow = 1, rel_widths = rel_widths, align = "h", axis = "b")
      ))

      cowplot::plot_grid(
        cowplot::ggdraw() +
          cowplot::draw_text(
            plot_title,
            x = 0.5,
            y = 0.5,
            hjust = 0.5,
            vjust = 0.5,
            size = 14,
            fontface = "bold"
          ),
        p_combined,
        nrow = 2,
        rel_heights = c(0.08, 1)
      )
    })
  } else {
    # No dendrogram
    plot_elements <- list(legend)
    rel_widths <- c(0.15)

    if (show_clusters && !is.null(p_clusters)) {
      plot_elements <- c(plot_elements, list(p_clusters))
      rel_widths <- c(rel_widths, 0.06)
    }

    plot_elements <- c(plot_elements, list(p_heatmap))
    rel_widths <- c(rel_widths, 1)

    if (show_feature_names) {
      plot_elements <- c(plot_elements, list(p_feature_names))
      rel_widths <- c(rel_widths, 0.25)
    }

    if (show_clusters && show_cluster_legend && !is.null(p_cluster_legend)) {
      plot_elements <- c(plot_elements, list(p_cluster_legend))
      rel_widths <- c(rel_widths, 0.08)
    }

    p_combined <- do.call(cowplot::plot_grid, c(
      plot_elements,
      list(nrow = 1, rel_widths = rel_widths, align = "h", axis = "b")
    ))

    p_final <- cowplot::plot_grid(
      cowplot::ggdraw() +
        cowplot::draw_text(
          plot_title,
          x = 0.5,
          y = 0.5,
          hjust = 0.5,
          vjust = 0.5,
          size = 14,
          fontface = "bold"
        ),
      p_combined,
      nrow = 2,
      rel_heights = c(0.08, 1)
    )
  }

  # ---- RETURN RESULTS ----
  # Attach clustering results if clustering was performed
  if (show_clusters && !is.null(clustering_results)) {
    attr(p_final, "clustering_results") <- clustering_results
  }
  return(p_final)
}

#' romicsFoldChangeBarplot()
#' @description Creates fold change barplots with significance annotations from a romics object
#' @param romics_object A romics object containing $data, $statistics, and $metadata layers
#' @param feature_list Character vector of features to include in barplot. If NULL, uses all features (limited to 100 max)
#' @param comparison Character vector of comparisons to include. Default "all" includes all available comparisons. Format: c("groupA_vs_groupB", "groupB_vs_groupC")
#' @param test_type Which test to use for fold changes. Options: "ttest", "wilcoxon", "auto" (default: "auto" - uses first available)
#' @param use_adjusted_p TRUE or FALSE, if TRUE uses adjusted p-values, if FALSE uses raw p-values (default: TRUE)
#' @param p_threshold Numeric threshold for p-values (default 0.05)
#' @param make_horizontal TRUE or FALSE, if TRUE creates horizontal barplot (features on y-axis), if FALSE creates vertical (features on x-axis) (default: FALSE)
#' @param order_by How to order features. Options: "provided" (default), "alphabetical", "pvalue", "fold_change"
#' @param significant_only TRUE or FALSE, if TRUE only shows features significant in at least one comparison (default: FALSE)
#' @param bar_colors Character vector of colors for different comparisons. If NULL, uses default palette
#' @param show_values TRUE or FALSE, if TRUE shows fold change values on bars (default: FALSE)
#' @param title Character string for plot title. If NULL, auto-generates title
#' @param fc_limits Numeric vector of length 2 specifying y-axis limits for fold changes c(min, max). If NULL, uses automatic limits
#' @param use_theme_ROP TRUE or FALSE, if TRUE applies theme_ROP() to the plot, if FALSE uses theme_classic() (default: TRUE)
#' @details Creates barplots showing fold changes with significance stars. Uses T-test or Wilcoxon test results from statistics layer.
#' @return ggplot2 object
#' @author Geremy Clair
#' @export
romicsFoldChangeBarplot <- function(romics_object,
                                    feature_list = NULL,
                                    comparison = "all",
                                    test_type = "auto",
                                    use_adjusted_p = FALSE,
                                    p_threshold = 0.05,
                                    make_horizontal = FALSE,
                                    order_by = "provided",
                                    significant_only = FALSE,
                                    bar_colors = NULL,
                                    show_values = FALSE,
                                    title = NULL,
                                    fc_limits = NULL,
                                    use_theme_ROP = TRUE) {



  # Handle feature list
  if (is.null(feature_list)) {
    # Use all features from the romics object
    all_features <- rownames(romics_object$data)
    n_features <- length(all_features)

    if (n_features > 100) {
      stop("Too many features to display (", n_features, " features). Please provide a feature_list with ≤100 features, or use significant_only=TRUE to filter features automatically.")
    }

    feature_list <- all_features
    message("No feature_list provided. Using all ", n_features, " features from romics object.")
  }

  # Create feature subset if feature_list is provided
  if (!is.character(feature_list)) {
    stop("feature_list must be a character vector")
  }

  message("Creating subset with ", length(feature_list), " features...")
  romics_object <- romicsFeatureSubset(romics_object, feature_list = feature_list)
  message("Subset created with ", nrow(romics_object$data), " features")

  # Check if we have any features left
  if (nrow(romics_object$data) == 0) {
    stop("No features found after subsetting. Check that feature_list matches rownames in romics_object$data")
  }

  # Extract statistics
  stats_df <- romics_object$statistics
  stat_cols <- colnames(stats_df)

  # Auto-detect available comparisons
  vs_columns <- stat_cols[grepl("_vs_.*_(p|padj)$", stat_cols)]
  if (length(vs_columns) == 0) {
    stop("No comparison columns found in statistics")
  }

  available_comparisons <- unique(gsub("_(p|padj|Ttest|Wilcoxon).*", "", vs_columns))
  available_comparisons <- available_comparisons[grepl("_vs_", available_comparisons)]

  message("Available comparisons: ", paste(available_comparisons, collapse = ", "))

  # Select comparisons to use
  if (length(comparison) == 1 && comparison == "all") {
    selected_comparisons <- available_comparisons
  } else {
    selected_comparisons <- comparison
    missing_comparisons <- selected_comparisons[!selected_comparisons %in% available_comparisons]
    if (length(missing_comparisons) > 0) {
      warning("Some requested comparisons not found: ", paste(missing_comparisons, collapse = ", "))
      selected_comparisons <- selected_comparisons[selected_comparisons %in% available_comparisons]
    }
  }

  if (length(selected_comparisons) == 0) {
    stop("No valid comparisons found")
  }

  message("Using comparisons: ", paste(selected_comparisons, collapse = ", "))

  # Auto-detect test type
  if (test_type == "auto") {
    # Check which tests are available
    ttest_available <- any(grepl("Ttest", stat_cols))
    wilcoxon_available <- any(grepl("Wilcoxon", stat_cols))

    if (ttest_available) {
      test_type <- "ttest"
      test_suffix <- "Ttest"
    } else if (wilcoxon_available) {
      test_type <- "wilcoxon"
      test_suffix <- "Wilcoxon"
    } else {
      stop("No T-test or Wilcoxon test results found in statistics")
    }
  } else if (test_type == "ttest") {
    test_suffix <- "Ttest"
  } else if (test_type == "wilcoxon") {
    test_suffix <- "Wilcoxon"
  } else {
    stop("test_type must be 'auto', 'ttest', or 'wilcoxon'")
  }

  message("Using test type: ", test_type)

  # Extract fold changes and p-values for each comparison
  plot_data_list <- list()
  p_suffix <- if (use_adjusted_p) "_padj" else "_p"
  p_type <- if (use_adjusted_p) "padj" else "p"

  for (comp in selected_comparisons) {
    # Get fold change column
    fc_patterns <- c(paste0("log(", gsub("_vs_", "/", comp), ")"),
                     paste0(comp, "_log2FC"),
                     paste0(comp, "_logFC"))

    fc_col <- NULL
    for (pattern in fc_patterns) {
      if (pattern %in% stat_cols) {
        fc_col <- pattern
        break
      }
    }

    if (is.null(fc_col)) {
      warning("No fold change column found for comparison: ", comp)
      next
    }

    # Get p-value column
    p_col <- paste0(comp, "_", test_suffix, p_suffix)

    if (!p_col %in% stat_cols) {
      warning("No p-value column found for comparison: ", comp, " with test: ", test_type)
      next
    }

    # Extract data
    comp_data <- data.frame(
      Feature = rownames(stats_df),
      Comparison = comp,
      FoldChange = stats_df[[fc_col]],
      PValue = stats_df[[p_col]],
      stringsAsFactors = FALSE
    )

    plot_data_list[[comp]] <- comp_data
  }

  if (length(plot_data_list) == 0) {
    stop("No valid data extracted for any comparison")
  }

  # Combine all data
  plot_data <- do.call(rbind, plot_data_list)

  # Filter for significant features only if requested
  if (significant_only) {
    significant_features <- plot_data %>%
      filter(PValue < p_threshold) %>%
      pull(Feature) %>%
      unique()

    if (length(significant_features) == 0) {
      stop("No significant features found with p < ", p_threshold)
    }

    # Additional check for too many significant features
    if (length(significant_features) > 100) {
      warning("Found ", length(significant_features), " significant features (>100). Consider using a more stringent p_threshold.")
    }

    plot_data <- plot_data %>%
      filter(Feature %in% significant_features)

    message("Filtered to ", length(significant_features), " significant features")
  }

  # Final check for number of features to plot
  n_features_to_plot <- length(unique(plot_data$Feature))
  if (n_features_to_plot > 100) {
    warning("Plotting ", n_features_to_plot, " features (>100). Plot may be crowded. Consider filtering features.")
  }

  # Order features
  if (order_by == "provided" && !is.null(feature_list)) {
    # Use provided order, filter for available features
    available_features <- intersect(feature_list, unique(plot_data$Feature))
    plot_data$Feature <- factor(plot_data$Feature, levels = available_features)
  } else if (order_by == "alphabetical") {
    plot_data$Feature <- factor(plot_data$Feature, levels = sort(unique(plot_data$Feature)))
  } else if (order_by == "pvalue") {
    # Order by best p-value across comparisons
    feature_order <- plot_data %>%
      group_by(Feature) %>%
      summarise(MinPValue = min(PValue, na.rm = TRUE), .groups = "drop") %>%
      arrange(MinPValue) %>%
      pull(Feature)
    plot_data$Feature <- factor(plot_data$Feature, levels = feature_order)
  } else if (order_by == "fold_change") {
    # Order by fold change of first comparison
    first_comp <- selected_comparisons[1]
    feature_order <- plot_data %>%
      filter(Comparison == first_comp) %>%
      arrange(desc(abs(FoldChange))) %>%
      pull(Feature)
    plot_data$Feature <- factor(plot_data$Feature, levels = feature_order)
  } else {
    # Default: use current order
    plot_data$Feature <- factor(plot_data$Feature, levels = unique(plot_data$Feature))
  }

  # Add significance stars - FIXED VERSION
  plot_data$Significance <- case_when(
    is.na(plot_data$PValue) ~ "",
    plot_data$PValue < 0.001 ~ "***",
    plot_data$PValue < 0.01 ~ "**",
    plot_data$PValue < p_threshold ~ "*",
    TRUE ~ ""
  )

  # Debug significance assignment
  sig_counts <- table(plot_data$Significance, useNA = "ifany")
  message("Significance star counts: ", paste(names(sig_counts), "=", sig_counts, collapse = ", "))

  # Set up colors
  if (is.null(bar_colors)) {
    n_comparisons <- length(selected_comparisons)
    if (n_comparisons <= 8) {
      bar_colors <- brewer.pal(max(n_comparisons, 3), "Set2")[1:n_comparisons]
    } else {
      bar_colors <- rainbow(n_comparisons)
    }
    names(bar_colors) <- selected_comparisons
  } else {
    if (length(bar_colors) < length(selected_comparisons)) {
      warning("Not enough colors provided, recycling colors")
      bar_colors <- rep(bar_colors, length.out = length(selected_comparisons))
    }
    names(bar_colors) <- selected_comparisons[1:length(bar_colors)]
  }

  # Create the base plot - REMOVED BLACK BORDERS
  p <- ggplot(plot_data, aes(x = Feature, y = FoldChange, fill = Comparison)) +
    geom_bar(stat = "identity", position = "dodge", color = NA) +  # REMOVED color = "black"
    scale_fill_manual(values = bar_colors) +
    geom_hline(yintercept = 0, color = "black", size = 0.5) +
    labs(
      x = "Features",
      y = "Log Fold Change",
      fill = "Comparison"
    )

  # Apply theme
  if (use_theme_ROP) {
    # Check if theme_ROP function exists
    if (exists("theme_ROP", mode = "function")) {
      p <- p + theme_ROP()
      message("Applied theme_ROP() to plot")
    } else {
      warning("theme_ROP() function not found. Using theme_classic() instead.")
      p <- p + theme_classic()
    }
  } else {
    p <- p + theme_classic()
  }

  # Add specific theme adjustments for text rotation (after main theme)
  p <- p + theme(
    axis.text.x = element_text(angle = if (make_horizontal) 0 else 45,
                               hjust = if (make_horizontal) 0.5 else 1),
    plot.title = element_text(hjust = 0.5)  # Center title
  )

  # Apply fold change limits if provided
  if (!is.null(fc_limits)) {
    if (length(fc_limits) != 2 || fc_limits[1] >= fc_limits[2]) {
      warning("fc_limits should be a numeric vector of length 2 with min < max. Using automatic limits.")
    } else {
      p <- p + ylim(fc_limits[1], fc_limits[2])
      message("Applied fold change limits: [", fc_limits[1], ", ", fc_limits[2], "]")
    }
  }

  # Add significance stars - IMPROVED VERSION
  if (any(plot_data$Significance != "")) {
    # Filter for significant data and add debug info
    star_data <- plot_data %>%
      filter(Significance != "") %>%
      mutate(
        # Calculate star position above/below bars
        star_y = case_when(
          FoldChange >= 0 ~ FoldChange + 0.15 * max(abs(plot_data$FoldChange), na.rm = TRUE),
          FoldChange < 0 ~ FoldChange - 0.15 * max(abs(plot_data$FoldChange), na.rm = TRUE),
          TRUE ~ 0
        )
      )

    # Debug star data
    message("Stars to plot: ", nrow(star_data))
    if (nrow(star_data) > 0) {
      star_summary <- star_data %>%
        group_by(Significance) %>%
        summarise(count = n(), .groups = "drop")
      message("Star breakdown: ", paste(star_summary$Significance, "=", star_summary$count, collapse = ", "))
    }

    # Adjust star position if fc_limits are set
    if (!is.null(fc_limits) && nrow(star_data) > 0) {
      star_data$star_y <- pmax(pmin(star_data$star_y, fc_limits[2] * 0.95), fc_limits[1] * 0.95)
    }

    # Add stars to plot
    p <- p + geom_text(data = star_data,
                       aes(x = Feature, y = star_y, label = Significance),
                       position = position_dodge(width = 0.9),
                       size = 4, inherit.aes = FALSE, color = "black")
  }

  # Add values on bars if requested
  if (show_values) {
    p <- p + geom_text(aes(label = round(FoldChange, 2)),
                       position = position_dodge(width = 0.9),
                       vjust = ifelse(plot_data$FoldChange >= 0, -0.5, 1.5),
                       size = 3, color = "black")
  }

  # Handle horizontal orientation
  if (make_horizontal) {
    p <- p + coord_flip() +
      labs(x = "Log Fold Change", y = "Features") +  # Swap axis labels
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
            axis.text.y = element_text(angle = 0, hjust = 1))
  }

  # Add title
  if (is.null(title)) {
    title <- paste("Fold Change Analysis -", test_type, "test")
    if (significant_only) {
      title <- paste(title, "(Significant features only)")
    }
  }
  p <- p + ggtitle(title)

  # Add legend for significance
  p <- p + labs(caption = paste("Significance levels: * p <", p_threshold,
                                ", ** p < 0.01, *** p < 0.001"))

  # Print summary
  n_features <- length(unique(plot_data$Feature))
  n_sig_features <- length(unique(plot_data$Feature[plot_data$Significance != ""]))

  message("\nBarplot Summary:")
  message("Total features plotted: ", n_features)
  message("Significant features (", p_type, " < ", p_threshold, "): ", n_sig_features)
  message("Using test: ", test_type)
  message("Using ", ifelse(use_adjusted_p, "adjusted", "raw"), " p-values")
  message("Orientation: ", ifelse(make_horizontal, "Horizontal (features on y-axis)", "Vertical (features on x-axis)"))
  message("Theme: ", ifelse(use_theme_ROP && exists("theme_ROP"), "theme_ROP()", "theme_classic()"))

  return(p)
}

#' romicsGLMheatmap()
#' @description Generate a heatmap displaying the presence/absence pattern of features across samples with GLM binomial test significance annotations
#' @param romics_object A romics_object created using romicsCreateObject() containing calculated GLM binomial test statistics.
#' @param filter_p logical (TRUE or FALSE) to indicate if features should be filtered by GLM test p-values.
#' @param p_type 'p' or 'padj' to indicate the type of p-value to use for filtering and significance determination.
#' @param p numeric value indicating the p-value threshold for significance (between 0 and 1).
#' @param mode 'vs' to show up/down regulation patterns or 'enrichment' to show only enrichment (positive directionality).
#' @param hclust_row logical (TRUE or FALSE) to indicate whether to perform hierarchical clustering on rows.
#' @param row_annotation logical (TRUE or FALSE) to indicate whether to display row annotations showing significance status for each comparison (red = significant, grey = not significant).
#' @param column_annotation logical (TRUE or FALSE) to indicate whether to display column annotation bar showing the factor level for each sample.
#' @param show_rownames logical (TRUE or FALSE) to indicate whether to display feature names on the left side of the heatmap.
#' @param factor character or NULL to specify which factor to use for column annotations. If NULL, uses the main factor from the romics_object.
#' @details Generate a ComplexHeatmap showing presence/absence of features with GLM binomial test results. Row annotations show significance status using a unified color scheme (red/grey), and column annotations show factor levels using ROP_colors palette.
#' @return A ComplexHeatmap object
#' @author Geremy Clair
#' @export
romicsGLMheatmap <- function(romics_object,
                             filter_p = TRUE,
                             p_type = c("p", "padj"),
                             p = 0.05,
                             mode = c("vs", "enrichment"),
                             hclust_row = TRUE,
                             row_annotation = TRUE,
                             column_annotation = TRUE,
                             show_rownames = FALSE,
                             factor = NULL) {

  if(!is.romicsObject(romics_object)) {
    stop("Your <romics_object> was not in the right format or was not created using the function 'romicsCreateObject()'.")
  }

  if(!filter_p %in% c(TRUE, FALSE)) {
    stop("<filter_p> has to be TRUE or FALSE, to indicate if you'd like to filter the analytes by GLM test p-values.")
  }

  if(missing(filter_p)) {
    filter_p <- TRUE
  }

  if(length(p) != 1 | p > 1 | p < 0) {
    stop("<p> has to be numerical, higher than 0, and lower than 1.")
  }

  if(missing(p)) {
    p <- 0.05
  }

  if(missing(p_type)) {
    p_type <- "p"
  }

  if(!p_type %in% c("p", "padj")) {
    stop("<p_type> has to be either 'p' or 'padj'.")
  }

  if(missing(hclust_row)) {
    hclust_row <- TRUE
  }

  if(missing(row_annotation)) {
    row_annotation <- TRUE
  }

  if(missing(column_annotation)) {
    column_annotation <- TRUE
  }

  if(missing(show_rownames)) {
    show_rownames <- FALSE
  }

  if(!show_rownames %in% c(TRUE, FALSE)) {
    stop("<show_rownames> has to be TRUE or FALSE.")
  }

  if(missing(mode)) {
    mode <- "vs"
  }

  # Handle factor parameter
  if(!is.null(factor)) {
    if(!factor %in% romicsFactorNames(romics_object)) {
      stop(paste0("Factor '", factor, "' not found in romics_object. Available factors: ",
                  paste(romicsFactorNames(romics_object), collapse = ", ")))
    }
    romics_object <- romicsChangeFactor(romics_object, main_factor = factor)
  }

  # Let's extract the values that are "present" (1) and absent (0)
  m <- as.data.frame(!romics_object$missingdata)
  for(i in 1:ncol(m)) {
    m[, i] <- as.numeric(m[, i])
  }
  m <- as.matrix(m)

  if(is.null(romics_object$statistics)) {
    stop("Your <romics_object> was not containing calculated statistics.")
  }

  # Fixed typo: glmBinomialTest instead of glmBionomialTest
  if(sum(grepl("_glmBinomialTest_p", colnames(romics_object$statistics))) == 0) {
    stop("Your <romics_object> was not containing GLM binomial test statistics.")
  }

  s <- romics_object$statistics
  s <- s[, grepl("directionality", colnames(s)) | grepl("glmBinomialTest", colnames(s))]

  if(filter_p == TRUE) {
    if(p_type == "p") {
      s <- s[, grepl("directionality", colnames(s)) | !grepl("glmBinomialTest_padj", colnames(s))]
    } else {
      s <- s[, grepl("directionality", colnames(s)) | grepl("glmBinomialTest_padj", colnames(s))]
    }
  }

  sp <- s[, grepl("glmBinomialTest_p", colnames(s))]
  sd <- s[, grepl("directionality", colnames(s))]

  if(mode == "enrichment") {
    sb <- (sp < p) + (sd == 1) == 2
  } else {
    sb <- sp < p
  }

  sb_filter <- data.frame(keep = rowSums(sb) > 0)
  m <- m[sb_filter[, 1], , drop = FALSE]
  sd <- sd[sb_filter[, 1], , drop = FALSE]
  sp <- sp[sb_filter[, 1], , drop = FALSE]

  # Create binary significance indicator (significant vs not significant)
  sig_indicator <- sd
  if(mode == "vs") {
    sig_indicator[sd == 1 | sd == -1] <- "significant"
  } else {
    sig_indicator[sd == 1] <- "significant"
    sig_indicator[sd == -1] <- "not significant"
  }
  sig_indicator[sp > p] <- "not significant"

  sidenote <- data.frame(sig_indicator)
  colnames(sidenote) <- gsub("_directionality", "", colnames(sidenote))
  sidenote <- sidenote[do.call(order, as.list(sidenote)), , drop = FALSE]
  m <- m[order(match(rownames(m), rownames(sidenote))), , drop = FALSE]

  # Extract factor levels for column annotation
  factor_levels <- suppressMessages(romicsExtractFactor(romics_object))

  # Match factor levels to columns in m
  col_annotation <- factor_levels[match(colnames(m), names(factor_levels))]
  col_annotation_df <- data.frame(Factor = as.character(col_annotation))
  rownames(col_annotation_df) <- colnames(m)

  # Create column annotation
  if(column_annotation == TRUE) {
    ha_column <- ComplexHeatmap::HeatmapAnnotation(
      df = col_annotation_df,
      col = list(Factor = setNames(
        ROP_colors[1:length(unique(col_annotation_df$Factor))],
        unique(col_annotation_df$Factor)
      ))
    )
  } else {
    ha_column <- NULL
  }

  # Determine row names side based on row_annotation
  # If row annotations are shown, put row names on the left
  # If no row annotations, put row names on the right (default)
  if(row_annotation == TRUE && show_rownames == TRUE) {
    row_names_side <- "left"
  } else {
    row_names_side <- "right"
  }

  # Create heatmap
  if(hclust_row == TRUE) {
    p <- ComplexHeatmap::Heatmap(m,
                                 cluster_rows = TRUE,
                                 cluster_columns = FALSE,
                                 row_names_gp = grid::gpar(fontsize = 3),
                                 show_row_names = show_rownames,
                                 row_names_side = row_names_side,
                                 top_annotation = ha_column)
  } else {
    p <- ComplexHeatmap::Heatmap(m,
                                 cluster_rows = FALSE,
                                 cluster_columns = FALSE,
                                 row_names_gp = grid::gpar(fontsize = 3),
                                 show_row_names = show_rownames,
                                 row_names_side = row_names_side,
                                 top_annotation = ha_column)
  }

  # Add row annotations with proper names and uniform colors
  if(row_annotation == TRUE) {
    # Define unified color scheme
    annot_colors <- c("significant" = "red", "not significant" = "grey")

    # Build the annotation arguments
    annot_args <- list()
    col_args <- list()

    for (i in 1:ncol(sidenote)) {
      comparison_name <- colnames(sidenote)[i]
      annot_args[[comparison_name]] <- sidenote[, i]
      col_args[[comparison_name]] <- annot_colors
      print(paste0("The sidenote #", i, " represents ", comparison_name, "."))
    }

    # Add the color list
    annot_args$col <- col_args
    annot_args$show_annotation_name <- TRUE
    annot_args$annotation_name_side <- "top"

    # Create the row annotation with all columns at once
    ha_row <- do.call(ComplexHeatmap::rowAnnotation, annot_args)

    p <- p + ha_row
  }

  return(p)
}
