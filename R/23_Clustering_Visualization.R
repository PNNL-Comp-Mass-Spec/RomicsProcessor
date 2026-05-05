#' romicsClusteringLinePlot()
#' @description Creates line plots to visualize feature behavior within clusters across samples
#' @param romics_object A romics_object created using createRomicsObject()
#' @param cluster_column Character indicating which clustering column to use from statistics layer
#' (e.g., "Cmeans_cluster", "hclust_clusters", or any other cluster assignment column)
#' @param factor_name Character indicating which metadata factor to use for x-axis ordering (default: "main")
#' @param color_palette Character vector of colors or viridis palette name (default: "viridis")
#' @param line_alpha Numeric (0-1) indicating base transparency of lines (default: 0.3)
#' @param line_size Numeric indicating line width (default: 0.7)
#' @param scale_data Logical indicating whether to scale features to 0-1 range within each cluster for visualization (default: TRUE)
#' @param facet_scales Character indicating whether facet scales should be "fixed" or "free_y" (default: "fixed")
#' @param show_mean Logical indicating whether to overlay cluster mean line (default: TRUE)
#' @param mean_line_size Numeric indicating thickness of mean line (default: 1.5)
#' @param aggregate_by_factor Character indicating metadata factor to aggregate samples by before plotting (default: NULL, no aggregation). When provided, samples are averaged within each factor level
#' @param aggregate_order Numeric vector specifying the order to arrange aggregated conditions (default: NULL). Only used with aggregate_by_factor. Example: c(1,3,4,2,5) reorders the first condition to position 1, third to position 2, etc.
#' @details This function creates line plots showing the abundance trajectory of all features within each cluster.
#' Multiple lines (one per feature) are overlaid with transparency to show line density.
#' The viridis color palette is used by default to distinguish clusters while maintaining perceptual uniformity.
#' Features are standardized within each cluster if scale_data=TRUE to allow comparison across clusters with different scales.
#' @return Returns a ggplot object showing feature trajectories colored by cluster with line density visualization
#' @author Geremy Clair
#' @export
#' @examples
#' # Basic clustering line plot using Fuzzy C-means results
#' plot <- romicsClusteringLinePlot(
#'   romics_object,
#'   cluster_column = "Cmeans_cluster"
#' )
#' plot
#'
#' # With hierarchical clustering results and custom styling
#' plot <- romicsClusteringLinePlot(
#'   romics_object,
#'   cluster_column = "hclust_clusters",
#'   factor_name = "treatment",
#'   line_alpha = 0.2,
#'   scale_data = TRUE,
#'   show_mean = TRUE
#' )
#' plot
romicsClusteringLinePlot <- function(romics_object,
                                    cluster_column,
                                    factor_name = "main",
                                    color_palette = "viridis",
                                    line_alpha = 0.3,
                                    line_size = 0.7,
                                    scale_data = TRUE,
                                    facet_scales = "fixed",
                                    show_mean = TRUE,
                                    mean_line_size = 1.5,
                                    aggregate_by_factor = NULL,
                                    aggregate_order = NULL) {

  # Input validation
  if(!is.romicsObject(romics_object) || missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }

  if(missing(cluster_column)) {
    stop("cluster_column must be specified")
  }

  if(!cluster_column %in% colnames(romics_object$statistics)) {
    stop(paste("Cluster column", cluster_column, "not found in statistics layer"))
  }

  if(!is.character(factor_name) || length(factor_name) != 1) {
    stop("factor_name must be a single character string")
  }

  if(factor_name == "main") {
    factor_name <- romics_object$main_factor
  }

  if(!factor_name %in% rownames(romics_object$metadata)) {
    stop(paste("Factor", factor_name, "not found in metadata"))
  }

  if(!is.numeric(line_alpha) || line_alpha < 0 || line_alpha > 1) {
    stop("line_alpha must be numeric between 0 and 1")
  }

  if(!is.numeric(line_size) || line_size <= 0) {
    stop("line_size must be a positive numeric value")
  }

  if(!is.numeric(mean_line_size) || mean_line_size <= 0) {
    stop("mean_line_size must be a positive numeric value")
  }

  if(!is.logical(scale_data)) {
    stop("scale_data must be TRUE or FALSE")
  }

  if(!is.logical(show_mean)) {
    stop("show_mean must be TRUE or FALSE")
  }

  if(!facet_scales %in% c("fixed", "free_y")) {
    stop("facet_scales must be either 'fixed' or 'free_y'")
  }

  if(!is.null(aggregate_by_factor) && !is.character(aggregate_by_factor)) {
    stop("aggregate_by_factor must be NULL or a character string matching a metadata field")
  }

  if(!is.null(aggregate_order) && (!is.numeric(aggregate_order) || !all(aggregate_order == as.integer(aggregate_order)))) {
    stop("aggregate_order must be a numeric vector of integers")
  }

  message("Creating clustering line plot...")
  message("  Cluster column: ", cluster_column)
  message("  Factor: ", factor_name)
  message("  Color palette: ", color_palette)

  # Extract data and clusters
  data <- romics_object$data
  clusters <- romics_object$statistics[[cluster_column]]

  # Handle factor aggregation
  if(!is.null(aggregate_by_factor)) {
    message("  Aggregating samples by factor: ", aggregate_by_factor)
    if(!aggregate_by_factor %in% rownames(romics_object$metadata)) {
      stop("Factor '", aggregate_by_factor, "' not found in metadata")
    }
    factor_vals <- romicsExtractFactor(romics_object, factor = aggregate_by_factor)
    unique_factors <- unique(factor_vals)
    agg_data <- matrix(0, nrow = nrow(data), ncol = length(unique_factors))
    colnames(agg_data) <- as.character(unique_factors)
    rownames(agg_data) <- rownames(data)
    for(i in seq_along(unique_factors)) {
      factor_level <- unique_factors[i]
      mask <- factor_vals == factor_level
      agg_data[, i] <- rowMeans(data[, mask, drop = FALSE], na.rm = TRUE)
    }

    # Apply custom ordering if provided
    if(!is.null(aggregate_order)) {
      if(length(aggregate_order) != length(unique_factors)) {
        stop("aggregate_order must have length equal to number of unique conditions (",
             length(unique_factors), ")")
      }
      if(!all(sort(aggregate_order) == seq_along(unique_factors))) {
        stop("aggregate_order must be a permutation of 1:", length(unique_factors))
      }
      agg_data <- agg_data[, aggregate_order, drop = FALSE]
      unique_factors <- unique_factors[aggregate_order]
      message("  Reordered conditions to: ", paste(aggregate_order, collapse = ", "))
    }

    data <- agg_data
    factor_values <- rep(as.character(unique_factors), times = 1)
  } else {
    if(!is.null(aggregate_order)) {
      stop("aggregate_order can only be used with aggregate_by_factor")
    }
    factor_values <- romicsExtractFactor(romics_object, factor = factor_name)
  }

  # Order samples by factor (skip if aggregated, since columns are already ordered)
  if(!is.null(aggregate_by_factor)) {
    sample_order <- seq_len(ncol(data))
  } else {
    sample_order <- order(as.character(factor_values))
  }
  data_ordered <- data[, sample_order]
  factor_ordered <- factor_values[sample_order]

  # Get unique clusters
  unique_clusters <- unique(clusters)
  n_na <- sum(is.na(unique_clusters))
  unique_clusters <- unique_clusters[!is.na(unique_clusters)]
  n_clusters <- length(unique_clusters)

  if(n_na > 0) {
    message("Found ", n_clusters, " clusters (", sum(is.na(clusters)), " features with NA clusters excluded)")
  } else {
    message("Found ", n_clusters, " clusters")
  }

  # Build plot data efficiently
  plot_rows <- list()
  row_idx <- 1

  # Process each cluster
  for(cluster in unique_clusters) {
    # Get features in this cluster
    cluster_mask <- clusters == cluster
    cluster_data <- data[cluster_mask, sample_order, drop = FALSE]

    if(nrow(cluster_data) == 0) next

    # Scale within cluster if requested
    if(scale_data) {
      for(feat in rownames(cluster_data)) {
        row_data <- as.numeric(cluster_data[feat, ])
        min_val <- min(row_data, na.rm = TRUE)
        max_val <- max(row_data, na.rm = TRUE)
        if(max_val > min_val) {
          cluster_data[feat, ] <- (row_data - min_val) / (max_val - min_val)
        } else {
          cluster_data[feat, ] <- 0.5
        }
      }
    }

    # Convert to long format (build rows list)
    for(i in seq_len(nrow(cluster_data))) {
      feat_name <- rownames(cluster_data)[i]
      values <- as.numeric(cluster_data[i, ])
      for(j in seq_along(values)) {
        if(!is.na(values[j])) {
          plot_rows[[row_idx]] <- list(
            feature = feat_name,
            sample_index = j,
            abundance = values[j],
            cluster = as.character(cluster),
            factor = as.character(factor_ordered[j])
          )
          row_idx <- row_idx + 1
        }
      }
    }
  }

  # Convert list to data frame (much faster than rbind loop)
  if(length(plot_rows) > 0) {
    plot_data <- data.frame(
      feature = sapply(plot_rows, `[[`, "feature"),
      sample_index = as.numeric(sapply(plot_rows, `[[`, "sample_index")),
      abundance = as.numeric(sapply(plot_rows, `[[`, "abundance")),
      cluster = sapply(plot_rows, `[[`, "cluster"),
      factor = sapply(plot_rows, `[[`, "factor"),
      stringsAsFactors = FALSE
    )
  } else {
    plot_data <- data.frame(feature = character(), sample_index = numeric(),
                            abundance = numeric(), cluster = character(),
                            factor = character())
  }

  # Create sample/condition mapping subtitle
  subtitle_text <- NULL
  if(!is.null(aggregate_by_factor)) {
    unique_conditions <- unique(factor_values)
    condition_mapping <- paste(seq_along(unique_conditions), ": ", unique_conditions, collapse = ", ", sep = "")
    subtitle_text <- paste("Condition order:", condition_mapping)
  } else {
    # Show sample order mapping for non-aggregated data
    sample_names <- colnames(romics_object$data)
    sample_order <- order(as.character(factor_values))
    ordered_samples <- sample_names[sample_order]
    # Truncate if too many samples
    if(length(ordered_samples) > 10) {
      sample_mapping <- paste(c(paste(1:10, ": ", ordered_samples[1:10], sep = ""),
                                paste("... (", length(ordered_samples) - 10, " more samples)", sep = "")),
                             collapse = ", ")
    } else {
      sample_mapping <- paste(seq_along(ordered_samples), ": ", ordered_samples, collapse = ", ", sep = "")
    }
    subtitle_text <- paste("Sample order:", sample_mapping)
  }

  # Create the plot with individual features in light gray background
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = sample_index, y = abundance, group = feature)) +
    ggplot2::geom_line(color = "lightgray", alpha = line_alpha * 0.5, linewidth = line_size * 0.5) +
    ggplot2::facet_wrap(~cluster, scales = facet_scales, nrow = ceiling(sqrt(n_clusters))) +
    ggplot2::labs(
      title = "Feature Clustering Line Plot",
      subtitle = subtitle_text,
      x = paste0("Sample Order (", if(!is.null(aggregate_by_factor)) aggregate_by_factor else factor_name, ")"),
      y = if(scale_data) "Scaled Abundance (0-1)" else "Abundance",
      color = "Cluster"
    ) +
    theme_ROP() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      strip.text = ggplot2::element_text(size = 11, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9, hjust = 0, margin = ggplot2::margin(b = 10))
    )

  # Add color palette
  if(length(color_palette) == 1 && color_palette %in% c("viridis", "plasma", "inferno", "magma", "cividis")) {
    p <- p + ggplot2::scale_color_viridis_d(option = color_palette, end = 0.9)
  } else if(length(color_palette) > 1) {
    # Use provided color palette
    p <- p + ggplot2::scale_color_manual(values = color_palette)
  } else {
    # Default to viridis
    p <- p + ggplot2::scale_color_viridis_d(end = 0.9)
  }

  # Add mean lines if requested
  if(show_mean) {
    mean_data <- data.frame()
    for(cluster in unique_clusters) {
      cluster_mask <- clusters == cluster
      cluster_data_subset <- data[cluster_mask, sample_order, drop = FALSE]

      if(scale_data) {
        for(feat in rownames(cluster_data_subset)) {
          row_data <- as.numeric(cluster_data_subset[feat, ])
          min_val <- min(row_data, na.rm = TRUE)
          max_val <- max(row_data, na.rm = TRUE)
          if(max_val > min_val) {
            cluster_data_subset[feat, ] <- (row_data - min_val) / (max_val - min_val)
          } else {
            cluster_data_subset[feat, ] <- 0.5
          }
        }
      }

      # Calculate mean
      means <- colMeans(cluster_data_subset, na.rm = TRUE)
      for(j in seq_along(means)) {
        mean_data <- rbind(mean_data, data.frame(
          sample_index = j,
          abundance = means[j],
          cluster = as.character(cluster),
          stringsAsFactors = FALSE
        ))
      }
    }

    p <- p + ggplot2::geom_line(
      data = mean_data,
      ggplot2::aes(x = sample_index, y = abundance, group = cluster, color = cluster),
      linewidth = mean_line_size * 3,
      alpha = 0.5,
      inherit.aes = FALSE
    )
  }

  message("Plot created successfully!")
  return(p)
}

#' romicsClusteringHeatmap()
#' @description Creates a heatmap visualization of feature clusters, optionally aggregated by a factor
#' @param romics_object A romics_object created using createRomicsObject()
#' @param cluster_column Character indicating which clustering column to use from statistics layer
#' (e.g., "Cmeans_cluster", "hclust_clusters", or any other cluster assignment column)
#' @param aggregate_by_factor Character indicating metadata factor to aggregate samples by before plotting (default: NULL, no aggregation)
#' @param aggregate_order Numeric vector specifying the order to arrange aggregated conditions (default: NULL). Only used with aggregate_by_factor
#' @param scale Logical indicating whether to scale data before plotting (default: TRUE)
#' @param scale_method Character indicating scaling method: "z-score" or "min-max" (default: "z-score", only used if scale=TRUE)
#' @param hclust_within Logical indicating whether to hierarchically cluster features within each cluster (default: TRUE)
#' @param hclust_columns Logical indicating whether to hierarchically cluster samples/conditions (columns) (default: FALSE)
#' @param hclust_dist_method Character indicating distance method for hierarchical clustering: "euclidean", "manhattan", "correlation" (default: "euclidean")
#' @param hclust_agglom_method Character indicating agglomeration method for hierarchical clustering: "ward.D2", "complete", "average", "single" (default: "ward.D2")
#' @param color_palette Character indicating color palette name (default: "viridis"). Supported palettes:
#' "viridis" (sequential), "RdBu" (diverging, good for scaled data), "RdYlBu" (diverging),
#' "RdYlGn" (diverging), "Spectral" (diverging). Other palette names default to viridis.
#' @details Creates a ComplexHeatmap showing feature abundances arranged by clusters. When aggregated, shows averaged values across samples within each factor level.
#' Features are annotated with their cluster ID. Optional within-cluster hierarchical clustering organizes features by similarity.
#' Optional column clustering groups similar samples/conditions together.
#' For scaled data, diverging palettes like "RdBu" help distinguish positive and negative deviations from the mean.
#' @return Returns a ComplexHeatmap object
#' @author Geremy Clair
#' @export
#' @examples
#' # Basic cluster heatmap with aggregation
#' hm <- romicsClusteringHeatmap(
#'   romics_object,
#'   cluster_column = "Cmeans_cluster",
#'   aggregate_by_factor = "condition"
#' )
#' hm
#'
#' # With scaling and within-cluster clustering
#' hm <- romicsClusteringHeatmap(
#'   romics_object,
#'   cluster_column = "Cmeans_cluster",
#'   aggregate_by_factor = "condition",
#'   scale = TRUE,
#'   hclust_within = TRUE
#' )
#' hm
#'
#' # With reordered conditions and diverging color palette for scaled data
#' hm <- romicsClusteringHeatmap(
#'   romics_object,
#'   cluster_column = "Cmeans_cluster",
#'   aggregate_by_factor = "condition",
#'   aggregate_order = c(1, 3, 2),
#'   color_palette = "RdBu"
#' )
#' hm
romicsClusteringHeatmap <- function(romics_object,
                                   cluster_column,
                                   aggregate_by_factor = NULL,
                                   aggregate_order = NULL,
                                   scale = TRUE,
                                   scale_method = "z-score",
                                   hclust_within = TRUE,
                                   hclust_columns = FALSE,
                                   hclust_dist_method = "euclidean",
                                   hclust_agglom_method = "ward.D2",
                                   color_palette = "viridis") {

  # Input validation
  if(!is.romicsObject(romics_object) || missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }

  if(missing(cluster_column)) {
    stop("cluster_column must be specified")
  }

  if(!cluster_column %in% colnames(romics_object$statistics)) {
    stop(paste("Cluster column", cluster_column, "not found in statistics layer"))
  }

  if(!is.logical(scale)) {
    stop("scale must be TRUE or FALSE")
  }

  if(!scale_method %in% c("z-score", "min-max")) {
    stop("scale_method must be 'z-score' or 'min-max'")
  }

  if(!is.logical(hclust_within)) {
    stop("hclust_within must be TRUE or FALSE")
  }

  if(!is.logical(hclust_columns)) {
    stop("hclust_columns must be TRUE or FALSE")
  }

  if(!hclust_dist_method %in% c("euclidean", "manhattan", "correlation")) {
    stop("hclust_dist_method must be one of: 'euclidean', 'manhattan', 'correlation'")
  }

  if(!hclust_agglom_method %in% c("ward.D2", "complete", "average", "single")) {
    stop("hclust_agglom_method must be one of: 'ward.D2', 'complete', 'average', 'single'")
  }

  if(!is.null(aggregate_by_factor) && !is.character(aggregate_by_factor)) {
    stop("aggregate_by_factor must be NULL or a character string matching a metadata field")
  }

  if(!is.null(aggregate_order) && (!is.numeric(aggregate_order) || !all(aggregate_order == as.integer(aggregate_order)))) {
    stop("aggregate_order must be a numeric vector of integers")
  }

  message("Creating cluster heatmap...")
  message("  Cluster column: ", cluster_column)
  message("  Scale: ", scale, if(scale) paste(" (", scale_method, ")", sep = "") else "")

  # Extract data and clusters
  data <- romics_object$data
  clusters <- romics_object$statistics[[cluster_column]]

  # Handle factor aggregation
  if(!is.null(aggregate_by_factor)) {
    message("  Aggregating samples by factor: ", aggregate_by_factor)
    if(!aggregate_by_factor %in% rownames(romics_object$metadata)) {
      stop("Factor '", aggregate_by_factor, "' not found in metadata")
    }
    factor_vals <- romicsExtractFactor(romics_object, factor = aggregate_by_factor)
    unique_factors <- unique(factor_vals)
    agg_data <- matrix(0, nrow = nrow(data), ncol = length(unique_factors))
    colnames(agg_data) <- as.character(unique_factors)
    rownames(agg_data) <- rownames(data)
    for(i in seq_along(unique_factors)) {
      factor_level <- unique_factors[i]
      mask <- factor_vals == factor_level
      agg_data[, i] <- rowMeans(data[, mask, drop = FALSE], na.rm = TRUE)
    }

    # Apply custom ordering if provided
    if(!is.null(aggregate_order)) {
      if(length(aggregate_order) != length(unique_factors)) {
        stop("aggregate_order must have length equal to number of unique conditions")
      }
      if(!all(sort(aggregate_order) == seq_along(unique_factors))) {
        stop("aggregate_order must be a permutation of 1:", length(unique_factors))
      }
      agg_data <- agg_data[, aggregate_order, drop = FALSE]
      message("  Reordered conditions to: ", paste(aggregate_order, collapse = ", "))
    }

    data <- agg_data
  }

  # Apply scaling if requested
  if(scale) {
    message("Scaling data using ", scale_method, "...")
    if(scale_method == "z-score") {
      data <- t(scale(t(data)))
    } else if(scale_method == "min-max") {
      data <- t(apply(data, 1, function(x) {
        min_val <- min(x, na.rm = TRUE)
        max_val <- max(x, na.rm = TRUE)
        if(max_val > min_val) {
          (x - min_val) / (max_val - min_val)
        } else {
          rep(0.5, length(x))
        }
      }))
    }
  }

  # Replace NAs with -20 for clustering (same as romicsHeatmap method)
  # This allows distance calculation while pushing missing values far from the distribution
  data_for_clustering <- data
  data_for_clustering[is.na(data_for_clustering)] <- -20

  # Prepare data for heatmap - order by cluster
  unique_clusters <- sort(unique(clusters[!is.na(clusters)]))
  n_clusters <- length(unique_clusters)

  # Build heatmap data with cluster annotations
  heatmap_data <- NULL
  cluster_annotations <- NULL
  row_order <- NULL

  for(clust in unique_clusters) {
    cluster_mask <- clusters == clust
    cluster_data <- data_for_clustering[cluster_mask, , drop = FALSE]

    # Hierarchical clustering within cluster if requested
    if(hclust_within && nrow(cluster_data) > 1) {
      hc <- stats::hclust(stats::dist(cluster_data))
      cluster_data <- data[cluster_mask, , drop = FALSE][hc$order, , drop = FALSE]
    }

    heatmap_data <- rbind(heatmap_data, cluster_data)
    cluster_annotations <- c(cluster_annotations, rep(as.character(clust), nrow(cluster_data)))
    row_order <- c(row_order, which(cluster_mask)[if(hclust_within && nrow(cluster_data) > 1) hc$order else seq_len(sum(cluster_mask))])
  }

  # Set default color palette
  if(is.null(color_palette)) {
    color_palette <- if(scale) "RdBu" else "viridis"
  }

  # Create color mapping for clusters (same as line plot)
  cluster_colors <- NULL
  if(color_palette == "viridis") {
    cluster_colors <- viridis::viridis(n_clusters, end = 0.9)
  } else if(color_palette %in% c("RdBu", "RdYlBu", "RdYlGn", "Spectral")) {
    # Use RColorBrewer palettes
    if(requireNamespace("RColorBrewer", quietly = TRUE)) {
      cluster_colors <- RColorBrewer::brewer.pal(min(n_clusters, 11), color_palette)
      if(n_clusters > length(cluster_colors)) {
        cluster_colors <- colorRampPalette(cluster_colors)(n_clusters)
      }
    } else {
      cluster_colors <- viridis::viridis(n_clusters, end = 0.9)
    }
  } else {
    cluster_colors <- viridis::viridis(n_clusters, end = 0.9)
  }

  # Create cluster color mapping
  cluster_color_map <- setNames(cluster_colors, as.character(unique_clusters))
  cluster_anno_colors <- cluster_color_map[cluster_annotations]

  # Create row annotation
  row_anno <- ComplexHeatmap::rowAnnotation(
    Cluster = cluster_annotations,
    col = list(Cluster = cluster_color_map),
    simple_anno_size = ggplot2::unit(5, "mm")
  )

  # Determine heatmap colors based on palette parameter
  if(color_palette == "viridis") {
    hm_colors <- viridis::viridis(256)
  } else if(color_palette %in% c("RdBu", "RdYlBu", "RdYlGn", "Spectral")) {
    # Use RColorBrewer palettes for diverging colors (useful when scaled)
    if(requireNamespace("RColorBrewer", quietly = TRUE)) {
      pal <- RColorBrewer::brewer.pal(11, color_palette)
      hm_colors <- circlize::colorRamp2(seq(-2, 2, length.out = 11), pal)
    } else {
      # Fallback to viridis if RColorBrewer not available
      hm_colors <- viridis::viridis(256)
    }
  } else {
    # Default to viridis for any other palette name
    hm_colors <- viridis::viridis(256)
  }

  # Calculate column dendrogram if requested
  col_dend <- NULL
  if(hclust_columns && ncol(heatmap_data) > 1) {
    message("Clustering samples/columns...")
    # Calculate distance based on selected method
    if(hclust_dist_method == "correlation") {
      col_dist <- stats::as.dist(1 - stats::cor(heatmap_data))
    } else {
      col_dist <- stats::dist(t(heatmap_data), method = hclust_dist_method)
    }
    col_hc <- stats::hclust(col_dist, method = hclust_agglom_method)
    col_dend <- as.dendrogram(col_hc)
  }

  # Create heatmap
  message("Creating ComplexHeatmap...")
  hm <- ComplexHeatmap::Heatmap(
    heatmap_data,
    name = if(scale) "Scaled\nAbundance" else "Abundance",
    col = hm_colors,
    left_annotation = row_anno,
    cluster_rows = FALSE,
    cluster_columns = if(hclust_columns && !is.null(col_dend)) col_dend else FALSE,
    show_row_names = FALSE,
    column_names_side = "bottom",
    column_names_rot = 45,
    column_names_centered = TRUE
  )

  message("Heatmap created successfully!")
  return(hm)
}

