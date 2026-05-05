# Helper function to aggregate data by factor
.aggregateRomicsByFactor <- function(romics_object, factor_name) {
  factor_vals <- romicsExtractFactor(romics_object, factor = factor_name)
  unique_factors <- unique(factor_vals)

  agg_data <- matrix(0, nrow = nrow(romics_object$data), ncol = length(unique_factors))
  colnames(agg_data) <- as.character(unique_factors)
  rownames(agg_data) <- rownames(romics_object$data)

  for(i in seq_along(unique_factors)) {
    factor_level <- unique_factors[i]
    mask <- factor_vals == factor_level
    agg_data[, i] <- rowMeans(romics_object$data[, mask, drop = FALSE], na.rm = TRUE)
  }

  list(data = agg_data, factor_map = setNames(as.character(unique_factors),
                                               as.character(colnames(agg_data))))
}

#' romicsFuzzyCMeans()
#' @description Performs Fuzzy C-means clustering on features based on their abundance behavior across samples
#' @param romics_object A romics_object created using createRomicsObject()
#' @param n_clusters Numeric indicating the number of clusters to create (default: 3)
#' @param m Numeric indicating the fuzziness parameter (default: 1.25). Higher values increase fuzziness. Must be > 1
#' @param iter_max Numeric indicating maximum number of iterations (default: 100)
#' @param standardize Logical indicating whether to standardize features before clustering (default: TRUE). Recommended for comparing features with different scales
#' @param method Character indicating preprocessing method: "none", "z-score" (default), or "min-max"
#' @param aggregate_by_factor Character indicating metadata factor to aggregate samples by before clustering (default: NULL, no aggregation). When provided, samples are averaged within each factor level
#' @details This function performs Fuzzy C-means clustering to group features by their abundance patterns across samples.
#' Features are assigned to clusters with membership degrees (0-1) indicating their strength of membership.
#' Requires the 'Mfuzz' package to be installed. If not present, the function will prompt the user to install it.
#' The data is standardized before clustering to ensure all features contribute equally regardless of their absolute abundance values.
#' @return Returns a romics_object with new columns in the statistics layer:
#' - Cmeans_cluster: The assigned cluster for each feature (highest membership)
#' - Cmeans_membership_[n]: Membership values (0-1) for each cluster
#' @author Geremy Clair
#' @export
#' @examples
#' # Basic Fuzzy C-means with 5 clusters
#' romics_object <- romicsFuzzyCMeans(romics_object, n_clusters = 5)
#'
#' # More fuzzy clustering with custom parameters
#' romics_object <- romicsFuzzyCMeans(
#'   romics_object,
#'   n_clusters = 4,
#'   m = 1.5,
#'   method = "z-score"
#' )
romicsFuzzyCMeans <- function(romics_object,
                             n_clusters = 3,
                             m = 1.25,
                             iter_max = 100,
                             standardize = TRUE,
                             method = "z-score",
                             aggregate_by_factor = NULL) {
  arguments <- as.list(match.call())

  # Input validation
  if(!is.romicsObject(romics_object) || missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }

  if(!is.numeric(n_clusters) || n_clusters < 2 || n_clusters != as.integer(n_clusters)) {
    stop("n_clusters must be an integer >= 2")
  }

  if(!is.numeric(m) || m <= 1) {
    stop("m (fuzziness parameter) must be numeric and > 1")
  }

  if(!is.numeric(iter_max) || iter_max < 1 || iter_max != as.integer(iter_max)) {
    stop("iter_max must be a positive integer")
  }

  if(!is.logical(standardize)) {
    stop("standardize must be TRUE or FALSE")
  }

  if(!method %in% c("none", "z-score", "min-max")) {
    stop("method must be one of: 'none', 'z-score', 'min-max'")
  }

  if(!is.null(aggregate_by_factor) && !is.character(aggregate_by_factor)) {
    stop("aggregate_by_factor must be NULL or a character string matching a metadata field")
  }

  # Check if Mfuzz is installed
  if(!requireNamespace("Mfuzz", quietly = TRUE)) {
    stop("The 'Mfuzz' package is required to run this function.\n",
         "Please install it using:\n",
         "  if (!requireNamespace('BiocManager', quietly = TRUE))\n",
         "    install.packages('BiocManager')\n",
         "  BiocManager::install('Mfuzz')")
  }

  message("Running Fuzzy C-means clustering...")
  message("  Number of clusters: ", n_clusters)
  message("  Fuzziness parameter (m): ", m)
  message("  Method: ", method)

  # Handle factor aggregation
  factor_map <- NULL
  if(!is.null(aggregate_by_factor)) {
    message("  Aggregating samples by factor: ", aggregate_by_factor)
    if(!aggregate_by_factor %in% rownames(romics_object$metadata)) {
      stop("Factor '", aggregate_by_factor, "' not found in metadata")
    }
    agg_result <- .aggregateRomicsByFactor(romics_object, aggregate_by_factor)
    data <- agg_result$data
    factor_map <- agg_result$factor_map
  } else {
    # Extract data
    data <- romics_object$data
  }
  n_features <- nrow(data)
  n_samples <- ncol(data)

  # Convert to matrix if needed
  if(!is.matrix(data)) {
    data <- as.matrix(data)
  }

  # Preprocess data
  if(standardize) {
    message("Preprocessing data with ", method, " standardization...")
    data_standardized <- switch(method,
      "z-score" = {
        # Z-score normalization: (x - mean) / sd
        means <- rowMeans(data, na.rm = TRUE)
        sds <- matrixStats::rowSds(data, na.rm = TRUE)
        # Handle zero standard deviation
        sds[sds == 0] <- 1
        (data - means) / sds
      },
      "min-max" = {
        # Min-max normalization: (x - min) / (max - min)
        mins <- matrixStats::rowMins(data, na.rm = TRUE)
        maxs <- matrixStats::rowMaxs(data, na.rm = TRUE)
        ranges <- maxs - mins
        # Handle zero range
        ranges[ranges == 0] <- 1
        (data - mins) / ranges
      },
      "none" = data
    )
  } else {
    data_standardized <- data
  }

  # Replace NAs with 0 for Mfuzz processing
  data_standardized[is.na(data_standardized)] <- 0

  # Create ExpressionSet for Mfuzz
  # Mfuzz expects features in rows, samples in columns
  eset <- Biobase::ExpressionSet(assayData = data_standardized)

  # Run Mfuzz clustering
  message("Performing Fuzzy C-means clustering...")
  tryCatch({
    # Mfuzz clustering
    cl <- Mfuzz::mfuzz(eset, c = n_clusters, m = m, iter.max = iter_max)

    # Extract cluster assignments and membership values
    cluster_assignments <- cl$cluster
    membership_matrix <- cl$membership

    # Initialize/check statistics layer
    if(is.null(romics_object$statistics)) {
      romics_object$statistics <- data.frame(matrix(nrow = nrow(data), ncol = 0))
      rownames(romics_object$statistics) <- rownames(data)
    }

    # Ensure stats has same rows as data
    if(nrow(romics_object$statistics) != nrow(data)) {
      warning("Statistics layer row count mismatch, recreating...")
      romics_object$statistics <- data.frame(matrix(nrow = nrow(data), ncol = 0))
      rownames(romics_object$statistics) <- rownames(data)
    }

    # Add cluster assignments
    romics_object$statistics$Cmeans_cluster <- as.character(cluster_assignments)

    # Add membership values for each cluster
    for(i in 1:ncol(membership_matrix)) {
      col_name <- paste0("Cmeans_membership_C", i)
      romics_object$statistics[[col_name]] <- membership_matrix[, i]
    }

    message("Clustering complete!")
    message("  Cluster assignments added to statistics layer as 'Cmeans_cluster'")
    message("  Membership values added as 'Cmeans_membership_C1', 'Cmeans_membership_C2', etc.")

    # Print cluster summary
    cluster_counts <- table(cluster_assignments)
    message("\nCluster summary:")
    for(i in 1:length(cluster_counts)) {
      message(paste0("  Cluster ", i, ": ", cluster_counts[i], " features"))
    }

  }, error = function(e) {
    stop("Error during Fuzzy C-means clustering: ", e$message)
  })

  # Update processing steps
  romics_object <- romicsUpdateSteps(romics_object, arguments)

  return(romics_object)
}


#' romicsClusteringElbowEval()
#' @description Evaluates the optimal number of clusters for feature clustering using the elbow method.
#' Supports both k-means and fuzzy c-means clustering evaluation.
#' @param romics_object A romics_object created using createRomicsObject()
#' @param method Character indicating clustering method: "kmeans" or "fuzzy_cmeans" (default: "kmeans")
#' @param max_clusters Numeric indicating maximum number of clusters to evaluate (default: 10)
#' @param standardize Logical indicating whether to standardize data before clustering (default: TRUE)
#' @param standardization Character indicating standardization method: "z-score", "min-max", or "none" (default: "z-score")
#' @param iter_max Numeric indicating maximum number of iterations for k-means (default: 100)
#' @param m Numeric fuzziness parameter for fuzzy c-means (default: 1.25, only used if method="fuzzy_cmeans")
#' @param nstart Numeric for k-means: number of random initializations (default: 25)
#' @param aggregate_by_factor Character indicating metadata factor to aggregate samples by before clustering (default: NULL, no aggregation). When provided, samples are averaged within each factor level
#' @details For k-means: Calculates within-cluster sum of squares (WCSS) for each cluster number.
#' For fuzzy c-means: Calculates within-cluster sum of squared distances weighted by membership degrees.
#' Lower values indicate tighter clusters. The "elbow" where the curve flattens suggests the optimal k.
#' @return Returns a ggplot object showing the elbow curve for cluster evaluation
#' @author Geremy Clair
#' @export
#' @examples
#' # Evaluate k-means clustering
#' plot <- romicsClusteringElbowEval(
#'   romics_object,
#'   method = "kmeans",
#'   max_clusters = 10
#' )
#' plot
#'
#' # Evaluate fuzzy c-means clustering
#' plot <- romicsClusteringElbowEval(
#'   romics_object,
#'   method = "fuzzy_cmeans",
#'   max_clusters = 8,
#'   m = 1.5
#' )
#' plot
romicsClusteringElbowEval <- function(romics_object,
                                      method = "kmeans",
                                      max_clusters = 10,
                                      standardize = TRUE,
                                      standardization = "z-score",
                                      iter_max = 100,
                                      m = 1.25,
                                      nstart = 25,
                                      aggregate_by_factor = NULL) {

  # Input validation
  if(!is.romicsObject(romics_object) || missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }

  if(!method %in% c("kmeans", "fuzzy_cmeans")) {
    stop("method must be either 'kmeans' or 'fuzzy_cmeans'")
  }

  if(!is.numeric(max_clusters) || max_clusters < 2 || max_clusters > nrow(romics_object$data)) {
    stop("max_clusters must be numeric between 2 and the number of features")
  }

  if(!is.logical(standardize)) {
    stop("standardize must be TRUE or FALSE")
  }

  if(!standardization %in% c("z-score", "min-max", "none")) {
    stop("standardization must be 'z-score', 'min-max', or 'none'")
  }

  if(!is.numeric(iter_max) || iter_max <= 0) {
    stop("iter_max must be a positive numeric value")
  }

  if(!is.numeric(m) || m <= 1) {
    stop("m must be a numeric value greater than 1")
  }

  if(!is.numeric(nstart) || nstart <= 0) {
    stop("nstart must be a positive numeric value")
  }

  if(!is.null(aggregate_by_factor) && !is.character(aggregate_by_factor)) {
    stop("aggregate_by_factor must be NULL or a character string matching a metadata field")
  }

  message("Evaluating optimal number of clusters using ", method, " method...")
  message("  Max clusters: ", max_clusters)
  message("  Standardization: ", standardization)

  # Handle factor aggregation
  if(!is.null(aggregate_by_factor)) {
    message("  Aggregating samples by factor: ", aggregate_by_factor)
    if(!aggregate_by_factor %in% rownames(romics_object$metadata)) {
      stop("Factor '", aggregate_by_factor, "' not found in metadata")
    }
    agg_result <- .aggregateRomicsByFactor(romics_object, aggregate_by_factor)
    data <- as.matrix(agg_result$data)
  } else {
    # Prepare data
    data <- as.matrix(romics_object$data)
  }

  # Standardize if requested
  if(standardize && standardization != "none") {
    if(standardization == "z-score") {
      data <- t(scale(t(data)))
    } else if(standardization == "min-max") {
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

  # Perform clustering evaluation
  if(method == "kmeans") {
    wcss <- vector(length = max_clusters)
    for(k in 1:max_clusters) {
      message("  Evaluating k=", k, "...")
      km <- stats::kmeans(t(data), centers = k, iter.max = iter_max, nstart = nstart)
      wcss[k] <- km$tot.withinss
    }
    plot_data <- data.frame(
      k = 1:max_clusters,
      value = wcss
    )
    y_label <- "Within-Cluster Sum of Squares (WCSS)"
    title <- "Elbow Method: K-Means Clustering Evaluation"

  } else if(method == "fuzzy_cmeans") {
    # Check if Mfuzz is available
    if(!requireNamespace("Mfuzz", quietly = TRUE)) {
      stop("The 'Mfuzz' package is required for fuzzy c-means. ",
           "Install it with: if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); ",
           "BiocManager::install('Mfuzz')")
    }

    wcss_fuzzy <- vector(length = max_clusters)
    wcss_fuzzy[1] <- NA  # c=1 is not meaningful for fuzzy c-means
    for(c_val in 2:max_clusters) {
      message("  Evaluating c=", c_val, "...")
      # Create ExpressionSet with proper row/column names
      eset <- Biobase::ExpressionSet(assayData = as.matrix(data))
      rownames(Biobase::exprs(eset)) <- rownames(data)
      colnames(Biobase::exprs(eset)) <- colnames(data)
      # Standardize if needed (data is already z-scored but Mfuzz expects raw data)
      # Skip standardization since we already standardized above
      # Run fuzzy c-means
      fc_result <- Mfuzz::mfuzz(eset, c = c_val, m = m, iter.max = iter_max)
      # Calculate within-cluster sum of squared distances
      # membership is features x clusters, centers is clusters x features
      sq_dists <- matrix(0, nrow(data), nrow(fc_result$centers))
      for(i in seq_len(nrow(data))) {
        for(c in seq_len(nrow(fc_result$centers))) {
          sq_dists[i, c] <- sum((data[i,] - fc_result$centers[c,])^2)
        }
      }
      wcss_fuzzy[c_val] <- sum((fc_result$membership^m) * sq_dists)
    }
    plot_data <- data.frame(
      k = 1:max_clusters,
      value = wcss_fuzzy
    )
    y_label <- "Weighted Within-Cluster Sum of Squares"
    title <- "Elbow Method: Fuzzy C-Means Clustering Evaluation"
  }

  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = k, y = value)) +
    ggplot2::geom_line(linewidth = 1, color = "steelblue") +
    ggplot2::geom_point(size = 3, color = "steelblue") +
    ggplot2::labs(
      title = title,
      x = "Number of Clusters",
      y = y_label
    ) +
    theme_ROP() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 11)
    )

  message("Evaluation complete!")
  return(p)
}

#' romicsClusteringSilhouetteEval()
#' @description Evaluates the optimal number of clusters for feature clustering using silhouette analysis.
#' Supports both k-means and fuzzy c-means clustering. Silhouette coefficient ranges from -1 to 1,
#' with values closer to 1 indicating well-separated clusters.
#' @param romics_object A romics_object created using createRomicsObject()
#' @param method Character indicating clustering method: "kmeans" or "fuzzy_cmeans" (default: "kmeans")
#' @param max_clusters Numeric indicating maximum number of clusters to evaluate (default: 10)
#' @param standardize Logical indicating whether to standardize data before clustering (default: TRUE)
#' @param standardization Character indicating standardization method: "z-score", "min-max", or "none" (default: "z-score")
#' @param iter_max Numeric indicating maximum number of iterations for k-means (default: 100)
#' @param m Numeric fuzziness parameter for fuzzy c-means (default: 1.25, only used if method="fuzzy_cmeans")
#' @param nstart Numeric for k-means: number of random initializations (default: 25)
#' @param aggregate_by_factor Character indicating metadata factor to aggregate samples by before clustering (default: NULL, no aggregation). When provided, samples are averaged within each factor level
#' @details Calculates the mean silhouette coefficient for each cluster number.
#' Higher silhouette values indicate better-defined clusters. The cluster number with the
#' highest average silhouette coefficient is typically optimal.
#' @return Returns a ggplot object showing silhouette scores for each cluster number
#' @author Geremy Clair
#' @export
#' @examples
#' # Evaluate k-means clustering with silhouette method
#' plot <- romicsClusteringSilhouetteEval(
#'   romics_object,
#'   method = "kmeans",
#'   max_clusters = 10
#' )
#' plot
#'
#' # Evaluate fuzzy c-means clustering
#' plot <- romicsClusteringSilhouetteEval(
#'   romics_object,
#'   method = "fuzzy_cmeans",
#'   max_clusters = 8
#' )
#' plot
romicsClusteringSilhouetteEval <- function(romics_object,
                                           method = "kmeans",
                                           max_clusters = 10,
                                           standardize = TRUE,
                                           standardization = "z-score",
                                           iter_max = 100,
                                           m = 1.25,
                                           nstart = 25,
                                           aggregate_by_factor = NULL) {

  # Input validation
  if(!is.romicsObject(romics_object) || missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }

  if(!method %in% c("kmeans", "fuzzy_cmeans")) {
    stop("method must be either 'kmeans' or 'fuzzy_cmeans'")
  }

  if(!is.numeric(max_clusters) || max_clusters < 2 || max_clusters > nrow(romics_object$data)) {
    stop("max_clusters must be numeric between 2 and the number of features")
  }

  if(!is.logical(standardize)) {
    stop("standardize must be TRUE or FALSE")
  }

  if(!standardization %in% c("z-score", "min-max", "none")) {
    stop("standardization must be 'z-score', 'min-max', or 'none'")
  }

  if(!is.numeric(iter_max) || iter_max <= 0) {
    stop("iter_max must be a positive numeric value")
  }

  if(!is.numeric(m) || m <= 1) {
    stop("m must be a numeric value greater than 1")
  }

  if(!is.numeric(nstart) || nstart <= 0) {
    stop("nstart must be a positive numeric value")
  }

  if(!is.null(aggregate_by_factor) && !is.character(aggregate_by_factor)) {
    stop("aggregate_by_factor must be NULL or a character string matching a metadata field")
  }

  message("Evaluating optimal number of clusters using silhouette method with ", method, "...")
  message("  Max clusters: ", max_clusters)
  message("  Standardization: ", standardization)

  # Handle factor aggregation
  if(!is.null(aggregate_by_factor)) {
    message("  Aggregating samples by factor: ", aggregate_by_factor)
    if(!aggregate_by_factor %in% rownames(romics_object$metadata)) {
      stop("Factor '", aggregate_by_factor, "' not found in metadata")
    }
    agg_result <- .aggregateRomicsByFactor(romics_object, aggregate_by_factor)
    data <- as.matrix(agg_result$data)
  } else {
    # Prepare data
    data <- as.matrix(romics_object$data)
  }

  # Standardize if requested
  if(standardize && standardization != "none") {
    if(standardization == "z-score") {
      data <- t(scale(t(data)))
    } else if(standardization == "min-max") {
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

  # Check for cluster package
  if(!requireNamespace("cluster", quietly = TRUE)) {
    stop("The 'cluster' package is required for silhouette analysis. Install it with: install.packages('cluster')")
  }

  # Perform clustering evaluation
  silhouette_scores <- vector(length = max_clusters)

  if(method == "kmeans") {
    silhouette_scores[1] <- NA  # k=1 is not meaningful for silhouette analysis
    for(k in 2:max_clusters) {
      message("  Evaluating k=", k, "...")
      km <- stats::kmeans(data, centers = k, iter.max = iter_max, nstart = nstart)
      sil <- cluster::silhouette(km$cluster, dist(data))
      silhouette_scores[k] <- mean(sil[, 3])
    }
    plot_data <- data.frame(
      k = 1:max_clusters,
      silhouette = silhouette_scores
    )
    title <- "Silhouette Analysis: K-Means Clustering Evaluation"

  } else if(method == "fuzzy_cmeans") {
    # Check if Mfuzz is available
    if(!requireNamespace("Mfuzz", quietly = TRUE)) {
      stop("The 'Mfuzz' package is required for fuzzy c-means. ",
           "Install it with: if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); ",
           "BiocManager::install('Mfuzz')")
    }

    silhouette_scores[1] <- NA  # c=1 is not meaningful for fuzzy c-means
    for(c_val in 2:max_clusters) {
      message("  Evaluating c=", c_val, "...")
      # Create ExpressionSet with proper row/column names
      eset <- Biobase::ExpressionSet(assayData = as.matrix(data))
      rownames(Biobase::exprs(eset)) <- rownames(data)
      colnames(Biobase::exprs(eset)) <- colnames(data)
      # Standardize if needed (data is already z-scored but Mfuzz expects raw data)
      # Skip standardization since we already standardized above
      # Run fuzzy c-means
      fc_result <- Mfuzz::mfuzz(eset, c = c_val, m = m, iter.max = iter_max)
      # Get hard cluster assignments from membership matrix (highest membership)
      clusters <- apply(fc_result$membership, 1, which.max)
      # Calculate silhouette score (use data not t(data) for feature-level clustering)
      sil <- cluster::silhouette(clusters, dist(data))
      silhouette_scores[c_val] <- mean(sil[, 3])
    }
    plot_data <- data.frame(
      k = 1:max_clusters,
      silhouette = silhouette_scores
    )
    title <- "Silhouette Analysis: Fuzzy C-Means Clustering Evaluation"
  }

  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = k, y = silhouette)) +
    ggplot2::geom_line(linewidth = 1, color = "darkgreen") +
    ggplot2::geom_point(size = 3, color = "darkgreen") +
    ggplot2::labs(
      title = title,
      x = "Number of Clusters",
      y = "Mean Silhouette Coefficient"
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
    ggplot2::ylim(c(-1, 1)) +
    theme_ROP() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 11)
    )

  message("Evaluation complete!")
  return(p)
}

