#' Helper function to find matching statistical columns
#' @param romics_object A romics_object with statistics
#' @param factor_levels Vector of factor levels to compare
#' @param pval_type Either "p" or "padj"
#' @return Data frame with combinations and their corresponding column names
findStatisticalColumns <- function(romics_object, factor_levels, pval_type) {
  if (is.null(romics_object$statistics) || ncol(romics_object$statistics) == 0) {
    return(data.frame(level1 = character(0), level2 = character(0),
                      column_name = character(0), test_name = character(0)))
  }

  stat_columns <- colnames(romics_object$statistics)

  # Generate all possible combinations
  by2combinations <- data.frame(t(combn(factor_levels, 2)))
  colnames(by2combinations) <- c("level1", "level2")

  results <- data.frame(level1 = character(0), level2 = character(0),
                        column_name = character(0), test_name = character(0))

  # For each combination, look for matching columns in both orders
  for (i in 1:nrow(by2combinations)) {
    level1 <- as.character(by2combinations[i, 1])
    level2 <- as.character(by2combinations[i, 2])

    # Create regex patterns for both orders
    # Pattern: level1_vs_level2_TESTNAME_pval_type
    pattern1 <- paste0("^", gsub("([.*+?^${}()|\\[\\]\\\\])", "\\\\\\1", level1),
                       "_vs_", gsub("([.*+?^${}()|\\[\\]\\\\])", "\\\\\\1", level2),
                       "_([^_]+)_", pval_type, "$")

    # Pattern: level2_vs_level1_TESTNAME_pval_type
    pattern2 <- paste0("^", gsub("([.*+?^${}()|\\[\\]\\\\])", "\\\\\\1", level2),
                       "_vs_", gsub("([.*+?^${}()|\\[\\]\\\\])", "\\\\\\1", level1),
                       "_([^_]+)_", pval_type, "$")

    # Check for matches in original order
    matches1 <- grep(pattern1, stat_columns, value = TRUE)
    if (length(matches1) > 0) {
      for (match in matches1) {
        test_name <- gsub(pattern1, "\\1", match)
        results <- rbind(results, data.frame(level1 = level1, level2 = level2,
                                             column_name = match, test_name = test_name))
      }
    }

    # Check for matches in reverse order
    matches2 <- grep(pattern2, stat_columns, value = TRUE)
    if (length(matches2) > 0) {
      for (match in matches2) {
        test_name <- gsub(pattern2, "\\1", match)
        results <- rbind(results, data.frame(level1 = level1, level2 = level2,
                                             column_name = match, test_name = test_name))
      }
    }
  }

  return(results)
}

#' singleFeaturePlot()
#' @description Creates a single feature plot from a romics_object with optional statistical comparison brackets. Supports multiple plot types including jitter, boxplot, violin, and combinations.
#' @param romics_object A romics_object created using romicsCreateObject().
#' @param feature Character. Name of the feature to plot. Must match exactly (case-sensitive). Default: "feature"
#' @param plot_type Character. Plot type: "jitter", "boxplot", "violin", "jb" (jitter+boxplot), or "jv" (jitter+violin). Default: "jb"
#' @param factor Character. Factor name from romics_object to use for grouping. Use "main" for the main factor. Default: "main"
#' @param ylim Numeric vector of length 2. Y-axis limits c(min, max). Default: NULL (auto-scaling)
#' @param title Character. Plot title. Use "auto" for feature name or provide custom title. Default: "auto"
#' @param pval Character. P-value type for statistical brackets: "none", "p", or "padj". Default: "none"
#' @param y_bracket_pos Numeric. Vertical spacing between statistical brackets. Default: 0.5
#' @param test_priority Character vector. Priority order for statistical tests when multiple are available. Default: c("Ttest", "ttest", "wilcoxon", "glm_binomial", "ANOVA")
#' @param significance_threshold Numeric. P-value threshold for highlighting significant results in red. Default: 0.05
#' @return A ggplot2 object
#' @details
#' The function automatically detects available statistical comparisons in the romics_object's statistics layer
#' and adds brackets showing p-values. When multiple statistical tests are available for the same comparison,
#' the function uses the test_priority order to select which test to display.
#'
#' Significant p-values (below significance_threshold) are displayed in red, while non-significant p-values
#' are displayed in black.
#'
#' @author Geremy Clair
#' @export
singleFeaturePlot <- function(romics_object, feature = "feature", plot_type = "jb", factor = "main",
                              ylim = NULL, title = "auto", pval = "none", y_bracket_pos = 0.5,
                              test_priority = c("Ttest", "ttest", "wilcoxon", "glm_binomial", "ANOVA"),
                              significance_threshold = 0.05) {
  # Load required libraries
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for this function")
  }
  library(ggplot2)

  # Input validation
  if (!is.romicsObject(romics_object) || missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }
  if (!plot_type %in% c("jitter", "boxplot", "violin", "jb", "jv")) {
    stop("plot_type has to be in 'jitter','boxplot','violin','jb',or 'jv'")
  }

  # Handle factor parameter
  if (missing(factor) || factor == "main") {
    factor <- romics_object$main_factor
  }
  if (!factor %in% romicsFactorNames(romics_object)) {
    cat("The chosen <factor> is not in the list of factors for the selected romics object.\n")
    cat("Available factors:\n")
    print(romicsFactorNames(romics_object))
    stop("Invalid factor specified")
  }

  # Validate other parameters
  if (length(pval) != 1 || !pval %in% c("none", "p", "padj")) {
    stop("<pval> has to be 'none','p',or 'padj'")
  }
  if (!is.numeric(y_bracket_pos) || y_bracket_pos <= 0 || length(y_bracket_pos) != 1) {
    stop("<y_bracket_pos> has to be a positive numerical vector of length 1")
  }
  if (!is.character(feature) || length(feature) != 1) {
    stop("<feature> should be a character vector of length 1")
  }
  if (!is.numeric(significance_threshold) || length(significance_threshold) != 1 ||
      significance_threshold <= 0 || significance_threshold >= 1) {
    stop("<significance_threshold> should be a numeric value between 0 and 1")
  }

  # Check for exact feature match
  if (!feature %in% rownames(romics_object$data)) {
    cat("Feature '", feature, "' not found in romics_object data.\n", sep = "")
    cat("Available features (first 10):\n")
    available_features <- rownames(romics_object$data)
    print(head(available_features, 10))
    if (length(available_features) > 10) {
      cat("... and", length(available_features) - 10, "more features\n")
    }
    cat("\nNote: Feature name must match exactly (case-sensitive)\n")
    stop("Feature not found")
  }

  # Update factor if needed
  if (factor != romics_object$main_factor) {
    romicsChangeFactor(romics_object, main_factor = factor)
  }

  # Extract data using exact feature name
  intensity <- as.numeric(romics_object$data[feature, ])
  group_vector <- as.character(romics_object$metadata[factor, ])
  fill_vector <- as.character(romics_object$metadata["colors_romics", ])

  # Create data frame with all samples
  data_all <- data.frame(
    intensity = intensity,
    group = group_vector,
    fill = fill_vector,
    sample = colnames(romics_object$data),
    stringsAsFactors = FALSE
  )

  # Get all unique factor levels
  all_levels <- unique(group_vector[!is.na(group_vector)])
  data_all$group <- factor(data_all$group, levels = all_levels)
  data_plot <- data_all[!is.na(data_all$intensity), ]

  # Get color mapping
  fill_colors <- setNames(unique(data_all$fill[!is.na(data_all$group)]),
                          unique(data_all$group[!is.na(data_all$group)]))

  # Create base plot
  plot <- ggplot(data_plot, aes(x = group, y = intensity)) +
    theme_ROP() +
    scale_x_discrete(drop = FALSE)

  if (!is.null(ylim) && is.numeric(ylim) && length(ylim) == 2) {
    plot <- plot + scale_y_continuous(limits = ylim)
  }

  # Add plot elements based on plot_type
  if (plot_type %in% c("boxplot", "jb")) {
    plot <- plot +
      geom_boxplot(aes(fill = group), alpha = 0.25) +
      scale_fill_manual(values = fill_colors, drop = FALSE)
  }

  if (plot_type %in% c("violin", "jv")) {
    plot <- plot +
      geom_violin(aes(fill = group), alpha = 0.25) +
      scale_fill_manual(values = fill_colors, drop = FALSE)
  }

  if (plot_type %in% c("jitter", "jb", "jv")) {
    plot <- plot +
      geom_jitter(aes(color = group), position = position_jitter(0.25), size = 3) +
      scale_color_manual(values = fill_colors, drop = FALSE)
  }

  if (plot_type %in% c("jitter", "violin", "jv")) {
    plot <- plot +
      stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1),
                   geom = "pointrange", color = "black", shape = 3)
  }

  # Add title
  if (title == "auto") {
    plot <- plot + ggtitle(feature)
  } else {
    plot <- plot + ggtitle(title)
  }

  # IMPROVED P-VALUE BRACKET GENERATION WITH SIGNIFICANCE HIGHLIGHTING
  if (pval != "none" && length(all_levels) > 1) {
    # Find all available statistical columns
    stat_matches <- findStatisticalColumns(romics_object, all_levels, pval)

    if (nrow(stat_matches) > 0) {
      # If multiple tests are available, prioritize based on test_priority
      if (length(unique(stat_matches$test_name)) > 1) {
        stat_matches$priority <- match(stat_matches$test_name, test_priority)
        stat_matches$priority[is.na(stat_matches$priority)] <- length(test_priority) + 1
        stat_matches <- stat_matches[order(stat_matches$priority), ]
        # Keep only the highest priority test for each comparison
        stat_matches <- stat_matches[!duplicated(paste(stat_matches$level1, stat_matches$level2)), ]
      }

      # Add brackets with significance-based coloring
      y_increment <- y_bracket_pos
      y_bracket_pos <- max(data_plot$intensity, na.rm = TRUE) + y_increment

      for (i in 1:nrow(stat_matches)) {
        pval_value <- romics_object$statistics[feature, stat_matches$column_name[i]]
        if (!is.na(pval_value)) {
          # Determine color based on significance
          bracket_color <- if (pval_value < significance_threshold) "red" else "black"

          plot <- plot +
            ggsignif::geom_signif(xmin = stat_matches$level1[i],
                                  xmax = stat_matches$level2[i],
                                  y_position = y_bracket_pos,
                                  annotations = formatC(as.numeric(pval_value), format = "e", digits = 2),
                                  textsize = 3,
                                  manual = TRUE)
          y_bracket_pos <- y_bracket_pos + y_increment
        }
      }

      # Print info about which tests were used
      used_tests <- unique(stat_matches$test_name)
      n_significant <- sum(sapply(stat_matches$column_name, function(col) {
        pval_val <- romics_object$statistics[feature, col]
        !is.na(pval_val) && pval_val < significance_threshold
      }))

      message(paste("P-value brackets added using:", paste(used_tests, collapse = ", ")))
      message(paste("Significant comparisons (p <", significance_threshold, "):", n_significant,
                    "out of", nrow(stat_matches), "(shown in red)"))
    } else {
      warning(paste0("No statistical comparisons found for factor <", factor,
                     "> with p-value type '", pval, "'. No brackets will be displayed."))
    }
  }

  return(plot)
}

#' multipleFeaturePlot()
#' @description Creates a grid of feature plots from a romics_object with optional statistical comparison brackets.
#' @param romics_object A romics_object created using romicsCreateObject().
#' @param features Character vector. Names of features to plot. Must match exactly (case-sensitive).
#' @param plot_type Character. Plot type: "jitter", "boxplot", "violin", "jb" (jitter+boxplot), or "jv" (jitter+violin). Default: "jb"
#' @param factor Character. Factor name from romics_object to use for grouping. Use "main" for the main factor. Default: "main"
#' @param ncol Integer. Number of columns in the plot grid. Default: NULL (automatic calculation)
#' @param limits Numeric vector of length 2 or "free". Y-axis limits c(min, max) or "free" for variable scales. Default: "free"
#' @param title Character. Overall plot title. Default: NULL (no title)
#' @param pval Character. P-value type for statistical brackets: "none", "p", or "padj". Default: "none"
#' @param y_bracket_pos Numeric. Vertical spacing between statistical brackets. Default: 0.5
#' @param test_priority Character vector. Priority order for statistical tests. Default: c("Ttest", "ttest", "wilcoxon", "glm_binomial", "ANOVA")
#' @param significance_threshold Numeric. P-value threshold for highlighting significant results in red. Default: 0.05
#' @param feature_titles Logical. Whether to display feature names as subplot titles. Default: TRUE
#' @param legend_position Character. Position of the legend: "right", "bottom", "none", etc. Default: "bottom"
#' @return A ggplot2 object with multiple feature plots arranged in a grid
#' @details
#' The function creates a grid of plots, each showing a different feature from the romics_object.
#' Statistical comparisons are added as brackets if requested. The grid layout can be customized.
#'
#' @author Brittney Gorman
#' @export
multipleFeaturePlot <- function(romics_object, features, plot_type = "jb", factor = "main",
                                ncol = NULL, limits = "free", title = NULL, pval = "none",
                                y_bracket_pos = 0.5, test_priority = c("Ttest", "ttest", "wilcoxon",
                                                                       "glm_binomial", "ANOVA"), significance_threshold = 0.05,
                                feature_titles = TRUE, legend_position = "bottom") {
  # Load required libraries
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for this function")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("patchwork package is required for this function")
  }
  library(ggplot2)
  library(patchwork)

  # Input validation
  if (!is.romicsObject(romics_object) || missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }

  if (missing(features) || !is.character(features) || length(features) < 1) {
    stop("features must be a character vector of feature names with at least one element")
  }

  # Check if all features exist in the romics_object
  missing_features <- features[!features %in% rownames(romics_object$data)]
  if (length(missing_features) > 0) {
    cat("The following features were not found in the romics_object:\n")
    print(missing_features)
    cat("Available features (first 10):\n")
    available_features <- rownames(romics_object$data)
    print(head(available_features, 10))
    if (length(available_features) > 10) {
      cat("... and", length(available_features) - 10, "more features\n")
    }
    stop("One or more features not found")
  }

  # Determine number of rows and columns for the grid
  n_features <- length(features)
  if (is.null(ncol)) {
    # Auto calculate columns: square root rounded to nearest integer
    ncol <- round(sqrt(n_features))
  }
  nrow <- ceiling(n_features / ncol)

  # Create a list to store individual plots
  plot_list <- list()

  # Generate each individual plot
  for (i in seq_along(features)) {
    feature <- features[i]

    # Create plot title based on feature name if requested
    plot_title <- if (feature_titles) feature else ""

    # Create single feature plot
    single_plot <- singleFeaturePlot(
      romics_object = romics_object,
      feature = feature,
      plot_type = plot_type,
      factor = factor,
      ylim = if (limits != "free") limits else NULL,
      title = plot_title,
      pval = pval,
      y_bracket_pos = y_bracket_pos,
      test_priority = test_priority,
      significance_threshold = significance_threshold
    )

    # Remove legend from individual plots (will add a combined one later)
    single_plot <- single_plot + theme(legend.position = "none")

    # Store the plot
    plot_list[[i]] <- single_plot
  }

  # Combine plots using patchwork
  combined_plot <- wrap_plots(plot_list, ncol = ncol)

  # Add a shared legend if requested
  if (legend_position != "none") {
    # Extract legend from a temporary plot
    temp_plot <- singleFeaturePlot(
      romics_object = romics_object,
      feature = features[1],
      plot_type = plot_type,
      factor = factor
    )

    # Apply shared legend and customize layout
    combined_plot <- combined_plot +
      plot_layout(guides = "collect") +
      plot_annotation(
        title = title,
        theme = theme(
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          legend.position = legend_position
        )
      )
  }

  # Return the combined plot
  return(combined_plot)
}


#' multipleFeatureComparisonPlot()
#' @description Creates a single plot comparing multiple features from a romics_object across groups.
#' @param romics_object A romics_object created using romicsCreateObject().
#' @param features Character vector. Names of features to plot. Must match exactly (case-sensitive).
#' @param point_type Character. How to represent individual data points: "scatter", "boxplot", "violin", "bar". Default: "scatter"
#' @param factor Character. Factor name from romics_object to use for grouping. Use "main" for the main factor. Default: "main"
#' @param group_by Character. How to organize the plot: "feature" (features together) or "factor" (factors together). Default: "feature"
#' @param normalize Logical or character. Whether/how to normalize feature values:
#'        FALSE (no normalization), "zscore", "minmax", or "percent_max". Default: FALSE
#' @param group_spacing Numeric. Space between groups (features or factors). Default: 0.5
#' @param point_size Numeric. Size of points for scatter representation. Default: 2
#' @param element_width Numeric. Width of boxplots/violins. Default: 0.7
#' @param element_alpha Numeric. Transparency for elements (0-1). Default: 0.7
#' @param title Character. Plot title. Default: "Feature Comparison"
#' @param feature_labels Character vector. Custom labels for features (same length as features). Default: NULL (use feature names)
#' @param x_axis_label Character. Label for the x-axis. Default: NULL (auto-set based on group_by)
#' @param y_axis_label Character. Label for the y-axis. Default: "Intensity" (changes based on normalization)
#' @param legend_title Character. Title for the legend. Default: NULL (auto-set based on group_by)
#' @param palette Character vector. Custom color palette for features/groups. Default: NULL (uses built-in palette)
#' @param show_dividers Logical. Whether to show dividing lines between groups. Default: TRUE
#' @param divider_color Character. Color for the dividing lines. Default: "gray70"
#' @param divider_linetype Character. Line type for dividers: "solid", "dashed", etc. Default: "dashed"
#' @return A ggplot2 object
#' @details
#' This function allows plotting multiple features within the same plot for direct comparison.
#' Features can be normalized to enable comparison on the same scale. Various representation types are
#' supported for individual data points, including scatter, boxplots, violins, and bars.
#'
#' @author Your Name
#' @export
multipleFeatureComparisonPlot <- function(romics_object, features, point_type = "scatter", factor = "main",
                                       group_by = "feature", normalize = FALSE,
                                       group_spacing = 0.5, point_size = 2, element_width = 0.7,
                                       element_alpha = 0.7, title = "Feature Comparison", feature_labels = NULL,
                                       x_axis_label = NULL, y_axis_label = "Intensity", legend_title = NULL,
                                       palette = NULL, show_dividers = TRUE, divider_color = "gray70",
                                       divider_linetype = "dashed",group_order = NULL) {
  # Load required libraries
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for this function")
  }
  library(ggplot2)
  library(reshape2)  # For data reshaping

  # Input validation
  if (!is.romicsObject(romics_object) || missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }

  if (missing(features) || !is.character(features) || length(features) < 1) {
    stop("features must be a character vector of feature names with at least one element")
  }

  if (!point_type %in% c("scatter", "boxplot", "violin", "bar")) {
    stop("point_type must be one of: 'scatter', 'boxplot', 'violin', 'bar'")
  }

  if (!group_by %in% c("feature", "factor")) {
    stop("group_by must be either 'feature' or 'factor'")
  }

  if (!is.logical(normalize) && !normalize %in% c("zscore", "minmax", "percent_max")) {
    stop("normalize must be FALSE or one of: 'zscore', 'minmax', 'percent_max'")
  }

  # Check if all features exist in the romics_object
  missing_features <- features[!features %in% rownames(romics_object$data)]
  if (length(missing_features) > 0) {
    cat("The following features were not found in the romics_object:\n")
    print(missing_features)
    stop("One or more features not found")
  }

  # Handle factor parameter
  if (missing(factor) || factor == "main") {
    factor <- romics_object$main_factor
  }

  if (!factor %in% romicsFactorNames(romics_object)) {
    cat("The chosen <factor> is not in the list of factors for the selected romics object.\n")
    cat("Available factors:\n")
    print(romicsFactorNames(romics_object))
    stop("Invalid factor specified")
  }

  # Extract data for the selected features
  feature_data <- romics_object$data[features, , drop = FALSE]

  # Extract group information
  group_vector <- as.character(romics_object$metadata[factor, ])
  group_levels <- unique(group_vector[!is.na(group_vector)])

  # Get unique group levels
  group_levels <- unique(group_vector[!is.na(group_vector)])

  # Apply custom group ordering if provided
  if (!is.null(group_order)) {
    # Check if all groups are in the ordering
    missing_groups <- setdiff(group_levels, group_order)
    extra_groups <- setdiff(group_order, group_levels)

    if (length(missing_groups) > 0) {
      warning("Some groups are not included in group_order and will be added at the end: ",
              paste(missing_groups, collapse = ", "))
      group_order <- c(group_order, missing_groups)
    }

    if (length(extra_groups) > 0) {
      warning("Some values in group_order don't exist in the data: ",
              paste(extra_groups, collapse = ", "))
      group_order <- intersect(group_order, c(group_levels, extra_groups))
    }

    # Apply the custom order
    group_levels <- group_order
  }


  # Create a data frame for plotting
  plot_data <- data.frame(
    sample = colnames(romics_object$data),
    group = group_vector,
    stringsAsFactors = FALSE
  )

  # Add feature intensities to the data frame
  for (feature in features) {
    plot_data[[feature]] <- as.numeric(romics_object$data[feature, ])
  }

  # Convert to long format for ggplot2
  plot_data_long <- reshape2::melt(
    plot_data,
    id.vars = c("sample", "group"),
    measure.vars = features,
    variable.name = "feature",
    value.name = "intensity"
  )

  # Remove NA values
  plot_data_long <- plot_data_long[!is.na(plot_data_long$intensity), ]

  # Apply normalization if requested
  if (normalize != FALSE) {
    # Update y-axis label based on normalization method
    if (normalize == "zscore") {
      y_axis_label <- "Z-score"
      # Apply z-score normalization by feature
      plot_data_long <- do.call(rbind, lapply(split(plot_data_long, plot_data_long$feature), function(df) {
        mean_val <- mean(df$intensity, na.rm = TRUE)
        sd_val <- sd(df$intensity, na.rm = TRUE)
        if (sd_val > 0) {
          df$intensity <- (df$intensity - mean_val) / sd_val
        }
        return(df)
      }))
    } else if (normalize == "minmax") {
      y_axis_label <- "Min-Max Normalized Intensity"
      # Apply min-max normalization by feature
      plot_data_long <- do.call(rbind, lapply(split(plot_data_long, plot_data_long$feature), function(df) {
        min_val <- min(df$intensity, na.rm = TRUE)
        max_val <- max(df$intensity, na.rm = TRUE)
        if (max_val > min_val) {
          df$intensity <- (df$intensity - min_val) / (max_val - min_val)
        }
        return(df)
      }))
    } else if (normalize == "percent_max") {
      y_axis_label <- "% of Maximum"
      # Apply percent of maximum normalization by feature
      plot_data_long <- do.call(rbind, lapply(split(plot_data_long, plot_data_long$feature), function(df) {
        max_val <- max(df$intensity, na.rm = TRUE)
        if (max_val > 0) {
          df$intensity <- (df$intensity / max_val) * 100
        }
        return(df)
      }))
    }
  }

  # Apply custom feature labels if provided
  if (!is.null(feature_labels)) {
    if (length(feature_labels) != length(features)) {
      warning("feature_labels length doesn't match features length. Using original feature names.")
    } else {
      # Create mapping from original feature names to custom labels
      feature_mapping <- setNames(feature_labels, features)
      # Replace feature levels with labels
      plot_data_long$feature <- factor(plot_data_long$feature,
                                       levels = features,
                                       labels = feature_labels)
    }
  } else {
    # Ensure feature is a factor with proper ordering
    plot_data_long$feature <- factor(plot_data_long$feature, levels = features)
  }

  # Ensure group is a factor with proper ordering
  plot_data_long$group <- factor(plot_data_long$group, levels = group_levels)

  # Define primary and secondary grouping based on group_by parameter
  if (group_by == "feature") {
    primary_group <- "feature"
    secondary_group <- "group"
    if (is.null(x_axis_label)) x_axis_label <- "Features"
    if (is.null(legend_title)) legend_title <- "Group"
  } else { # group_by == "factor"
    primary_group <- "group"
    secondary_group <- "feature"
    if (is.null(x_axis_label)) x_axis_label <- "Groups"
    if (is.null(legend_title)) legend_title <- "Feature"
  }

  # Calculate x positions for plotting
  primary_levels <- levels(plot_data_long[[primary_group]])
  secondary_levels <- levels(plot_data_long[[secondary_group]])

  n_primary <- length(primary_levels)
  n_secondary <- length(secondary_levels)

  # Create position map for x-axis positioning
  position_df <- expand.grid(
    primary = seq_len(n_primary),
    secondary = seq_len(n_secondary)
  )

  # Calculate x positions with spacing between primary groups
  base_width <- element_width / n_secondary
  half_group_width <- (n_secondary * base_width) / 2

  position_df$x_pos <- ((position_df$primary - 1) * (1 + group_spacing)) +
    ((position_df$secondary - 0.5) * base_width) - half_group_width + 0.5

  # Map these positions to the actual primary and secondary names
  primary_map <- data.frame(
    value = primary_levels,
    num = seq_along(primary_levels),
    stringsAsFactors = FALSE
  )

  secondary_map <- data.frame(
    value = secondary_levels,
    num = seq_along(secondary_levels),
    stringsAsFactors = FALSE
  )

  # Add numeric indices to the plotting data
  colnames(primary_map)[1] <- primary_group
  colnames(secondary_map)[1] <- secondary_group

  plot_data_long$primary_num <- match(plot_data_long[[primary_group]], primary_map[[primary_group]])
  plot_data_long$secondary_num <- match(plot_data_long[[secondary_group]], secondary_map[[secondary_group]])

  # Merge to get x positions
  plot_data_pos <- merge(plot_data_long,
                         position_df,
                         by.x = c("primary_num", "secondary_num"),
                         by.y = c("primary", "secondary"))

  # Set up color palette
  color_group <- if (group_by == "feature") "group" else "feature"
  color_levels <- levels(plot_data_long[[color_group]])
  n_colors <- length(color_levels)

  if (is.null(palette)) {
    # Default color palettes based on number of groups
    if (n_colors <= 8) {
      palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                   "#0072B2", "#D55E00", "#CC79A7", "#999999")[1:n_colors]
    } else {
      palette <- colorRampPalette(
        c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
      )(n_colors)
    }
  } else if (length(palette) < n_colors) {
    warning("Provided palette has fewer colors than needed. Recycling colors.")
    palette <- rep(palette, length.out = n_colors)
  }

  # Create the tick positions and labels for x-axis
  tick_positions <- sapply(1:n_primary, function(p) {
    ((p - 1) * (1 + group_spacing)) + 0.5
  })

  # Create the base plot with common aesthetics
  p <- ggplot(plot_data_pos) +
    theme_bw() +
    labs(
      title = title,
      x = x_axis_label,
      y = y_axis_label
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = if(n_primary > 6) 45 else 0, hjust = if(n_primary > 6) 1 else 0.5),
      legend.position = "right"
    ) +
    scale_x_continuous(
      breaks = tick_positions,
      labels = primary_levels,
      expand = c(0.05, 0.05)
    )

  # Add representation based on point_type
  if (point_type == "scatter") {
    # For scatter plots, we use direct mapping without aes_string
    if (color_group == "group") {
      p <- p + geom_point(
        aes(x = x_pos, y = intensity, color = group),
        data = plot_data_pos,
        position = position_jitterdodge(
          jitter.width = element_width / (4 * n_secondary),
          dodge.width = element_width / n_secondary
        ),
        size = point_size, alpha = element_alpha
      )
    } else { # color_group == "feature"
      p <- p + geom_point(
        aes(x = x_pos, y = intensity, color = feature),
        data = plot_data_pos,
        position = position_jitterdodge(
          jitter.width = element_width / (4 * n_secondary),
          dodge.width = element_width / n_secondary
        ),
        size = point_size, alpha = element_alpha
      )
    }
    p <- p + scale_color_manual(values = palette, name = legend_title)

  } else if (point_type == "boxplot") {
    # Create a grouping variable for each unique combination
    plot_data_pos$grouping <- interaction(plot_data_pos[[primary_group]],
                                          plot_data_pos[[secondary_group]],
                                          drop = TRUE,
                                          sep = "_")

    # For boxplots, use the group aesthetic directly
    if (color_group == "group") {
      p <- p + geom_boxplot(
        aes(x = x_pos, y = intensity, fill = group, group = grouping),
        data = plot_data_pos,

        varwidth = FALSE,
        alpha = element_alpha,
        outlier.size = point_size * 0.7,
        outlier.alpha = element_alpha
      )
    } else { # color_group == "feature"
      p <- p + geom_boxplot(
        aes(x = x_pos, y = intensity, fill = feature, group = grouping),
        data = plot_data_pos,

        varwidth = FALSE,
        alpha = element_alpha,
        outlier.size = point_size * 0.7,
        outlier.alpha = element_alpha
      )
    }
    p <- p + scale_fill_manual(values = palette, name = legend_title)

  } else if (point_type == "violin") {
    # Create a grouping variable for each unique combination
    plot_data_pos$grouping <- interaction(plot_data_pos[[primary_group]],
                                          plot_data_pos[[secondary_group]],
                                          drop = TRUE,
                                          sep = "_")

    # For violins, use the group aesthetic directly
    if (color_group == "group") {
      p <- p + geom_violin(
        aes(x = x_pos, y = intensity, fill = group, group = grouping),
        data = plot_data_pos,
        width = element_width / n_secondary,
        alpha = element_alpha,
        scale = "width",
        trim = TRUE
      )
    } else { # color_group == "feature"
      p <- p + geom_violin(
        aes(x = x_pos, y = intensity, fill = feature, group = grouping),
        data = plot_data_pos,
        width = element_width / n_secondary,
        alpha = element_alpha,
        scale = "width",
        trim = TRUE
      )
    }
    p <- p + scale_fill_manual(values = palette, name = legend_title)

  } else if (point_type == "bar") {
    # Bar plot implementation with direct mapping instead of aes_string
    formula_str <- paste("intensity ~", primary_group, "+", secondary_group)
    summary_formula <- as.formula(formula_str)

    summary_data <- aggregate(summary_formula, data = plot_data_long,
                              FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                                  se = sd(x, na.rm = TRUE)/sqrt(length(x))))
    summary_data <- do.call(data.frame, summary_data)
    colnames(summary_data)[ncol(summary_data)-1:0] <- c("mean", "se")

    summary_data$primary_num <- match(summary_data[[primary_group]], primary_map[[primary_group]])
    summary_data$secondary_num <- match(summary_data[[secondary_group]], secondary_map[[secondary_group]])

    summary_pos <- merge(summary_data,
                         position_df,
                         by.x = c("primary_num", "secondary_num"),
                         by.y = c("primary", "secondary"))

    # Create a grouping variable for bars
    summary_pos$grouping <- interaction(summary_pos[[primary_group]],
                                        summary_pos[[secondary_group]],
                                        drop = TRUE,
                                        sep = "_")

    if (color_group == "group") {
      p <- p + geom_bar(
        aes(x = x_pos, y = mean, fill = group),
        data = summary_pos,
        stat = "identity",
        width = element_width / n_secondary,
        varwidth = FALSE,
        alpha = element_alpha
      )
    } else { # color_group == "feature"
      p <- p + geom_bar(
        aes(x = x_pos, y = mean, fill = feature),
        data = summary_pos,
        stat = "identity",
        width = element_width / n_secondary,
        varwidth = FALSE,
        alpha = element_alpha
      )
    }

    p <- p + geom_errorbar(
      aes(x = x_pos, ymin = mean - se, ymax = mean + se),
      data = summary_pos,
      width = (element_width / n_secondary) * 0.4,
      color = "black"
    )

    p <- p + scale_fill_manual(values = palette, name = legend_title)
  }

  # Add vertical dividers between primary groups
  if (show_dividers && n_primary > 1) {
    separator_positions <- sapply(1:(n_primary-1), function(i) {
      i * (1 + group_spacing) - group_spacing/2
    })

    p <- p + geom_vline(xintercept = separator_positions,
                        linetype = divider_linetype,
                        color = divider_color,
                        size = 0.5)
  }

  return(p)
}


#' bubblePlotFeatures()
#' @description will generate a bubble plot from a list of features the features will be clustered using a hclust, mean values will be represented as colors and percentage of completeness as dot size.
#' @param romics_object A romics_object created using romicsCreateObject()
#' @param factor has to be either 'main" or a factor of the romics_object, the list of factor can be retrieved using the function romicsFactorNames()
#' @param feature_list A list of features contained in the Romics_object can be found using the function romicsFeatureNames()
#' @param scale_feature has to be TRUE or FALSE to indicate if the calculated means have to be scaled by feature or not
#' @param gradient_colors has to be a vector of lenght 3 containing three colors, they will be used for the low, med, and high values of the gradient respectively.
#' @param cluster_features TRUE or FALSE to indicate if features should be hierarchically clustered and displayed with a dendrogram (default: TRUE)
#' @details This function will generate a ggplot style figure.
#' @return A ggplot styled figure
#' @author Geremy Clair
#' @export
bubblePlotFeatures<-function(romics_object, factor="main", scale_feature=T, feature_list=c("a","b"), gradient_colors=c("#fff1ed","#fc896a","#69010e"), cluster_features=TRUE){
    arguments<-as.list(match.call())

    if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("<romics_object> is missing or is not in the appropriate format.")}
    if(missing(factor)){factor="main"}
    if(!factor %in% c("main",romicsFactorNames(romics_object))){stop("The selected <factor> is not in the list of factors of the <romics_object>.")}
    if(!scale_feature %in% c(TRUE,FALSE)){stop("<scale_feature> has to be TRUE or FALSE")}
    if(!cluster_features %in% c(TRUE,FALSE)){stop("<cluster_features> has to be TRUE or FALSE")}
    if(missing(feature_list)){feature_list<-rownames(romics_object$data)}
    if(factor !="main"){
    romics_object<-romicsChangeFactor(romics_object, main_factor = factor)
    }

    romics_object<- romics_object[!names(romics_object)=="statistics"]
    class(romics_object)<-"romics_object"
    romics_object<-romicsMean(romics_object)
    romics_object<-romicsPercentComplete(romics_object)

    df<-romics_object$statistics
    # Preserve original feature_list order for later use when cluster_features=FALSE
    original_feature_order <- feature_list
    feature_list<-feature_list[feature_list %in% rownames(df)]
    df<-df[rownames(df) %in% feature_list,]
    df$Feature<- rownames(df)
    rownames(df)<-NULL

    df_mean<-df[grepl("_mean|Feature",colnames(df))]
    if(scale_feature==T){
      df_mean[,grepl("_mean",colnames(df_mean))]<-t(scale(t(df_mean[,grepl("_mean",colnames(df_mean))])))
    }
    colnames(df_mean)<-gsub("_mean","",colnames(df_mean))

    df_completeness<-df[grepl("_percentage_completeness|Feature",colnames(df))]
    colnames(df_completeness)<-gsub("_percentage_completeness","",colnames(df_completeness))

    #groups<-gsub("_mean","",colnames(df)[grepl("_mean",colnames(df))])

    d<-data.frame(matrix(ncol=4,nrow=nrow(df)))
    colnames(d)<-c("Feature","Group","Mean","Percentage_completeness")

    rownames(df)<-df$Feature
    mat<-as.matrix(df[,colnames(df)!="Feature"]  )

    groups<-colnames(df_mean)[colnames(df_mean)!="Feature"]

    for (i in 1:length(groups)){
      d$Feature<-df$Feature
      d$Group<-groups[i]
      d$Mean<-df_mean[,colnames(df_mean)==groups[i]]
      d$Percentage_completeness<-df_completeness[,colnames(df_completeness)==groups[i]]
      if(i==1){d_final<-d}else(d_final<-rbind(d_final,d))
    }

    d_final<-as.data.frame(d_final)

    # Perform clustering if requested
    if(cluster_features) {
      clust<- hclust(dist(mat))
      dend<-ggtree::ggtree(clust,hang=-1,)
      # Order features by clustering
      d_final$Feature<-factor(d_final$Feature, levels=clust$labels[clust$order])
    } else {
      # Order features by input order (filter to only features that exist in the data)
      input_order <- original_feature_order[original_feature_order %in% unique(d_final$Feature)]
      d_final$Feature<-factor(d_final$Feature, levels=input_order)
    }

    p<-ggplot(d_final,aes(y=Feature,x=Group,size = Percentage_completeness,color=Mean))+
      geom_point()+
      theme_ROP()

    if(missing(gradient_colors)){
      if(scale_feature==T){p<-p+scale_color_viridis_c(name = 'scaled (mean in group)')}else{
        p<-p+scale_color_viridis_c(name = "mean abundance")}
    }else{
      if(scale_feature==T){p<-p+scale_color_gradient2(low= gradient_colors[1],mid = gradient_colors[2],high = gradient_colors[3], name = 'scaled (mean in group)')}else{
        p<-p+scale_color_viridis_c(name = "mean abundance")}
    }

    p<-p+ theme(axis.line  = element_blank()) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      ylab('') + theme(axis.ticks = element_blank()) +scale_y_discrete(position = "right")

    # Add dendrogram only if clustering is enabled
    if(cluster_features) {
      p<-cowplot::plot_grid(dend, p, nrow = 1, rel_widths = c(0.5,2), align = 'h')
    }
    p
  }
