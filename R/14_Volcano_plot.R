#' romicsVolcano()
#' @description generate one or multiple volcano plots from t.test, wilcox.test, or LMM run using the functions romicsTtest(), romicsWilcoxTest(), or romicsLMM() respectively
#' @param p_type 'p' or 'padj' to indicate the type of pvalue to consider for the volcano plots.
#' @param p numeric value to indicate the maximum pvalue to use for the coloring of the volcano plot
#' @param colors numeric vector of lenght 3 used to color the features lower, not significantly changing and higher in the considered comparisons.
#' @param stat_type 't.test', 'wilcox.test', 'LMM' or 'auto' to indicate what statistics to use for the volcano plots generation. If 'auto', will detect available tests.
#' @param plot either 'all' if all paired comparisons have to be displayed OR a vector of numeric values comprised between 1 and the maximum number of possible plots generated.
#' @param plotly logical (TRUE or FALSE) to indicate whether to return interactive plotly plots (TRUE) or static ggplot plots (FALSE).
#' @param label_features logical (TRUE or FALSE) to indicate whether to label top features on the plot.
#' @param top_features numeric value indicating the number of top upregulated and downregulated features to label (only those passing the filters).
#' @param top_features_by 'abundance' or 'p' to indicate the criteria for selecting top features. 'abundance' selects by fold change magnitude, 'p' selects by p-value significance.
#' @param xlim numeric vector of length 2 specifying the x-axis limits (e.g., c(-2, 2)). If NULL, limits are determined automatically.
#' @param ylim numeric vector of length 2 specifying the y-axis limits (e.g., c(0, 10)), or a single numeric value to set max y-limit (min is set to 0). If NULL and auto_ylim=TRUE, limits are set based on data max.
#' @param auto_ylim logical (TRUE or FALSE) to automatically set ylim based on the maximum -log10(p) value in the data (default: TRUE). Ignored if ylim is manually specified.
#' @param size numeric value for point size in the plots. Default is 2.
#' @param alpha numeric value for point opacity in the plots (0-1). Default is 0.5.
#' @details generate one or multiple volcano plots from t.test, wilcox.test, or LMM run using the functions romicsTtest(), romicsWilcoxTest(), or romicsLMM() respectively. Automatically detects whether data is clustered or not. If xlim or ylim are specified, a warning is issued if any data points fall outside these limits.
#' @return This function will print different plots requested
#' @author Geremy Clair
#' @export
romicsVolcano <- function(romics_object, p_type = "p", p = 0.05, min_fold_change = 0.6,
                          colors = c("#2cbcb2", "#242021", "#d44e28"),
                          stat_type = "auto", plot = "all", plotly = FALSE,
                          label_features = FALSE, top_features = 10, top_features_by = "abundance",
                          xlim = NULL, ylim = NULL, auto_ylim = TRUE, size = 2, alpha = 0.5) {

  if(!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }
  if(missing(p_type)){p_type = "p"}
  if(!p_type %in% c("p", "padj")){
    stop("'p_type' has to be either 'p' or 'padj'")
  }
  if(missing(colors)){
    colors = c("#2cbcb2", "#242021", "#d44e28")
  }
  if(!is.character(colors) | length(colors) != 3){
    warning("'colors' should be a color character vector of lenght 3. The defaults colors were used")
    colors = c("#2cbcb2", "#242021", "#d44e28")
  }
  if(missing(p)){p = 0.05}
  if(!is.numeric(p) | p > 1 | p < 0){
    stop("'p' should be numeric and comprised between 0 and 1.")
  }
  if(missing(min_fold_change)){min_fold_change = 0.6}
  if(!is.numeric(min_fold_change) | min_fold_change < 0){
    stop("'min_fold_change' should be numeric and comprised higher than 0.")
  }
  if(missing(plotly)){plotly = FALSE}
  if(!is.logical(plotly)){
    warning("'plotly' was not logical (TRUE/FALSE). FALSE was used by default (ggplot output).")
    plotly = FALSE
  }
  if(missing(label_features)){label_features = FALSE}
  if(!is.logical(label_features)){
    warning("'label_features' was not logical (TRUE/FALSE). FALSE was used by default.")
    label_features = FALSE
  }
  if(missing(top_features)){top_features = 10}
  if(!is.numeric(top_features) | top_features < 0 | top_features != round(top_features)){
    warning("'top_features' should be a positive integer. 10 was used by default.")
    top_features = 10
  }
  if(missing(top_features_by)){top_features_by = "abundance"}
  if(!top_features_by %in% c("abundance", "p")){
    warning("'top_features_by' should be either 'abundance' or 'p'. 'abundance' was used by default.")
    top_features_by = "abundance"
  }
  if(missing(plot)){plot = "all"}
  if(missing(xlim)){xlim = NULL}
  if(!is.null(xlim)){
    if(!is.numeric(xlim) | length(xlim) != 2){
      warning("'xlim' should be a numeric vector of length 2 (e.g., c(-2, 2)). NULL was used by default.")
      xlim = NULL
    } else if(xlim[1] >= xlim[2]){
      warning("'xlim' should have xlim[1] < xlim[2]. NULL was used by default.")
      xlim = NULL
    }
  }
  if(missing(ylim)){ylim = NULL}
  if(!is.null(ylim)){
    if(is.numeric(ylim) & length(ylim) == 1){
      # Allow single numeric value as max y-limit
      ylim <- c(0, ylim)
    } else if(!is.numeric(ylim) | length(ylim) != 2){
      warning("'ylim' should be a numeric vector of length 2 (e.g., c(0, 10)) or a single numeric value for max. NULL was used by default.")
      ylim = NULL
    } else if(ylim[1] >= ylim[2]){
      warning("'ylim' should have ylim[1] < ylim[2]. NULL was used by default.")
      ylim = NULL
    }
  }
  if(missing(auto_ylim)){auto_ylim = TRUE}
  if(!is.logical(auto_ylim)){
    warning("'auto_ylim' was not logical (TRUE/FALSE). TRUE was used by default.")
    auto_ylim = TRUE
  }

  # Extract stats
  stats <- romics_object$statistics
  # Remove duplicate columns
  stats <- stats[, unique(colnames(stats))]

  # Auto-detect stat_type if not specified
  if(stat_type == "auto"){
    if(sum(grepl("_Ttest_p", colnames(stats))) > 0) {
      stat_type = "t.test"
      print("'stat_type' was set to 't.test' (auto-detected)")
    } else if(sum(grepl("_Wilcox_test_p", colnames(stats))) > 0) {
      stat_type = "wilcox.test"
      print("'stat_type' was set to 'wilcox.test' (auto-detected)")
    } else if(sum(grepl("_LMMtest_p", colnames(stats))) > 0) {
      stat_type = "LMM"
      print("'stat_type' was set to 'LMM' (auto-detected)")
    } else {
      stop("No statistical tests found in the romics_object")
    }
  }

  # Verify if test type exists
  if(sum(grepl("_Ttest_p", colnames(stats))) +
     sum(grepl("_Wilcox_test_p", colnames(stats))) +
     sum(grepl("_LMMtest_p", colnames(stats))) <= 0) {
    stop("Prior to plot the volcano plot(s) either T.tests, Wilcox.tests, or LMM tests have to be run")
  }

  # If stat_type missing use t.test by default unless does not exist
  if(missing(stat_type)){
    stat_type = "t.test"
    print("'stat_type' was missing 't.test' were used by default")
  }
  if(!stat_type %in% c("t.test", "wilcox.test", "LMM")){
    stop("'stat_type' has to be 't.test', 'wilcox.test', or 'LMM'.")
  }

  # Remove columns with p or padj (depending on the p_type)
  if (p_type == "p") {
    stats <- stats[, !grepl("_padj$", colnames(stats))]
  }
  if (p_type == "padj") {
    stats <- stats[, !grepl("_p$", colnames(stats))]
  }

  # Check if data is clustered (has "_within_" pattern)
  is_clustered <- grepl("_within_", colnames(stats))
  clustered_data <- any(is_clustered)

  # Extract appropriate test columns based on stat_type
  if(stat_type == "t.test"){
    test_col <- stats[, grepl("_Ttest_", colnames(stats)), drop = FALSE]
    colnames(test_col) <- sub("_Ttest_.*", "", colnames(test_col))
    colnames(test_col) <- sub("_vs_", "\\/", colnames(test_col))
  } else if(stat_type == "wilcox.test") {
    test_col <- stats[, grepl("_Wilcox_test_", colnames(stats)), drop = FALSE]
    colnames(test_col) <- sub("_Wilcox_test_.*", "", colnames(test_col))
    colnames(test_col) <- sub("_vs_", "\\/", colnames(test_col))
  } else if(stat_type == "LMM"){
    test_col <- stats[, grepl("_LMMtest_", colnames(stats)), drop = FALSE]
    # For LMM, extract everything before _LMMtest_
    colnames(test_col) <- sub("_LMMtest_.*", "", colnames(test_col))
    colnames(test_col) <- sub("_vs_", "\\/", colnames(test_col))
  }

  # Transform the test columns to calculate the -log10(p(test))
  # Replace any p-values <= 0 or >= 1 with valid boundaries before log transformation
  test_col[test_col <= 0 | is.na(test_col)] <- 1e-320  # Smallest positive number R can represent
  test_col[test_col >= 1] <- 1
  test_col <- log10(test_col) * -1
  # Cap at 320 to avoid infinite values
  test_col[test_col > 320 | is.infinite(test_col)] <- 320

  # Identify the fold_change columns
  if(stat_type == "LMM" | clustered_data){
    fc_col <- stats[, grepl("log\\(.*\\/.*\\)", colnames(stats)), drop = FALSE]
  } else {
    fc_col <- stats[, grepl("log\\(.*\\/.*\\)", colnames(stats)), drop = FALSE]
  }

  # Check if all fold changes are logged
  if(sum(grepl("log\\(", colnames(fc_col))) == ncol(fc_col)){
    fc_log = TRUE
  } else if(sum(grepl("log\\(", colnames(fc_col))) == 0){
    fc_log = FALSE
  } else {
    warning("some of the fold-changes in the statistics layer of the 'romics_object' were calculated both prior and after log_transform.")
    warning("Only the log transformed will be used to generate the Volcano plots")
    fc_col <- fc_col[, grepl("log\\(", colnames(fc_col)), drop = FALSE]
  }

  # If not logged then log transform the fold change columns
  if(fc_log == FALSE){
    fc_col = log2(fc_col)
    log_type = 2
    min_fold_change <- log2(min_fold_change)
  } else {
    if(romicsLogCheck(romics_object) & grepl("fun\\|log2", romics_object$steps[grepl("fun\\|log", romics_object$steps)])){
      log_type = 2
    } else {
      log_type = 10
    }
  }

  # Format the colnames so they are identical to the pvalues ones
  if(clustered_data){
    # For clustered: "log(CHR/CTL)_within_Leiden_clust_01" -> "CHR/CTL_within_Leiden_clust_01"
    colnames(fc_col) <- sub("log\\(", "", colnames(fc_col))
    colnames(fc_col) <- sub("\\)_", "_", colnames(fc_col))
    colnames(fc_col) <- sub("_vs_", "\\/", colnames(fc_col))
  } else {
    # For non-clustered: "log(CHR/CTL)" -> "CHR/CTL"
    colnames(fc_col) <- sub("log", "", colnames(fc_col))
    colnames(fc_col) <- sub("\\(", "", colnames(fc_col))
    colnames(fc_col) <- sub("\\)$", "", colnames(fc_col))
  }

  minus_log_p <- log10(p) * -1

  if(sum(colnames(fc_col) %in% colnames(test_col)) != ncol(fc_col)){
    warning("Some of the fold-change columns were not having a equivalent statistical test to generate a Volcano plot.")
  }

  if(plot == "all"){plot <- 1:ncol(fc_col)}
  if(!is.numeric(plot) & sum(!plot %in% 1:ncol(fc_col)) != 0){
    stop(paste0("'plot' as to be either 'all' or a numeric vector with values comprised between 1 and ", ncol(fc_col), "."))
  } else {
    for(i in plot){
      # Ensure we're getting a data frame, not a vector
      fc_col_i <- fc_col[, i, drop = FALSE]
      test_col_match <- colnames(fc_col_i)[1]
      test_col_i <- test_col[, colnames(test_col) == test_col_match, drop = FALSE]

      # Check if match exists
      if(ncol(test_col_i) == 0) {
        warning(paste0("No matching p-value column found for ", test_col_match, ". Skipping this comparison."))
        next
      }

      df <- data.frame(
        ID = rownames(fc_col_i),
        fc = as.numeric(fc_col_i[, 1]),
        p = as.numeric(test_col_i[, 1]),
        stringsAsFactors = FALSE
      )

      # Calculate ylim per comparison if auto_ylim is enabled and ylim not manually specified
      ylim_current <- ylim
      if(is.null(ylim_current) && auto_ylim){
        max_p_value <- max(df$p, na.rm = TRUE)
        ylim_current <- c(0, max_p_value * 1.1)
      }

      # Extract comparison and cluster info from column name
      col_name <- colnames(fc_col_i)[1]

      # Check if this is a clustered result (contains "_within_")
      if(grepl("_within_", col_name)){
        # Extract comparison and cluster information
        comparison <- sub("_within_.*", "", col_name)
        cluster_info <- sub(".*_within_", "", col_name)
        title_suffix <- paste0(" (", cluster_info, ")")
      } else {
        # For standard tests
        comparison <- col_name
        title_suffix <- ""
      }

      class <- rep("non_significant", nrow(df))
      class[df$p > minus_log_p & df$fc < (min_fold_change * -1)] <- "down"
      class[df$p > minus_log_p & df$fc > (min_fold_change)] <- "up"
      class <- paste0(class, "_in_", comparison)
      df$class <- class

      # Check for out-of-bounds points if limits are specified
      if(!is.null(xlim)){
        n_out_x <- sum(df$fc < xlim[1] | df$fc > xlim[2], na.rm=TRUE)
        if(n_out_x > 0){
          warning(paste0("For comparison '", comparison, "': ", n_out_x, " feature(s) have fold-change values outside xlim [", xlim[1], ", ", xlim[2], "]"))
        }
      }
      if(!is.null(ylim_current)){
        n_out_y <- sum(df$p < ylim_current[1] | df$p > ylim_current[2], na.rm=TRUE)
        if(n_out_y > 0){
          warning(paste0("For comparison '", comparison, "': ", n_out_y, " feature(s) have p-value(s) outside ylim [", ylim_current[1], ", ", ylim_current[2], "]"))
        }
      }

      # Prepare labels for top features
      df$label <- ""
      if(label_features && top_features > 0) {
        # Get top upregulated features (passing filters)
        up_features <- df[grepl("up_in_", df$class), ]
        if(nrow(up_features) > 0) {
          if(top_features_by == "abundance") {
            up_features <- up_features[order(up_features$fc, decreasing = TRUE), ]
          } else {
            up_features <- up_features[order(up_features$p, decreasing = TRUE), ]
          }
          n_up_to_label <- min(top_features, nrow(up_features))
          if(n_up_to_label > 0) {
            df$label[df$ID %in% up_features$ID[1:n_up_to_label]] <- df$ID[df$ID %in% up_features$ID[1:n_up_to_label]]
          }
        }
        # Get top downregulated features (passing filters)
        down_features <- df[grepl("down_in_", df$class), ]
        if(nrow(down_features) > 0) {
          if(top_features_by == "abundance") {
            down_features <- down_features[order(down_features$fc, decreasing = FALSE), ]
          } else {
            down_features <- down_features[order(down_features$p, decreasing = TRUE), ]
          }
          n_down_to_label <- min(top_features, nrow(down_features))
          if(n_down_to_label > 0) {
            df$label[df$ID %in% down_features$ID[1:n_down_to_label]] <- df$ID[df$ID %in% down_features$ID[1:n_down_to_label]]
          }
        }
      }

      plot_title <- paste0("Volcano plot for ", comparison, title_suffix)
      yaxis_title <- paste0("-log10(", p_type, "_", stat_type, "_", comparison, ")")
      xaxis_title <- paste0("log", log_type, "(", comparison, ")")

      if(!plotly){
        fig <- ggplot(df, aes(x = fc, y = p, colour = class)) +
          geom_point(alpha = alpha, size = size) +
          theme_ROP() +
          ggtitle(plot_title) +
          xlab(xaxis_title) +
          ylab(yaxis_title) +
          scale_colour_manual(values = colors)

        # Add axis limits if specified
        if(!is.null(xlim) | !is.null(ylim_current)){
          fig <- fig + coord_cartesian(xlim = xlim, ylim = ylim_current, expand = FALSE)
        }

        # Add labels if requested
        if(label_features && any(df$label != "")) {
          fig <- fig + ggrepel::geom_text_repel(aes(label = label),
                                                color = "black",
                                                size = 3,
                                                max.overlaps = Inf,
                                                show.legend = FALSE)
        }
        plot(fig)
      } else {
        # Create text for hover
        hover_text <- paste("ID=", df$ID)
        if(label_features) {
          hover_text[df$label != ""] <- paste("ID=", df$ID[df$label != ""], "(Top feature)")
        }

        # Build xaxis list with optional range
        xaxis_list <- list(title = xaxis_title)
        if(!is.null(xlim)){
          xaxis_list$range <- xlim
        }
        # Build yaxis list with optional range
        yaxis_list <- list(title = yaxis_title)
        if(!is.null(ylim_current)){
          yaxis_list$range <- ylim_current
        }

        fig <- plotly::plot_ly(x = df$fc,
                       y = df$p,
                       color = df$class,
                       colors = colors,
                       type = "scatter", mode = "markers",
                       text = hover_text,
                       hovertemplate = "%{text}<extra></extra>",
                       marker = list(size = size * 4, opacity = alpha)) %>%
          plotly::layout(title = plot_title,
                         xaxis = xaxis_list,
                         yaxis = yaxis_list)

        # Add text annotations for labeled points
        if(label_features && any(df$label != "")) {
          labeled_points <- df[df$label != "", ]
          for(j in 1:nrow(labeled_points)) {
            fig <- fig %>% plotly::add_annotations(
              x = labeled_points$fc[j],
              y = labeled_points$p[j],
              text = labeled_points$label[j],
              showarrow = TRUE,
              arrowhead = 2,
              arrowsize = 0.5,
              ax = 20,
              ay = -20,
              font = list(size = 10)
            )
          }
        }
        print(fig)
      }
    }
  }
}


#' romicsVolcanoByCluster()
#' @description generate volcano plots for all comparisons within a specific cluster
#' @param romics_object a romics object containing statistics
#' @param within character string specifying the cluster to filter (e.g., "Leiden_clust_01"). If NULL, the first cluster matching cluster_factor is used.
#' @param cluster_factor character string specifying the clustering factor to use (e.g., "Leiden_clust", "kmeans_clust"). If NULL, defaults to "Leiden_clust" if available.
#' @param multipanel logical (TRUE or FALSE) to indicate whether to display all plots in a panel (TRUE) or one by one (FALSE)
#' @param stat_type 't.test', 'wilcox.test', 'LMM' or 'auto' to indicate what statistics to use. If 'auto', will detect available tests.
#' @param p_type 'p' or 'padj' to indicate the type of pvalue to consider for the volcano plots
#' @param p numeric value to indicate the maximum pvalue to use for the coloring of the volcano plot
#' @param min_fold_change numeric value indicating the minimum fold change threshold
#' @param colors character vector of length 3 used to color the features (down, non-significant, up)
#' @param plotly logical (TRUE or FALSE) to indicate whether to return interactive plotly plots (TRUE) or static ggplot plots (FALSE)
#' @param label_features logical (TRUE or FALSE) to indicate whether to label top features on the plot
#' @param top_features numeric value indicating the number of top upregulated and downregulated features to label
#' @param top_features_by 'abundance' or 'p' to indicate the criteria for selecting top features
#' @param n_cols numeric value indicating the number of columns in the multipanel layout (NULL for auto square-ish grid)
#' @param return_plot logical (TRUE or FALSE) to indicate whether to return the plot object instead of printing it (default: FALSE)
#' @param xlim numeric vector of length 2 specifying the x-axis limits (e.g., c(-2, 2)). If NULL, limits are determined automatically.
#' @param ylim numeric vector of length 2 specifying the y-axis limits (e.g., c(0, 10)), or a single numeric value to set max y-limit (min is set to 0). If NULL and auto_ylim=TRUE, limits are set based on data max.
#' @param auto_ylim logical (TRUE or FALSE) to automatically set ylim based on the maximum -log10(p) value in the data (default: TRUE). Ignored if ylim is manually specified.
#' @param size numeric value for point size in the plots. Default is 2.
#' @param alpha numeric value for point opacity in the plots (0-1). Default is 0.5.
#' @details generates volcano plots for all comparisons within a specific cluster from any statistical test (t.test, wilcox.test, or LMM). If xlim or ylim are specified, a warning is issued if any data points fall outside these limits.
#' @return Returns plot list invisibly if return_plot=FALSE, or returns the combined plot object if return_plot=TRUE
#' @author Geremy Clair
#' @export
romicsVolcanoByCluster <- function(romics_object, within = NULL, cluster_factor = NULL, multipanel = TRUE, stat_type = "auto",
                                   p_type = "p", p = 0.05, min_fold_change = 0.6,
                                   colors = c("#2cbcb2", "#242021", "#d44e28"),
                                   plotly = FALSE, label_features = FALSE, top_features = 10,
                                   top_features_by = "abundance", n_cols = NULL, return_plot = FALSE,
                                   xlim = NULL, ylim = NULL, auto_ylim = TRUE, size = 2, alpha = 0.5) {
  # Input validation
  if(!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }

  # Determine cluster_factor if not specified
  if(is.null(cluster_factor)) {
    # Default to "Leiden_clust" if available in statistics
    stats_cols <- colnames(romics_object$statistics)
    if(any(grepl("_within_Leiden_clust_", stats_cols))) {
      cluster_factor <- "Leiden_clust"
    } else {
      # Try to find the first available cluster factor
      cluster_patterns <- unique(gsub(".*_within_([^_]+)_[0-9]+.*", "\\1",
                                      stats_cols[grepl("_within_", stats_cols)]))
      if(length(cluster_patterns) > 0) {
        cluster_factor <- cluster_patterns[1]
      } else {
        stop("No clustered statistics found. Please run clustering analysis first.")
      }
    }
  }
  if(!is.character(cluster_factor)) {
    stop("'cluster_factor' should be a character string")
  }

  # If within is not specified, find the first cluster matching cluster_factor
  if(is.null(within)) {
    stats_cols <- colnames(romics_object$statistics)
    matching_clusters <- unique(gsub(paste0(".*_within_(", cluster_factor, "_[0-9]+).*"), "\\1",
                                      stats_cols[grepl(paste0("_within_", cluster_factor, "_"), stats_cols)]))
    if(length(matching_clusters) == 0) {
      stop(paste0("No clusters found for cluster_factor '", cluster_factor,
                  "'. Available factors may be: ", paste(unique(gsub(".*_within_([^_]+)_.*", "\\1",
                                      stats_cols[grepl("_within_", stats_cols)])), collapse = ", ")))
    }
    within <- matching_clusters[1]
    message(paste0("Using cluster: ", within))
  }
  if(!is.character(within)) {
    stop("'within' should be a character string")
  }
  if(missing(multipanel)){multipanel = TRUE}
  if(!is.logical(multipanel)){
    warning("'multipanel' was not logical (TRUE/FALSE). TRUE was used by default.")
    multipanel = TRUE
  }
  if(missing(stat_type)){stat_type = "auto"}
  if(!stat_type %in% c("t.test", "wilcox.test", "LMM", "auto")){
    stop("'stat_type' has to be 't.test', 'wilcox.test', 'LMM', or 'auto'.")
  }
  if(missing(p_type)){p_type = "p"}
  if(!p_type %in% c("p", "padj")){
    stop("'p_type' has to be either 'p' or 'padj'")
  }
  if(missing(colors)){
    colors = c("#2cbcb2", "#242021", "#d44e28")
  }
  if(!is.character(colors) | length(colors) != 3){
    warning("'colors' should be a color character vector of length 3. The defaults colors were used")
    colors = c("#2cbcb2", "#242021", "#d44e28")
  }
  if(missing(p)){p = 0.05}
  if(!is.numeric(p) | p > 1 | p < 0){
    stop("'p' should be numeric and comprised between 0 and 1.")
  }
  if(missing(min_fold_change)){min_fold_change = 0.6}
  if(!is.numeric(min_fold_change) | min_fold_change < 0){
    stop("'min_fold_change' should be numeric and comprised higher than 0.")
  }
  if(missing(plotly)){plotly = FALSE}
  if(!is.logical(plotly)){
    warning("'plotly' was not logical (TRUE/FALSE). FALSE was used by default (ggplot output).")
    plotly = FALSE
  }
  if(missing(label_features)){label_features = FALSE}
  if(!is.logical(label_features)){
    warning("'label_features' was not logical (TRUE/FALSE). FALSE was used by default.")
    label_features = FALSE
  }
  if(missing(top_features)){top_features = 10}
  if(!is.numeric(top_features) | top_features < 0 | top_features != round(top_features)){
    warning("'top_features' should be a positive integer. 10 was used by default.")
    top_features = 10
  }
  if(missing(top_features_by)){top_features_by = "abundance"}
  if(!top_features_by %in% c("abundance", "p")){
    warning("'top_features_by' should be either 'abundance' or 'p'. 'abundance' was used by default.")
    top_features_by = "abundance"
  }
  if(!is.null(n_cols)){
    if(!is.numeric(n_cols) | n_cols < 1 | n_cols != round(n_cols)){
      warning("'n_cols' should be a positive integer or NULL (auto). NULL was used by default.")
      n_cols = NULL
    }
  }
  if(missing(return_plot)){return_plot = FALSE}
  if(!is.logical(return_plot)){
    warning("'return_plot' was not logical (TRUE/FALSE). FALSE was used by default.")
    return_plot = FALSE
  }
  if(missing(xlim)){xlim = NULL}
  if(!is.null(xlim)){
    if(!is.numeric(xlim) | length(xlim) != 2){
      warning("'xlim' should be a numeric vector of length 2 (e.g., c(-2, 2)). NULL was used by default.")
      xlim = NULL
    } else if(xlim[1] >= xlim[2]){
      warning("'xlim' should have xlim[1] < xlim[2]. NULL was used by default.")
      xlim = NULL
    }
  }
  if(missing(ylim)){ylim = NULL}
  if(!is.null(ylim)){
    if(is.numeric(ylim) & length(ylim) == 1){
      # Allow single numeric value as max y-limit
      ylim <- c(0, ylim)
    } else if(!is.numeric(ylim) | length(ylim) != 2){
      warning("'ylim' should be a numeric vector of length 2 (e.g., c(0, 10)) or a single numeric value for max. NULL was used by default.")
      ylim = NULL
    } else if(ylim[1] >= ylim[2]){
      warning("'ylim' should have ylim[1] < ylim[2]. NULL was used by default.")
      ylim = NULL
    }
  }
  if(missing(auto_ylim)){auto_ylim = TRUE}
  if(!is.logical(auto_ylim)){
    warning("'auto_ylim' was not logical (TRUE/FALSE). TRUE was used by default.")
    auto_ylim = TRUE
  }

  # Extract stats
  stats <- romics_object$statistics
  stats <- stats[, unique(colnames(stats))]

  # Auto-detect stat_type if not specified
  if(stat_type == "auto"){
    has_lmm_clustered    <- sum(grepl(paste0("_within_", within, "_LMMtest_"),    colnames(stats))) > 0
    has_ttest_clustered  <- sum(grepl(paste0("_within_", within, "_Ttest_"),      colnames(stats))) > 0
    has_wilcox_clustered <- sum(grepl(paste0("_within_", within, "_Wilcox_test_"), colnames(stats))) > 0

    if(has_lmm_clustered) {
      stat_type = "LMM"
      print(paste0("'stat_type' was set to 'LMM' (auto-detected for cluster '", within, "')"))
    } else if(has_ttest_clustered) {
      stat_type = "t.test"
      print(paste0("'stat_type' was set to 't.test' (auto-detected for cluster '", within, "')"))
    } else if(has_wilcox_clustered) {
      stat_type = "wilcox.test"
      print(paste0("'stat_type' was set to 'wilcox.test' (auto-detected for cluster '", within, "')"))
    } else {
      stop(paste0("No statistical tests found for cluster '", within, "'. Please check the cluster name or run appropriate tests first."))
    }
  }

  # Remove columns with p or padj (depending on the p_type)
  if (p_type == "p") {
    stats <- stats[, !grepl("_padj$", colnames(stats))]
  }
  if (p_type == "padj") {
    stats <- stats[, !grepl("_p$", colnames(stats))]
  }

  # Filter columns for the specified cluster based on stat_type
  # The 'within' parameter specifies the actual cluster name to search for
  if(stat_type == "t.test"){
    fc_pattern   <- paste0("log\\(.*\\/.*\\)_within_", within, "$|log\\(.*\\/.*\\)_within_", within, "_")
    test_pattern <- paste0("_within_", within, "_Ttest_")
  } else if(stat_type == "wilcox.test"){
    fc_pattern   <- paste0("log\\(.*\\/.*\\)_within_", within, "$|log\\(.*\\/.*\\)_within_", within, "_")
    test_pattern <- paste0("_within_", within, "_Wilcox_test_")
  } else if(stat_type == "LMM"){
    fc_pattern   <- paste0("log\\(.*\\/.*\\)_within_", within, "$|log\\(.*\\/.*\\)_within_", within, "_")
    test_pattern <- paste0("_within_", within, "_LMMtest_")
  }

  fc_cols   <- grepl(fc_pattern,   colnames(stats))
  test_cols <- grepl(test_pattern, colnames(stats))

  if(sum(fc_cols) == 0 | sum(test_cols) == 0) {
    stop(paste0("No statistics found for cluster '", within, "' with stat_type '", stat_type,
                "'. Please check the cluster name and test type."))
  }

  # Extract the filtered columns and convert to matrix to avoid list issues
  fc_col   <- as.matrix(stats[, fc_cols,   drop = FALSE])
  test_col <- as.matrix(stats[, test_cols, drop = FALSE])

  # Process test columns based on stat_type — strip _within_<cluster>_<test>_... suffix
  if(stat_type == "t.test"){
    colnames(test_col) <- sub(paste0("_within_", within, "_Ttest_.*"),      "", colnames(test_col))
  } else if(stat_type == "wilcox.test"){
    colnames(test_col) <- sub(paste0("_within_", within, "_Wilcox_test_.*"), "", colnames(test_col))
  } else if(stat_type == "LMM"){
    colnames(test_col) <- sub(paste0("_within_", within, "_LMMtest_.*"),     "", colnames(test_col))
  }
  colnames(test_col) <- sub("_vs_", "\\/", colnames(test_col))

  # Replace any p-values <= 0 or >= 1 with valid boundaries before log transformation
  test_col[test_col <= 0 | is.na(test_col)] <- 1e-320  # Smallest positive number R can represent
  test_col[test_col >= 1] <- 1
  test_col <- log10(test_col) * -1
  # Cap at 320 to avoid infinite values
  test_col[test_col > 320 | is.infinite(test_col)] <- 320

  # Check if fold changes are logged
  if(sum(grepl("log\\(", colnames(fc_col))) == ncol(fc_col)){
    fc_log = TRUE
  } else if(sum(grepl("log\\(", colnames(fc_col))) == 0){
    fc_log = FALSE
  } else {
    warning("some of the fold-changes were calculated both prior and after log_transform. Only the log transformed will be used.")
    fc_col <- fc_col[, grepl("log\\(", colnames(fc_col)), drop = FALSE]
    fc_log = TRUE
  }

  # Log transform if needed
  if(fc_log == FALSE){
    fc_col = log2(fc_col)
    log_type = 2
    min_fold_change_adj <- log2(min_fold_change)
  } else {
    if(romicsLogCheck(romics_object) & grepl("fun\\|log2", romics_object$steps[grepl("fun\\|log", romics_object$steps)])){
      log_type = 2
    } else {
      log_type = 10
    }
    min_fold_change_adj <- min_fold_change
  }

  # Format fold-change column names (convert back to data frame for consistency with rest of code)
  fc_col <- as.data.frame(fc_col, stringsAsFactors = FALSE)
  test_col <- as.data.frame(test_col, stringsAsFactors = FALSE)

  colnames(fc_col) <- sub("log\\(", "", colnames(fc_col))
  colnames(fc_col) <- sub("\\)_within_.*", "", colnames(fc_col))
  colnames(fc_col) <- sub("_vs_", "\\/", colnames(fc_col))

  minus_log_p <- log10(p) * -1

  # Create list to store plots
  plot_list <- list()

  # Generate plots for each comparison
  for(i in 1:ncol(fc_col)){
    fc_col_i       <- fc_col[, i, drop = FALSE]
    test_col_match <- colnames(fc_col_i)[1]
    test_col_i     <- test_col[, colnames(test_col) == test_col_match, drop = FALSE]

    if(ncol(test_col_i) == 0) {
      warning(paste0("No matching p-value column found for ", test_col_match, ". Skipping this comparison."))
      next
    }

    df <- data.frame(
      ID = rownames(fc_col_i),
      fc = as.numeric(fc_col_i[, 1]),
      p  = as.numeric(test_col_i[, 1]),
      stringsAsFactors = FALSE
    )

    # Calculate ylim per comparison if auto_ylim is enabled and ylim not manually specified
    ylim_current <- ylim
    if(is.null(ylim_current) && auto_ylim){
      max_p_value <- max(df$p, na.rm = TRUE)
      ylim_current <- c(0, max_p_value * 1.1)
    }

    comparison <- colnames(fc_col_i)[1]

    class <- rep("non_significant", nrow(df))
    class[df$p > minus_log_p & df$fc < (min_fold_change_adj * -1)] <- "down"
    class[df$p > minus_log_p & df$fc >  min_fold_change_adj]        <- "up"
    class <- paste0(class, "_in_", comparison)
    df$class <- class

    # Check for out-of-bounds points if limits are specified
    if(!is.null(xlim)){
      n_out_x <- sum(df$fc < xlim[1] | df$fc > xlim[2], na.rm=TRUE)
      if(n_out_x > 0){
        warning(paste0("For comparison '", comparison, "' in cluster '", within, "': ", n_out_x, " feature(s) have fold-change values outside xlim [", xlim[1], ", ", xlim[2], "]"))
      }
    }
    if(!is.null(ylim_current)){
      n_out_y <- sum(df$p < ylim_current[1] | df$p > ylim_current[2], na.rm=TRUE)
      if(n_out_y > 0){
        warning(paste0("For comparison '", comparison, "' in cluster '", within, "': ", n_out_y, " feature(s) have p-value(s) outside ylim [", ylim_current[1], ", ", ylim_current[2], "]"))
      }
    }

    # Prepare labels for top features
    df$label <- ""
    if(label_features && top_features > 0) {
      up_features <- df[grepl("up_in_", df$class), ]
      if(nrow(up_features) > 0) {
        if(top_features_by == "abundance") {
          up_features <- up_features[order(up_features$fc, decreasing = TRUE), ]
        } else {
          up_features <- up_features[order(up_features$p, decreasing = TRUE), ]
        }
        n_up_to_label <- min(top_features, nrow(up_features))
        if(n_up_to_label > 0) {
          df$label[df$ID %in% up_features$ID[1:n_up_to_label]] <- df$ID[df$ID %in% up_features$ID[1:n_up_to_label]]
        }
      }
      down_features <- df[grepl("down_in_", df$class), ]
      if(nrow(down_features) > 0) {
        if(top_features_by == "abundance") {
          down_features <- down_features[order(down_features$fc, decreasing = FALSE), ]
        } else {
          down_features <- down_features[order(down_features$p, decreasing = TRUE), ]
        }
        n_down_to_label <- min(top_features, nrow(down_features))
        if(n_down_to_label > 0) {
          df$label[df$ID %in% down_features$ID[1:n_down_to_label]] <- df$ID[df$ID %in% down_features$ID[1:n_down_to_label]]
        }
      }
    }

    plot_title  <- paste0("Volcano plot for ", comparison)
    yaxis_title <- paste0("-log10(", p_type, "_", stat_type, ")")
    xaxis_title <- paste0("log", log_type, "(", comparison, ")")

    if(!plotly){
      fig <- ggplot(df, aes(x = fc, y = p, colour = class)) +
        geom_point(alpha = alpha, size = size) +
        theme_ROP() +
        ggtitle(plot_title) +
        xlab(xaxis_title) +
        ylab(yaxis_title) +
        scale_colour_manual(values = colors) +
        theme(legend.position = "bottom",
              plot.title  = element_text(hjust = 0.5, size = 10),
              axis.title  = element_text(size = 9))

      # Add axis limits if specified
      if(!is.null(xlim) | !is.null(ylim_current)){
        fig <- fig + coord_cartesian(xlim = xlim, ylim = ylim_current, expand = FALSE)
      }

      if(label_features && any(df$label != "")) {
        fig <- fig + ggrepel::geom_text_repel(aes(label = label),
                                              color        = "black",
                                              size         = 2.5,
                                              max.overlaps = Inf,
                                              show.legend  = FALSE)
      }

      plot_list[[i]] <- fig

    } else {
      hover_text <- paste("ID=", df$ID)
      if(label_features) {
        hover_text[df$label != ""] <- paste("ID=", df$ID[df$label != ""], "(Top feature)")
      }

      # Build xaxis list with optional range
      xaxis_list <- list(title = xaxis_title)
      if(!is.null(xlim)){
        xaxis_list$range <- xlim
      }
      # Build yaxis list with optional range
      yaxis_list <- list(title = yaxis_title)
      if(!is.null(ylim_current)){
        yaxis_list$range <- ylim_current
      }

      fig <- plot_ly(x = df$fc,
                     y = df$p,
                     color = df$class,
                     colors = colors,
                     type = "scatter", mode = "markers",
                     text = hover_text,
                     hovertemplate = "%{text}<extra></extra>",
                     marker = list(size = size * 4, opacity = alpha)) %>%
        plotly::layout(title  = plot_title,
                       xaxis  = xaxis_list,
                       yaxis  = yaxis_list)

      if(label_features && any(df$label != "")) {
        labeled_points <- df[df$label != "", ]
        for(j in 1:nrow(labeled_points)) {
          fig <- fig %>% plotly::add_annotations(
            x         = labeled_points$fc[j],
            y         = labeled_points$p[j],
            text      = labeled_points$label[j],
            showarrow = TRUE,
            arrowhead = 2,
            arrowsize = 0.5,
            ax        = 20,
            ay        = -20,
            font      = list(size = 10)
          )
        }
      }

      plot_list[[i]] <- fig
    }
  }

  # Remove NULL entries
  plot_list <- plot_list[!sapply(plot_list, is.null)]

  if(length(plot_list) == 0) {
    stop("No plots were generated. Please check your input parameters.")
  }

  # Overall title string
  overall_title_text <- paste0("Volcano plots within ", within)

  # Create combined plot
  if(!plotly){
    if(multipanel){
      n_plots <- length(plot_list)
      if(is.null(n_cols)){
        n_cols_used <- ceiling(sqrt(n_plots))
      } else {
        n_cols_used <- min(n_cols, n_plots)
      }
      n_rows <- ceiling(n_plots / n_cols_used)

      grid <- cowplot::plot_grid(
        plotlist = plot_list,
        ncol     = n_cols_used,
        nrow     = n_rows,
        align    = "hv"
      )

      title_grob <- cowplot::ggdraw() +
        cowplot::draw_label(
          overall_title_text,
          fontface = "bold",
          x        = 0.5,
          hjust    = 0.5,
          size     = 13
        )

      combined_plot <- cowplot::plot_grid(
        title_grob, grid,
        ncol        = 1,
        rel_heights = c(0.05, 1)
      )

      if(return_plot){
        return(combined_plot)
      } else {
        print(combined_plot)
        invisible(plot_list)
      }

    } else {
      # One by one display
      if(return_plot){
        # If returning plots one by one, just return the list
        return(plot_list)
      } else {
        for(i in seq_along(plot_list)) {
          cat(paste0("\n--- Plot ", i, " of ", length(plot_list), " ---\n"))
          print(plot_list[[i]])
          if(i < length(plot_list)) {
            invisible(readline(prompt = "Press [enter] to see next plot"))
          }
        }
        invisible(plot_list)
      }
    }

  } else {
    # For plotly
    if(multipanel){
      n_plots <- length(plot_list)
      if(is.null(n_cols)){
        n_cols_used <- ceiling(sqrt(n_plots))
      } else {
        n_cols_used <- min(n_cols, n_plots)
      }
      n_rows <- ceiling(n_plots / n_cols_used)

      combined_plot <- plotly::subplot(plot_list,
                                       nrows  = n_rows,
                                       ncols  = n_cols_used,
                                       shareX = FALSE,
                                       shareY = FALSE,
                                       margin = 0.08) %>%
        plotly::layout(
          title = list(
            text    = overall_title_text,
            x       = 0.5,
            xanchor = "center"
          ),
          showlegend = TRUE,
          height     = 300 * n_rows,
          width      = 400 * n_cols_used
        )

      if(return_plot){
        return(combined_plot)
      } else {
        print(combined_plot)
        invisible(plot_list)
      }

    } else {
      if(return_plot){
        return(plot_list)
      } else {
        for(i in seq_along(plot_list)) {
          cat(paste0("\n--- Plot ", i, " of ", length(plot_list), " ---\n"))
          print(plot_list[[i]])
          if(i < length(plot_list)) {
            invisible(readline(prompt = "Press [enter] to see next plot"))
          }
        }
        invisible(plot_list)
      }
    }
  }
}

#' romicsPathwayVolcano()
#' @description Generate pathway-colored volcano plots where features in user-defined lists are colored by pathway
#' @param romics_object A romics_object containing statistical test results
#' @param pathway_list A list of character vectors, each containing feature IDs. Names of list elements become pathway labels
#' @param comparison Character string specifying which comparison to plot (must match a column in statistics layer)
#' @param p_type Character. P-value type: "p" or "padj". Default: "p"
#' @param p Numeric. P-value threshold for significance line. Default: 0.05
#' @param min_fold_change Numeric. Fold-change threshold line. Default: 0
#' @param pathway_colors Character vector of colors for pathways. If NULL or shorter than pathway_list length, colors are auto-generated. Default: NULL
#' @param pathway_shapes Numeric vector of shapes for pathways (1=circle, 4=diamond, etc.). Default: 1 (circles)
#' @param background_color Character. Color for features not in any pathway. Default: "gray60"
#' @param background_shape Numeric. Shape for background features (1=circle). Default: 1
#' @param xlim Numeric vector of length 2 for x-axis limits. Default: NULL (auto)
#' @param ylim Numeric vector of length 2 for y-axis limits. Default: NULL (auto)
#' @param size Numeric. Point size. Default: 2
#' @param alpha Numeric. Point transparency for significant features (0-1). Default: 0.8
#' @param alpha_nonsignificant Numeric. Point transparency for non-significant features (0-1). Default: 0.4
#' @param label_features Logical. Whether to label pathway features with colors matching their pathway. Default: FALSE
#' @param label_pathway_features Logical. Alias for label_features (both enable pathway feature labeling). Default: FALSE
#' @param plotly Logical. Whether to return an interactive plotly plot (TRUE) or static ggplot (FALSE). Default: FALSE
#' @details Creates a volcano plot where features are colored by membership in user-defined pathways. All pathways are plotted on the same axes with p-value and fold-change threshold lines shown as gray dotted lines.
#' @return A ggplot2 object (if plotly=FALSE) or a plotly object (if plotly=TRUE)
#' @author Geremy Clair
#' @export
romicsPathwayVolcano <- function(romics_object,
                                pathway_list,
                                comparison,
                                p_type = "p",
                                p = 0.05,
                                min_fold_change = 0,
                                pathway_colors = NULL,
                                pathway_shapes = NULL,
                                background_color = "gray60",
                                background_shape = 1,
                                xlim = NULL,
                                ylim = NULL,
                                size = 2,
                                alpha = 0.8,
                                alpha_nonsignificant = 0.4,
                                label_features = FALSE,
                                label_pathway_features = FALSE,
                                plotly = FALSE) {

  if (!is.romicsObject(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }
  if (!is.list(pathway_list) || length(pathway_list) == 0) {
    stop("pathway_list must be a non-empty list of character vectors")
  }
  if (is.null(names(pathway_list))) {
    stop("pathway_list elements must be named (pathway names)")
  }
  if (!is.character(comparison) || length(comparison) != 1) {
    stop("comparison must be a single character string")
  }
  if (missing(plotly)) {
    plotly <- FALSE
  }
  if (!is.logical(plotly)) {
    warning("'plotly' should be logical (TRUE/FALSE). FALSE was used by default.")
    plotly <- FALSE
  }

  stats <- romics_object$statistics

  # Find p-value column matching comparison
  if (p_type == "padj") {
    p_col_pattern <- paste0("^", gsub("([.*+?^${}()|\\[\\]\\\\])", "\\\\\\1", comparison), "_.*_padj$")
  } else {
    p_col_pattern <- paste0("^", gsub("([.*+?^${}()|\\[\\]\\\\])", "\\\\\\1", comparison), "_.*_p$")
  }
  p_cols <- grep(p_col_pattern, colnames(stats), value = TRUE)

  if (length(p_cols) == 0) {
    stop(paste("No p-value columns found for comparison:", comparison))
  }

  p_col <- p_cols[1]

  # Find fold-change column
  # Convert _vs_ to / for matching fold-change column names
  comparison_fc <- gsub("_vs_", "/", comparison)
  # Escape special characters
  comparison_fc_escaped <- gsub("([.*+?^${}()|\\[\\]\\\\])", "\\\\\\1", comparison_fc)
  fc_col_pattern <- paste0("^log\\(", comparison_fc_escaped, "\\)$|^\\(", comparison_fc_escaped, "\\)$")
  fc_cols <- grep(fc_col_pattern, colnames(stats), value = TRUE)

  if (length(fc_cols) == 0) {
    stop(paste("No fold-change columns found for comparison:", comparison))
  }

  fc_col <- fc_cols[1]

  # Extract data
  pvals <- suppressWarnings(as.numeric(stats[[p_col]]))
  fc <- suppressWarnings(as.numeric(stats[[fc_col]]))

  # Log transform fold-change if needed
  if (!grepl("^log", fc_col)) {
    fc <- log2(fc)
    min_fold_change <- log2(min_fold_change)
  }

  # Store p threshold as p_threshold to avoid naming conflict with column
  p_threshold <- p

  # Handle extreme p-values before log transformation
  pvals[pvals <= 0 | is.na(pvals)] <- 1e-320
  pvals[pvals >= 1] <- 1
  log_pvals <- -log10(pvals)
  log_pvals[is.infinite(log_pvals)] <- 320

  # Create data frame
  df <- data.frame(
    ID = rownames(stats),
    fc = fc,
    p = log_pvals,
    pathway = "Background",
    stringsAsFactors = FALSE
  )

  # Assign pathways
  for (pathway_name in names(pathway_list)) {
    features_in_pathway <- pathway_list[[pathway_name]]
    df$pathway[df$ID %in% features_in_pathway] <- pathway_name
  }

  # Calculate significance (above p-value threshold AND above fold-change threshold)
  df$significant <- (df$p > -log10(p_threshold)) & (abs(df$fc) > min_fold_change)

  # Set up colors
  if (is.null(pathway_colors)) {
    pathway_colors <- RColorBrewer::brewer.pal(n = min(max(3, length(pathway_list)), 12), name = "Set1")[1:length(pathway_list)]
  } else if (length(pathway_colors) < length(pathway_list)) {
    warning("pathway_colors shorter than pathway_list. Using provided colors and cycling.")
    pathway_colors <- rep(pathway_colors, length.out = length(pathway_list))
  }

  # Set up shapes (convert to fillable shapes: 21=circle, 23=diamond, 22=square, 24=triangle)
  if (is.null(pathway_shapes)) {
    pathway_shapes <- rep(23, length(pathway_list))  # Diamond (23) by default for pathways
  } else if (length(pathway_shapes) < length(pathway_list)) {
    pathway_shapes <- rep(pathway_shapes, length.out = length(pathway_list))
  }
  # Convert to fillable shapes if needed (1->21 circle, 4->23 diamond, etc.)
  pathway_shapes[pathway_shapes == 1] <- 21  # circle outline -> filled circle
  pathway_shapes[pathway_shapes == 4] <- 23  # diamond outline -> filled diamond
  background_shape_filled <- if (background_shape == 1) 21 else if (background_shape == 4) 23 else background_shape

  # Create color and shape mappings
  color_map <- c(Background = background_color, setNames(pathway_colors, names(pathway_list)))
  shape_map <- c(Background = background_shape_filled, setNames(pathway_shapes, names(pathway_list)))

  # Create plot with fill color and alpha based on significance
  p <- ggplot2::ggplot(df, ggplot2::aes(x = fc, y = p, fill = pathway, shape = pathway, alpha = significant, color = pathway, text = ID)) +
    ggplot2::geom_point(size = size, stroke = 0.1) +
    ggplot2::scale_color_manual(values = color_map, guide = "none") +
    ggplot2::geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "gray50", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = c(-min_fold_change, min_fold_change), linetype = "dashed", color = "gray50", linewidth = 0.5) +
    ggplot2::scale_fill_manual(values = color_map) +
    ggplot2::scale_shape_manual(values = shape_map) +
    ggplot2::scale_alpha_manual(values = c("TRUE" = alpha, "FALSE" = alpha_nonsignificant), guide = "none") +
    ggplot2::labs(
      title = paste("Pathway Volcano Plot:", comparison),
      x = "log2(Fold Change)",
      y = paste0("-log10(", p_type, ")"),
      fill = "Pathway",
      shape = "Pathway"
    ) +
    theme_ROP() +
    ggplot2::theme(
      legend.position = "right",
      legend.box = "vertical"
    )

  # Apply axis limits if specified
  if (!is.null(xlim)) {
    p <- p + ggplot2::xlim(xlim)
  }
  if (!is.null(ylim)) {
    p <- p + ggplot2::ylim(ylim)
  }

  # Add feature labels if requested
  if (label_features || label_pathway_features) {
    # Label pathway features colored by pathway
    pathway_features <- df[df$pathway != "Background", ]
    if (nrow(pathway_features) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = pathway_features,
        ggplot2::aes(label = ID, color = pathway),
        size = 3,
        show.legend = FALSE
      )
    }
  }

  # Convert to plotly if requested
  if (plotly) {
    p <- plotly::ggplotly(p, tooltip = c("text", "x", "y"))
  }

  return(p)
}
