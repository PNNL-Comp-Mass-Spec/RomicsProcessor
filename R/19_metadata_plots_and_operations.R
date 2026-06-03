#' romicsMetadataRelations()
#' @description display relationships between different factors contained in the metadata layers of an Romics_object
#' @param romics_object A romics_object created using romicsCreateObject().
#' @param factors A character vector (length>2) containing the list of factors of the romics_object. the list of factors can be identified by using the function romicsFactorNames()
#' @param textsize Size of the text labels (default: 3)
#' @param colors Optional vector of colors for the strata. If NULL, colors are generated automatically.
#' @details generates a plot displaying relationships between different factors
#' @return a ggplot
#' @author Geremy Clair
#' @export
#'
romicsMetadataRelations <- function(romics_object = romics_object,
                                    factors = c("factor1", "factor2", "factor3"),
                                    textsize = 3,
                                    colors = NULL) {
  # Check for required packages and attach ggalluvial
  if (!requireNamespace("ggalluvial", quietly = TRUE)) {
    stop("The ggalluvial package is required for this function. Please install it with: install.packages('ggalluvial')")
  }
  # Attach ggalluvial so stat_stratum is available
  library(ggalluvial)

  # Input validation
  if(!is.romicsObject(romics_object) || missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }

  if(missing(factors)) {
    stop("factors parameter is required")
  }

  if(!all(factors %in% romicsFactorNames(romics_object))) {
    stop("Some factors specified are not present in the romics_object. Use romicsFactorNames() to see available factors.")
  }

  if(length(factors) < 2) {
    stop("At least 2 factors must be specified to show relationships")
  }

  if(!is.numeric(textsize) || length(textsize) != 1) {
    stop("textsize must be a single numeric value")
  }

  # Extract metadata for the selected factors
  meta <- as.data.frame(t(romics_object$metadata))

  # Create a data frame with only the requested factors
  m <- data.frame(matrix(ncol = length(factors), nrow = nrow(meta)))
  for(i in 1:length(factors)) {
    m[,i] <- meta[, factors[i]]
    colnames(m)[i] <- factors[i]
  }

  # Convert to lodes format
  ml <- ggalluvial::to_lodes_form(m)

  # Generate colors if not provided
  if(is.null(colors)) {
    unique_strata <- unique(ml$stratum)
    n_colors_needed <- length(unique_strata)

    # Try to use ROP_colors first (40 colors available)
    if(exists("ROP_colors") && n_colors_needed <= length(ROP_colors)) {
      colors <- ROP_colors[1:n_colors_needed]
    } else if(requireNamespace("viridis", quietly = TRUE)) {
      # Fall back to viridis if ROP_colors doesn't have enough colors
      colors <- viridis::viridis(n_colors_needed)
    } else {
      # Final fallback to a basic color palette
      default_colors <- c(
        "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
        "#D55E00", "#CC79A7", "#999999", "#000000", "#E69F00",
        "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00"
      )
      # Repeat the colors if needed
      colors <- rep(default_colors, length.out = n_colors_needed)
    }
  }

  # Create the plot
  ggplot(ml, aes(x = x, stratum = stratum,
                 alluvium = alluvium, fill = stratum, label = stratum)) +
    ggalluvial::geom_flow(stat = "alluvium", lode.guidance = "frontback",
                          color = "NA") +
    ggalluvial::geom_stratum(alpha = 0.5) +
    geom_text(stat = "stratum", size = textsize) +
    theme_ROP() +
    xlab("Metadata") +
    ylab("Frequency") +
    scale_fill_manual(values = colors)
}


#' factorCountLevels()
#' @description Indicates in a vector how many times each level of a factor is detected, with optional visualization
#' @param romics_object A romics_object created using romicsCreateObject().
#' @param factor A character vector (length=1) containing a factor of the romics_object. The list of factors can be identified by using the function romicsFactorNames().
#' @param plot Logical indicating whether to return a ggplot bar chart. Default is FALSE (returns numeric vector).
#' @param textsize Numeric value for text size in the plot. Default is 3. Only used when plot=TRUE.
#' @param fill_color Color for the bars in the plot. Default is "gray70". Only used when plot=TRUE.
#' @param order Logical indicating whether to order bars by count (descending). Default is FALSE. Only used when plot=TRUE.
#' @details Generates a count of how many times each level is seen in a factor. Can return either a named numeric vector or a ggplot bar chart.
#' @return If plot=FALSE: A named numeric vector with counts for each factor level. If plot=TRUE: A ggplot object showing the counts as a bar chart.
#' @examples \dontrun{
#' # Get counts as a named vector
#' counts <- factorCountLevels(romics_object = my_romics, factor = "Treatment")
#'
#' # Generate a plot
#' factorCountLevels(romics_object = my_romics, factor = "Treatment", plot = TRUE)
#'
#' # Customize the plot
#' factorCountLevels(romics_object = my_romics, factor = "Sex",
#'                   plot = TRUE, textsize = 4, fill_color = "gray50", order = TRUE)
#'
#' # Access specific level count
#' counts <- factorCountLevels(romics_object = my_romics, factor = "Sex")
#' counts["Male"]
#' }
#' @author Geremy Clair
#' @export
factorCountLevels <- function(romics_object = romics_object,
                              factor = "factor",
                              plot = FALSE,
                              textsize = 3,
                              fill_color = "gray70",
                              order = FALSE) {

  # Check for ggplot2 if plotting is requested
  if(plot && !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required for plotting. Please install it with: install.packages('ggplot2')")
  }

  # Input validation
  if(!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }

  if(missing(factor)) {
    factor <- romics_object$main_factor
  }

  if(!all(factor %in% romicsFactorNames(romics_object)) | length(factor) != 1) {
    stop("'factor' has to be a character vector of length=1 contained in the romics_object. Factors can be identified by using the function romicsFactorNames()")
  }

  if(!is.logical(plot) | length(plot) != 1) {
    stop("'plot' must be a single logical value (TRUE or FALSE)")
  }

  if(!is.numeric(textsize) | length(textsize) != 1) {
    stop("'textsize' has to be a single numerical value")
  }

  if(!is.logical(order) | length(order) != 1) {
    stop("'order' must be a single logical value (TRUE or FALSE)")
  }

  # Extract metadata
  meta <- as.data.frame(t(romics_object$metadata))
  meta <- meta[, colnames(meta) == factor]

  # Count factor levels
  f <- as.factor(meta)
  counts <- table(f)

  # Convert to named numeric vector
  result <- as.numeric(counts)
  names(result) <- names(counts)

  # Return plot if requested
  if(plot) {
    # Create data frame for plotting
    plot_data <- data.frame(
      Level = names(result),
      Count = result,
      stringsAsFactors = FALSE
    )

    # Order by count if requested
    if(order) {
      plot_data$Level <- factor(plot_data$Level,
                                levels = plot_data$Level[order(plot_data$Count, decreasing = TRUE)])
    }

    # Create the plot with no outline
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Level, y = Count)) +
      ggplot2::geom_bar(stat = "identity", fill = fill_color, color = NA) +
      ggplot2::geom_text(ggplot2::aes(label = Count),
                         vjust = -0.5,
                         size = textsize) +
      ggplot2::labs(
        title = paste("Count of Levels for Factor:", factor),
        x = factor,
        y = "Count"
      ) +
      ggplot2::ylim(0, max(plot_data$Count) * 1.1)  # Add 10% space for labels

    # Apply theme_ROP if available, otherwise use theme_minimal
    if(exists("theme_ROP") && is.function(theme_ROP)) {
      p <- p + theme_ROP()
    } else {
      p <- p + ggplot2::theme_minimal() +
        ggplot2::theme(
          panel.grid.major.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )
    }

    return(p)
  }

  # Return numeric vector if plot=FALSE
  return(result)
}


#' @title Create Stacked Barplot of Factor Distribution
#' @description Creates and displays stacked barplots showing the percentage distribution of one factor's levels across different levels of another factor from a romics_object.
#'
#' @param romics_object A romics_object created using romicsCreateObject().
#' @param split_by_factor A character vector (length=1) containing a factor of the romics_object that will define the x-axis categories. The list of factors can be identified by using the function romicsFactorNames()
#' @param frequency_for_factor A character vector (length=1) containing a factor of the romics_object that will be used to calculate frequency distributions within each level of split_by_factor. The list of factors can be identified by using the function romicsFactorNames()
#' @param title A character string specifying the plot title. If NULL (default), an automatic title will be generated.
#' @param xlab A character string specifying the x-axis label. If NULL (default), the split_by_factor name will be used.
#' @param ylab A character string specifying the y-axis label. Default is "Percentage (%)".
#' @param colors A character vector of colors for the different levels of frequency_for_factor. If NULL (default), ROP_colors will be used when available, otherwise rainbow colors will be used.
#' @param legend_title A character string specifying the legend title. If NULL (default), the frequency_for_factor name will be used.
#' @param show_percentages A logical value indicating whether to display percentage labels on the bars. Default is TRUE.
#' @param text_size A numeric value specifying the size of percentage labels on the bars. Default is 3.
#' @param text_color A character string specifying the color of percentage labels. Default is "black".
#' @param show_summary_table A logical value indicating whether to print the summary table. Default is TRUE.
#' @param export_summary_table A logical value indicating whether to export the summary table to the main environment. Default is FALSE.
#' @param name_summary_table A character string specifying the name of the exported summary table. Default is "summary_table_FactorStackedBarplot".
#' @details This function generates stacked barplots where each bar represents a level of the split_by_factor, and the stacked segments show the percentage distribution of the frequency_for_factor levels within each category. The function automatically converts extracted factors to factor type using as.factor(romicsExtractFactor()). Only complete cases (non-NA values) are included in the analysis. Colors are taken from ROP_colors by default when sufficient colors are available.
#' @return A ggplot object representing the stacked barplot. If export_summary_table is TRUE, the summary table is attached as an attribute (accessible via attr(plot, "summary_table")).
#' @author Geremy Clair
#' @export
romicsFactorStackedBarplot <- function(romics_object,
                                       split_by_factor,
                                       frequency_for_factor,
                                       title = NULL,
                                       xlab = NULL,
                                       ylab = "Percentage (%)",
                                       colors = NULL,
                                       legend_title = NULL,
                                       show_percentages = TRUE,
                                       text_size = 3,
                                       text_color = "black",
                                       show_summary_table = TRUE,
                                       export_summary_table = FALSE,
                                       name_summary_table = "summary_table_FactorStackedBarplot") {

  # Check for required packages
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("The dplyr package is required for this function. Please install it with: install.packages('dplyr')")
  }

  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("The tidyr package is required for this function. Please install it with: install.packages('tidyr')")
  }

  if (!is.romicsObject(romics_object)) stop("romics_object must be a romics object")
  if (missing(split_by_factor)      || !is.character(split_by_factor))      stop("split_by_factor must be a character string")
  if (missing(frequency_for_factor) || !is.character(frequency_for_factor)) stop("frequency_for_factor must be a character string")

  tryCatch({
    split_factor <- as.factor(romicsExtractFactor(romics_object, split_by_factor))
    freq_factor  <- as.factor(romicsExtractFactor(romics_object, frequency_for_factor))
  }, error = function(e) {
    stop("Error extracting factors. Please verify that the factor names exist in the romics_object metadata.")
  })

  if (length(split_factor) != length(freq_factor)) stop("The extracted factors must have the same length")

  plot_data <- data.frame(
    split_by      = as.character(split_factor),
    frequency_for = as.character(freq_factor),
    stringsAsFactors = FALSE
  )
  plot_data <- plot_data[complete.cases(plot_data), ]
  if (nrow(plot_data) == 0) stop("No complete cases found after removing NA values")

  split_levels        <- sort(unique(plot_data$split_by))
  factor_levels_alpha <- sort(unique(plot_data$frequency_for))        # alphabetical: for color assignment and legend
  factor_levels_stack <- rev(factor_levels_alpha)                     # reversed: for stacking order in the bars

  # ── Summary table ─────────────────────────────────────────────────────────────
  summary_data <- plot_data %>%
    dplyr::group_by(split_by, frequency_for) %>%
    dplyr::summarise(count = n(), .groups = "drop") %>%
    dplyr::group_by(split_by) %>%
    dplyr::mutate(total = sum(count), percentage = (count / total) * 100) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(split_by, frequency_for)

  summary_table <- as.data.frame(summary_data[, c("split_by", "frequency_for", "count", "percentage")])
  colnames(summary_table) <- c(split_by_factor, frequency_for_factor, "Count", "Percentage (%)")

  if (show_summary_table) {
    cat("Summary Statistics:\n==================\n")
    print(summary_table, row.names = FALSE)
    cat("\n")
  }
  # Store summary table for later attachment to plot
  summary_table_to_attach <- summary_table

  # ── Defaults ──────────────────────────────────────────────────────────────────
  if (is.null(title))        title        <- paste("Distribution of", frequency_for_factor, "by", split_by_factor)
  if (is.null(xlab))         xlab         <- split_by_factor
  if (is.null(legend_title)) legend_title <- frequency_for_factor

  # ── Colors ────────────────────────────────────────────────────────────────────
  # Colors are assigned in alphabetical order so user-supplied colors map
  # predictably: colors[1] -> first level alphabetically, colors[2] -> second, etc.
  n_levels <- length(factor_levels_alpha)
  make_default_colors <- function(n) {
    if (n <= 9) RColorBrewer::brewer.pal(max(3, n), "Set1")[seq_len(n)]
    else        colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(n)
  }
  if (!is.null(colors)) {
    if (length(colors) < n_levels) {
      warning(paste0("Not enough colors provided (", length(colors), " given, ",
                     n_levels, " needed). Falling back to default colors."))
      colors <- make_default_colors(n_levels)
    } else {
      colors <- colors[seq_len(n_levels)]
    }
  } else {
    colors <- make_default_colors(n_levels)
  }
  color_map <- setNames(colors, factor_levels_alpha)   # named by alphabetical levels

  # ── Build complete plot data ───────────────────────────────────────────────────
  # Use factor_levels_stack for the factor levels so bars stack in reversed order,
  # but color_map names still match because scale_fill_manual matches by name.
  plot_summary_data <- plot_data %>%
    dplyr::group_by(split_by, frequency_for) %>%
    dplyr::summarise(count = n(), .groups = "drop") %>%
    dplyr::group_by(split_by) %>%
    dplyr::mutate(total = sum(count), percentage = (count / total) * 100) %>%
    dplyr::ungroup() %>%
    tidyr::complete(split_by      = split_levels,
             frequency_for = factor_levels_stack,
             fill = list(count = 0L, total = 0, percentage = 0)) %>%
    dplyr::mutate(
      split_by      = factor(split_by,      levels = split_levels),
      frequency_for = factor(frequency_for, levels = factor_levels_stack)
    ) %>%
    dplyr::filter(percentage > 0) %>%
    dplyr::arrange(split_by, frequency_for)

  # ── Plot ──────────────────────────────────────────────────────────────────────
  p <- ggplot(plot_summary_data,
              aes(x = split_by, y = percentage, fill = frequency_for)) +
    geom_bar(stat = "identity", position = position_stack(reverse = TRUE),
             color = "white", linewidth = 0.5) +
    scale_fill_manual(
      values = color_map,
      name   = legend_title,
      breaks = factor_levels_alpha,   # legend displayed in alphabetical order
      drop   = FALSE
    ) +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    labs(title = title, x = xlab, y = ylab) +
    theme_minimal() +
    theme(
      plot.title         = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x        = element_text(angle = 45, hjust = 1),
      axis.title         = element_text(size = 12),
      legend.title       = element_text(size = 11, face = "bold"),
      legend.text        = element_text(size = 10),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank()
    )

  # ── Percentage labels ─────────────────────────────────────────────────────────
  if (show_percentages) {
    label_data <- plot_summary_data %>%
      dplyr::arrange(split_by, frequency_for) %>%
      dplyr::group_by(split_by) %>%
      dplyr::mutate(
        stack_top    = cumsum(percentage),
        stack_bottom = stack_top - percentage,
        label_y      = stack_bottom + percentage / 2
      ) %>%
      dplyr::ungroup() %>%
      dplyr::filter(percentage > 3) %>%
      dplyr::mutate(label_text = paste0(round(percentage, 1), "%"))

    if (nrow(label_data) > 0) {
      p <- p + geom_text(
        data        = label_data,
        aes(x = split_by, y = label_y, label = label_text),
        color       = text_color,
        size        = text_size,
        fontface    = "bold",
        inherit.aes = FALSE,
        vjust       = 0.5,
        hjust       = 0.5
      )
    }
  }

  # Attach summary table if export was requested
  if (export_summary_table) {
    attr(p, "summary_table") <- summary_table_to_attach
    attr(p, "summary_table_name") <- name_summary_table
    # Also assign to the global environment with the specified name
    assign(name_summary_table, summary_table_to_attach, envir = .GlobalEnv)
    message("Summary table '", name_summary_table, "' has been assigned to the global environment")
  }

  return(p)
}

#' romicsFactorSummaryTable()
#' @description Generates a summary table showing the distribution of samples across two factor levels.
#' This is the data table version of romicsFactorStackedBarplot, useful for exporting or further analysis.
#' @param romics_object A romics_object created using romicsCreateObject().
#' @param split_by_factor A character vector (length=1) containing a factor of the romics_object that will define the rows/categories. The list of factors can be identified by using the function romicsFactorNames()
#' @param frequency_for_factor A character vector (length=1) containing a factor of the romics_object that will be used to calculate frequency distributions within each level of split_by_factor. The list of factors can be identified by using the function romicsFactorNames()
#' @param percentage Logical. If TRUE, includes percentage columns in addition to counts. Default: TRUE.
#' @param wide_format Logical. If TRUE, returns data in wide format (split_by_factor levels as rows, frequency_for_factor levels as columns).
#'        If FALSE, returns long format (one row per split_by_factor-frequency_for_factor combination). Default: FALSE.
#' @param include_totals Logical. If TRUE and wide_format = TRUE, includes a "Total" column with row totals. Default: FALSE.
#' @param round_percentages Numeric. Number of decimal places for percentage values. Default: 1.
#' @return A data.frame containing the summary table in either long or wide format.
#' @author Geremy Clair
#' @export
romicsFactorSummaryTable <- function(romics_object,
                                     split_by_factor,
                                     frequency_for_factor,
                                     percentage = TRUE,
                                     wide_format = FALSE,
                                     include_totals = FALSE,
                                     round_percentages = 1) {

  # Load required library
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package is required for this function")
  }

  # Input validation
  if (!is.list(romics_object) || is.null(romics_object$metadata)) {
    stop("romics_object must be a list containing a 'metadata' component")
  }

  if(missing(split_by_factor) || !is.character(split_by_factor) || length(split_by_factor) != 1) {
    stop("<split_by_factor> must be a single character string")
  }

  if(missing(frequency_for_factor) || !is.character(frequency_for_factor) || length(frequency_for_factor) != 1) {
    stop("<frequency_for_factor> must be a single character string")
  }

  if(missing(percentage)) {percentage <- TRUE}
  if(missing(wide_format)) {wide_format <- FALSE}
  if(missing(include_totals)) {include_totals <- FALSE}
  if(missing(round_percentages)) {round_percentages <- 1}

  # Extract factors and convert to factor type
  tryCatch({
    split_factor <- as.factor(romicsExtractFactor(romics_object, split_by_factor))
    freq_factor <- as.factor(romicsExtractFactor(romics_object, frequency_for_factor))
  }, error = function(e) {
    stop("Error extracting factors. Please check that the factor names exist in the romics_object metadata.")
  })

  # Check if factors have the same length
  if (length(split_factor) != length(freq_factor)) {
    stop("The extracted factors must have the same length")
  }

  # Create data frame for analysis
  plot_data <- data.frame(
    split_by = split_factor,
    frequency_for = freq_factor,
    stringsAsFactors = FALSE
  )

  # Remove rows with NA values
  plot_data <- plot_data[complete.cases(plot_data), ]

  if (nrow(plot_data) == 0) {
    stop("No complete cases found after removing NA values")
  }

  # Calculate counts and percentages
  summary_data <- plot_data %>%
    dplyr::group_by(split_by, frequency_for) %>%
    dplyr::summarise(count = n(), .groups = "drop") %>%
    dplyr::group_by(split_by) %>%
    dplyr::mutate(
      total = sum(count),
      percentage = (count / total) * 100
    ) %>%
    dplyr::ungroup()

  # Round percentages if requested
  if(percentage && !is.null(round_percentages)) {
    summary_data$percentage <- round(summary_data$percentage, round_percentages)
  }

  # Rename columns to use original factor names
  colnames(summary_data)[colnames(summary_data) == "split_by"] <- split_by_factor
  colnames(summary_data)[colnames(summary_data) == "frequency_for"] <- frequency_for_factor

  # Return long format if requested
  if(!wide_format) {
    if(!percentage) {
      # Remove percentage and total columns if not requested
      summary_data <- summary_data[, c(split_by_factor, frequency_for_factor, "count")]
    }
    return(summary_data)
  }

  # Convert to wide format
  if(percentage) {
    # Create separate wide tables for counts and percentages
    count_wide <- tidyr::pivot_wider(
      summary_data[, c(split_by_factor, frequency_for_factor, "count")],
      names_from = all_of(frequency_for_factor),
      values_from = count,
      values_fill = 0
    )

    perc_wide <- tidyr::pivot_wider(
      summary_data[, c(split_by_factor, frequency_for_factor, "percentage")],
      names_from = all_of(frequency_for_factor),
      values_from = percentage,
      values_fill = 0
    )

    # Combine count and percentage columns
    # For each frequency_for_factor level, create "level (count, %)" format
    freq_levels <- unique(summary_data[[frequency_for_factor]])

    result_wide <- count_wide[, split_by_factor, drop = FALSE]

    for(level in freq_levels) {
      if(level %in% colnames(count_wide) && level %in% colnames(perc_wide)) {
        counts <- count_wide[[level]]
        percs <- perc_wide[[level]]
        # Create combined column with format "count (percentage%)"
        combined_col <- paste0(counts, " (", percs, "%)")
        result_wide[[level]] <- combined_col
      }
    }

    wide_result <- result_wide

  } else {
    # Just counts in wide format
    wide_result <- tidyr::pivot_wider(
      summary_data[, c(split_by_factor, frequency_for_factor, "count")],
      names_from = all_of(frequency_for_factor),
      values_from = count,
      values_fill = 0
    )
  }

  # Add totals if requested
  if(include_totals && wide_format) {
    # Get row totals from original summary_data
    totals <- summary_data %>%
      dplyr::group_by(.data[[split_by_factor]]) %>%
      dplyr::summarise(Total = first(total), .groups = "drop")

    # Merge with wide result
    wide_result <- merge(wide_result, totals, by = split_by_factor, all.x = TRUE)
  }

  return(wide_result)
}

#' romicsSplitFactor
#' @description This function will split different elements contained in a factor (separated by a common separating character string)
#' @param romics_object A romics_object created using romicsCreateObject().
#' @param factor A character vector (lenght=1) containing a factors of the romics_object. the list of factors can be identified by using the function romicsFactorNames()
#' @param  splitting_string A character vector (lenght=1) that indicates the separator of the the different factors to be used (when multiple factors to be created are contained in a given factor)
#' @details generate new factor from an existing factor, it will indicates if the new factors are <present> or <absent> from the original factor, if a splitting string is used, then it will create additional factors based on the original factor content.
#' @return a Romics_object with new factors
#' @author Geremy Clair, Naina Beishembieva
#' @export
#'
romicsSplitFactor<-function(romics_object=romics_object,factor=f,splitting_string=";"){
  arguments<-as.list(match.call())
  if(missing(splitting_string)){splitting_string<-";"}
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("<romics_object> is missing or is not in the appropriate format")}
  if(!factor %in% romicsFactorNames(romics_object)){stop("Your <factor> is not in the list of factors of the romics object, you can use the function romicsFactorNames to identify those.")}
  if(!is.character(splitting_string) & length(splitting_string)!=1){
  stop("the <splitting_string> needs to be a character vector of length 1.")
  }

  m<-romics_object$metadata
  f<-romicsExtractFactor(romics_object,factor = factor)
  f<-as.character(f)
  l<-strsplit(x = f,split = splitting_string)
  new_factor_list<-unlist(l)
  unique_new_f<-unique(new_factor_list)

  if(sum(unique_new_f %in% rownames(m))!=0){
    print("One or more of the newly created factors were already factors of the <romics_object>")
    print(unique_new_f[unique_new_f %in% rownames(m)])
    stop()
  }

  df <- matrix("Absent", ncol = ncol(m), nrow = length(unique_new_f))

  # Vectorized approach to fill the new dataframe `df`
  df <- sapply(1:ncol(m), function(j) {
    sapply(unique_new_f, function(uniq_val) {
      ifelse(uniq_val %in% l[[j]], uniq_val, "Absent")
    })
  })

  rownames(df) <- unique_new_f
  colnames(df) <- colnames(m)

  romics_object$metadata<-rbind(m,df)

  romics_object<-romicsUpdateSteps(romics_object,arguments)
  return(romics_object)
}

#' romicsCombineFactors
#' @description This function will use multiple factor of an romics object to create a new one.
#' @param romics_object A romics_object created using romicsCreateObject().
#' @param list_of_factors A character vector (length=>1) containing names of factors of the romics_object in the desired order. The list of factor names can be identified by using the function romicsFactorNames().
#' @param separator A character vector (length=1) that indicates the separator of the different factors to be used (when multiple factors to be created are contained in a given factor)
#' @param replace_missing A character string to replace missing values. Default: "unknown"
#' @param factor_name A character string for the name of the new factor. If not provided, uses the collapsed names of input factors.
#' @details Generate new factor from a list of existing factors in the specified order. The new factor will be the aggregation of the different factors content in the order provided.
#' @return A Romics_object with a new factor. The name of the factor will be the collapsed name of the factors aggregated (separated with the <separator>), the content of the factor will be the aggregation of the different factors content in the specified order.
#' @author Geremy Clair, Naina Beishembieva
#' @export
romicsCombineFactors <- function(romics_object = romics_object,
                                 list_of_factors = c("f1", "f2"),
                                 separator = ";",
                                 replace_missing = "unknown",
                                 factor_name = NULL) {

  arguments <- as.list(match.call())

  # Set defaults
  if (missing(separator)) {
    separator <- ";"
  }
  if (missing(replace_missing)) {
    replace_missing <- "unknown"
  }

  # Validation
  if (!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("<romics_object> is missing or is not in the appropriate format")
  }

  if ((!is.character(list_of_factors)) |
      (sum(list_of_factors %in% romicsFactorNames(romics_object)) != length(list_of_factors))) {
    stop("Your <list_of_factors> is not a character vector or is not part of the list of factors of the romics object. You can use the function romicsFactorNames() to identify the names of the factors.")
  }

  # Determine new factor name
  if (is.null(factor_name)) {
    new_factor_name <- paste(as.character(list_of_factors), collapse = separator)
  } else {
    new_factor_name <- factor_name
  }

  if (new_factor_name %in% romicsFactorNames(romics_object)) {
    stop("Your <new factor> is already part of the list of factors. You can use the function romicsFactorNames() to identify the names of the factors.")
  }

  # Get metadata
  m <- romics_object$metadata

  # FIXED: Subset in the ORDER specified by list_of_factors
  # This ensures the order is preserved
  m_subset <- m[list_of_factors, , drop = FALSE]

  # Handle missing values if specified
  if (!is.null(replace_missing)) {
    m_subset[is.na(m_subset)] <- replace_missing
  }

  # Combine factors in the correct order
  # Apply column-wise, combining rows in the order they appear
  new_f <- apply(m_subset, 2, function(col) paste(col, collapse = separator))
  new_f <- as.data.frame(t(new_f))
  rownames(new_f) <- new_factor_name

  # Add to metadata
  romics_object$metadata <- rbind(m, new_f)

  # Update steps
  romics_object <- romicsUpdateSteps(romics_object, arguments)

  return(romics_object)
}

#' romicsRemoveFactor
#' @description This function removes one or more factors from the metadata of a romics object.
#' @param romics_object A romics_object created using romicsCreateObject().
#' @param factors_to_remove A character vector containing the name(s) of factor(s) to remove. The list of factor names can be identified by using the function romicsFactorNames().
#' @param force A logical indicating whether to allow removal of the main factor. Default: FALSE (prevents accidental removal of main factor).
#' @details This function removes specified factors from the metadata layer of a romics object. By default, it prevents removal of the main factor to avoid breaking the object structure. Set force=TRUE to override this protection.
#' @return A romics_object with the specified factor(s) removed from the metadata.
#' @author Geremy Clair
#' @export
#' @examples \dontrun{
#' # Remove a single factor
#' romics_object <- romicsRemoveFactor(romics_object, factors_to_remove = "old_factor")
#'
#' # Remove multiple factors
#' romics_object <- romicsRemoveFactor(romics_object, factors_to_remove = c("factor1", "factor2"))
#'
#' # Remove main factor (use with caution!)
#' romics_object <- romicsRemoveFactor(romics_object, factors_to_remove = "condition", force = TRUE)
#' }
romicsRemoveFactor <- function(romics_object,
                               factors_to_remove,
                               force = FALSE) {

  arguments <- as.list(match.call())

  # Input validation
  if (!is.romicsObject(romics_object) || missing(romics_object)) {
    stop("<romics_object> is missing or is not in the appropriate format.")
  }

  if (missing(factors_to_remove)) {
    stop("<factors_to_remove> is required. Please specify one or more factor names to remove.")
  }

  if (!is.character(factors_to_remove)) {
    stop("<factors_to_remove> must be a character vector containing factor name(s).")
  }

  if (missing(force)) {
    force <- FALSE
  }

  # Get current factors
  current_factors <- romicsFactorNames(romics_object)

  if (length(current_factors) == 0) {
    stop("No factors found in the romics object metadata.")
  }

  # Check if all factors to remove exist
  missing_factors <- factors_to_remove[!factors_to_remove %in% current_factors]

  if (length(missing_factors) > 0) {
    stop(paste0("The following factor(s) were not found in the romics object: ",
                paste(missing_factors, collapse = ", "),
                "\nAvailable factors: ",
                paste(current_factors, collapse = ", ")))
  }

  # Check if trying to remove main factor without force
  if (!force && romics_object$main_factor %in% factors_to_remove) {
    stop(paste0("You are trying to remove the main factor ('",
                romics_object$main_factor,
                "'). This may break some functions that depend on it.\n",
                "If you really want to remove it, set force = TRUE.\n",
                "Consider using romicsChangeMainFactor() to set a new main factor first."))
  }

  # Check if removing all factors
  if (length(factors_to_remove) == length(current_factors)) {
    warning("You are removing ALL factors from the romics object. This may cause issues with downstream analyses.")
  }

  # Get metadata
  metadata <- romics_object$metadata

  # Remove specified factors
  factors_to_keep <- rownames(metadata)[!rownames(metadata) %in% factors_to_remove]

  if (length(factors_to_keep) == 0) {
    romics_object$metadata <- data.frame(matrix(ncol = ncol(metadata), nrow = 0))
    colnames(romics_object$metadata) <- colnames(metadata)
    message("All factors have been removed. The metadata now contains only column identifiers.")
  } else {
    romics_object$metadata <- metadata[factors_to_keep, , drop = FALSE]
  }

  # If main factor was removed (with force=TRUE), set new main factor or clear it
  if (romics_object$main_factor %in% factors_to_remove) {
    if (length(factors_to_keep) > 0) {
      romics_object$main_factor <- factors_to_keep[1]
      message(paste0("Main factor was removed. New main factor set to: '",
                     romics_object$main_factor, "'"))
    } else {
      romics_object$main_factor <- NULL
      message("Main factor was removed and no other factors remain.")
    }
  }

  # Update steps
  romics_object <- romicsUpdateSteps(romics_object, arguments)

  # Print summary
  n_removed <- length(factors_to_remove)
  message(paste0("Successfully removed ", n_removed, " factor(s): ",
                 paste(factors_to_remove, collapse = ", ")))

  remaining_factors <- romicsFactorNames(romics_object)
  if (length(remaining_factors) > 0) {
    message(paste0("Remaining factors (", length(remaining_factors), "): ",
                   paste(remaining_factors, collapse = ", ")))
  }

  return(romics_object)
}

#' @title Create Metadata Factor from Feature
#' @description Creates a new factor in the metadata based on a feature's intensity using various classification methods.
#' @param romics_object A romics_object created using romicsCreateObject().
#' @param feature A character string specifying the name of the feature to use. Must be present in rownames of the data layer.
#' @param method Classification method: "presence" (present/absent), "value" (above/below a fixed value), "percentage" (above/below a percentile). Default: "presence".
#' @param value Numeric value for splitting when method = "value", Required when method = "value".
#' @param percentage Numeric between 0 and 100 for percentage method (e.g., 75 for top 25% vs bottom 75%). Default: 50.
#' @param n_groups Integer for number of groups when using percentage method with equal splits (2, 3, or 4). Default: 2. When > 2, overrides percentage parameter.
#' @param factor_name Character string for the name of the new factor. If NULL, auto-generated based on method. Default: NULL.
#' @param labels Character vector of custom labels. Length depends on n_groups. If NULL, auto-generated. For binary splits, order is c(high, low). Default: NULL.
#' @param na_label Character string for samples where the feature is NA/not detected. Default: "absent".
#' @param na_as_low Logical. For binary methods, should NA values be treated as "low" rather than a separate category? Default: FALSE.
#' @param strict_inequality Logical. Use strict inequality (>) or >= for high group. Default: FALSE (uses >).
#' @param include_value_in_name Logical. Include value/percentage in auto-generated factor name. Default: TRUE.
#' @param verbose Logical. Print detailed summary. Default: TRUE.
#' @return Returns the romics_object with the newly created factor added to its metadata.
#' @author Geremy Clair
#' @export
romicsFactorFromFeature <- function(romics_object,
                                    feature,
                                    method = "presence",
                                    value = NULL,
                                    percentage = 50,
                                    n_groups = 2,
                                    factor_name = NULL,
                                    labels = NULL,
                                    na_label = "absent",
                                    na_as_low = FALSE,
                                    strict_inequality = FALSE,
                                    include_value_in_name = TRUE,
                                    verbose = TRUE) {

  arguments <- as.list(match.call())

  # Input validation
  if (!is.romicsObject(romics_object) || missing(romics_object)) {
    stop("<romics_object> is missing or is not in the appropriate format.")
  }

  if (missing(feature) || !is.character(feature) || length(feature) != 1) {
    stop("<feature> must be a single character string specifying the feature name.")
  }

  if (missing(method)) {method <- "presence"}
  valid_methods <- c("presence", "value", "percentage")
  if (!method %in% valid_methods) {
    stop(paste("<method> must be one of:", paste(valid_methods, collapse = ", ")))
  }

  if (method == "value" && is.null(value)) {
    stop("When method = 'value', you must provide a <value> parameter.")
  }

  if (method == "percentage") {
    if (missing(percentage)) {percentage <- 50}
    if (missing(n_groups))   {n_groups   <- 2}

    if (!is.numeric(percentage) || percentage <= 0 || percentage >= 100) {
      stop("<percentage> must be between 0 and 100 (exclusive).")
    }

    if (!is.numeric(n_groups) || !n_groups %in% c(2, 3, 4)) {
      stop("<n_groups> must be 2, 3, or 4.")
    }
  }

  if (missing(na_as_low))             {na_as_low             <- FALSE}
  if (missing(strict_inequality))     {strict_inequality     <- FALSE}
  if (missing(include_value_in_name)) {include_value_in_name <- TRUE}
  if (missing(verbose))               {verbose               <- TRUE}

  # Check if feature exists
  if (!feature %in% rownames(romics_object$data)) {
    stop(paste0("Feature '", feature, "' not found in the data layer.\n",
                "Available features: ", paste(head(rownames(romics_object$data), 10), collapse = ", "),
                ifelse(nrow(romics_object$data) > 10, "...", "")))
  }

  # Extract feature values
  feature_values <- as.numeric(romics_object$data[feature, ])
  non_na_values  <- feature_values[!is.na(feature_values)]

  # ── Early exit: all values are NA (applies to all methods) ───────────────────
  if (length(non_na_values) == 0) {
    if (verbose) {
      message(paste0("All values are NA for feature '", feature, "'. Creating all-absent factor."))
    }
    if (is.null(na_label)) { na_label <- paste0(feature, "_absent") }
    if (is.null(factor_name)) {
      factor_name <- switch(method,
                            "presence"   = paste0(feature, "_presence"),
                            "value"      = paste0(feature, "_value"),
                            "percentage" = paste0(feature, "_top", percentage, "pct")
      )
    }
    new_factor <- factor(rep(na_label, length(feature_values)), levels = na_label)

    if (factor_name %in% romicsFactorNames(romics_object)) {
      stop(paste0("Factor '", factor_name, "' already exists in the metadata.\n",
                  "Please choose a different factor_name or remove the existing factor first using romicsRemoveFactor()."))
    }

    new_factor_row <- data.frame(matrix(as.character(new_factor), nrow = 1))
    colnames(new_factor_row) <- colnames(romics_object$metadata)
    rownames(new_factor_row) <- factor_name
    romics_object$metadata <- rbind(romics_object$metadata, new_factor_row)
    romics_object$metadata[factor_name, ] <- factor(as.character(romics_object$metadata[factor_name, ]),
                                                    levels = levels(new_factor))
    romics_object <- romicsUpdateSteps(romics_object, arguments)
    return(romics_object)
  }

  # Initialize variables
  cutpoint  <- NULL
  cutpoints <- NULL
  new_factor <- character(length(feature_values))

  # ── Process based on method ───────────────────────────────────────────────────
  if (method == "presence") {
    # Presence/Absence classification
    if (is.null(labels)) {
      labels <- c(paste0(feature, "_present"), paste0(feature, "_absent"))
    } else if (length(labels) != 2) {
      stop("For method = 'presence', <labels> must be a vector of length 2: c(present_label, absent_label)")
    }

    new_factor[!is.na(feature_values)] <- labels[1]
    new_factor[is.na(feature_values)]  <- labels[2]

    if (is.null(factor_name)) {
      factor_name <- paste0(feature, "_presence")
    }

    new_factor <- factor(new_factor, levels = labels)

  } else if (method == "value") {
    # Value-based classification (fixed cutpoint)
    cutpoint <- value

    if (is.null(labels)) {
      labels <- c(paste0(feature, "_high"), paste0(feature, "_low"))
    } else if (length(labels) != 2) {
      stop("For method = 'value', <labels> must be a vector of length 2: c(high_label, low_label)")
    }

    # Apply value cutoff
    new_factor[!is.na(feature_values) & feature_values >  value] <- labels[1]
    new_factor[!is.na(feature_values) & feature_values <= value] <- labels[2]

    # Handle NA values
    if (na_as_low) {
      new_factor[is.na(feature_values)] <- labels[2]
      factor_levels <- labels
    } else {
      if (is.null(na_label)) {na_label <- paste0(feature, "_not_detected")}
      new_factor[is.na(feature_values)] <- na_label
      factor_levels <- c(labels, na_label)
    }

    if (is.null(factor_name)) {
      if (include_value_in_name) {
        value_str   <- format(value, scientific = FALSE, drop0trailing = TRUE)
        factor_name <- paste0(feature, "_value_", value_str)
      } else {
        factor_name <- paste0(feature, "_value")
      }
    }

    new_factor <- factor(new_factor, levels = factor_levels)

  } else if (method == "percentage") {

    if (n_groups == 2) {
      # Binary split: top X% vs rest — cutpoint is the (1 - percentage/100) quantile
      cutpoint <- quantile(non_na_values, probs = 1 - (percentage/100), na.rm = TRUE)

      if (is.null(labels)) {
        labels <- c(paste0(feature, "_high"),
                    paste0(feature, "_low"))
      } else if (length(labels) != 2) {
        stop("For n_groups = 2, <labels> must be a vector of length 2: c(high_label, low_label)")
      }

      new_factor[!is.na(feature_values) & feature_values >  cutpoint] <- labels[1]
      new_factor[!is.na(feature_values) & feature_values <= cutpoint] <- labels[2]

      # Handle NA values
      if (na_as_low) {
        new_factor[is.na(feature_values)] <- labels[2]
        factor_levels <- labels
      } else {
        if (is.null(na_label)) {na_label <- paste0(feature, "_absent")}
        new_factor[is.na(feature_values)] <- na_label
        factor_levels <- c(labels, na_label)
      }

      if (is.null(factor_name)) {
        factor_name <- paste0(feature, "_top", percentage, "pct")
      }

      new_factor <- factor(new_factor, levels = factor_levels)

    } else if (n_groups == 3) {
      # Tertile split
      cutpoints <- quantile(non_na_values, probs = c(1/3, 2/3), na.rm = TRUE)

      if (is.null(labels)) {
        labels <- c(paste0(feature, "_high"), paste0(feature, "_medium"), paste0(feature, "_low"))
      } else if (length(labels) != 3) {
        stop("For n_groups = 3, <labels> must be a vector of length 3: c(high_label, medium_label, low_label)")
      }

      new_factor[!is.na(feature_values) & feature_values >  cutpoints[2]] <- labels[1]
      new_factor[!is.na(feature_values) & feature_values >  cutpoints[1] & feature_values <= cutpoints[2]] <- labels[2]
      new_factor[!is.na(feature_values) & feature_values <= cutpoints[1]] <- labels[3]

      # Handle NA values
      if (na_as_low) {
        new_factor[is.na(feature_values)] <- labels[3]
        factor_levels <- labels
      } else {
        if (is.null(na_label)) {na_label <- paste0(feature, "_not_detected")}
        new_factor[is.na(feature_values)] <- na_label
        factor_levels <- c(labels, na_label)
      }

      if (is.null(factor_name)) {
        factor_name <- paste0(feature, "_tertile")
      }

      new_factor <- factor(new_factor, levels = factor_levels)

    } else if (n_groups == 4) {
      # Quartile split
      cutpoints <- quantile(non_na_values, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)

      if (is.null(labels)) {
        labels <- c(paste0(feature, "_Q4"), paste0(feature, "_Q3"),
                    paste0(feature, "_Q2"), paste0(feature, "_Q1"))
      } else if (length(labels) != 4) {
        stop("For n_groups = 4, <labels> must be a vector of length 4: c(Q4_label, Q3_label, Q2_label, Q1_label)")
      }

      new_factor[!is.na(feature_values) & feature_values >  cutpoints[3]] <- labels[1]
      new_factor[!is.na(feature_values) & feature_values >  cutpoints[2] & feature_values <= cutpoints[3]] <- labels[2]
      new_factor[!is.na(feature_values) & feature_values >  cutpoints[1] & feature_values <= cutpoints[2]] <- labels[3]
      new_factor[!is.na(feature_values) & feature_values <= cutpoints[1]] <- labels[4]

      # Handle NA values
      if (na_as_low) {
        new_factor[is.na(feature_values)] <- labels[4]
        factor_levels <- labels
      } else {
        if (is.null(na_label)) {na_label <- paste0(feature, "_not_detected")}
        new_factor[is.na(feature_values)] <- na_label
        factor_levels <- c(labels, na_label)
      }

      if (is.null(factor_name)) {
        factor_name <- paste0(feature, "_quartile")
      }

      new_factor <- factor(new_factor, levels = factor_levels)
    }
  }

  # Check if factor already exists
  if (factor_name %in% romicsFactorNames(romics_object)) {
    stop(paste0("Factor '", factor_name, "' already exists in the metadata.\n",
                "Please choose a different factor_name or remove the existing factor first using romicsRemoveFactor()."))
  }

  # Add to metadata
  new_factor_row <- data.frame(matrix(as.character(new_factor), nrow = 1))
  colnames(new_factor_row) <- colnames(romics_object$metadata)
  rownames(new_factor_row) <- factor_name

  romics_object$metadata <- rbind(romics_object$metadata, new_factor_row)

  # Ensure the new row is factor type
  romics_object$metadata[factor_name, ] <- factor(as.character(romics_object$metadata[factor_name, ]),
                                                  levels = levels(new_factor))

  # Update steps
  romics_object <- romicsUpdateSteps(romics_object, arguments)

  # Print summary
  if (verbose) {
    factor_table <- table(new_factor, useNA = "ifany")

    message(paste0("Factor '", factor_name, "' created successfully"))
    message(paste0("Method: ", method))

    if (method == "value") {
      message(paste0("Cutpoint: ", format(cutpoint, scientific = FALSE, digits = 3)))
    } else if (method == "percentage" && n_groups == 2) {
      message(paste0("Percentile: ", percentage, "% (cutpoint: ",
                     format(cutpoint, scientific = FALSE, digits = 3), ")"))
    } else if (method == "percentage" && n_groups > 2) {
      message(paste0("Groups: ", n_groups, " equal groups"))
      message(paste0("Cutpoints: ", paste(format(cutpoints, scientific = FALSE, digits = 3), collapse = ", ")))
    }

    message("Summary:")
    for (i in 1:length(factor_table)) {
      level_name <- names(factor_table)[i]
      if (is.na(level_name)) level_name <- "<NA>"
      message(paste0("  ", level_name, ": n=", factor_table[i]))
    }
  }

  return(romics_object)
}


#' romicsTransferFactor()
#' @description Transfer one or more metadata factors from one romics_object to another.
#' The function ensures that factors are properly aligned based on sample names.
#' If the source object is a subset, samples present only in the target will be
#' labeled as "not_present_in_subset".
#' @param source_romics_object A romics_object containing the factors to be transferred
#' @param target_romics_object A romics_object that will receive the transferred factors
#' @param factors Character vector specifying the names of factors to transfer from source to target
#' @param overwrite Boolean. If TRUE, existing factors in target will be overwritten. Default: FALSE
#' @param allow_partial Boolean. If TRUE, allows transfer even when source is a subset of target. Default: FALSE
#' @param missing_label Character string to use for samples not present in source. Default: "not_present_in_subset"
#' @return Returns the target_romics_object with the transferred factors added to its metadata
#' @details This function transfers metadata factors between romics_objects. By default, both objects
#' must contain the same samples. When allow_partial = TRUE, the source can be a subset of the target,
#' and missing samples will be labeled with missing_label.
#' @author Geremy Clair
#' @export
romicsTransferFactor <- function(source_romics_object, target_romics_object, factors,
                                 overwrite = FALSE, allow_partial = FALSE,
                                 missing_label = "not_present_in_subset") {
  arguments <- as.list(match.call())

  # Input validation
  if(!is.romicsObject(source_romics_object) | missing(source_romics_object)) {
    stop("source_romics_object is missing or is not in the appropriate format.")
  }
  if(!is.romicsObject(target_romics_object) | missing(target_romics_object)) {
    stop("target_romics_object is missing or is not in the appropriate format.")
  }
  if(missing(factors)) {
    stop("factors parameter is required. It should specify which factors to transfer.")
  }

  # Convert single factor to vector if needed
  if(!is.vector(factors)) {
    factors <- c(factors)
  }

  # Check if specified factors exist in source
  source_factors <- rownames(source_romics_object$metadata)
  missing_factors <- factors[!factors %in% source_factors]
  if(length(missing_factors) > 0) {
    stop(paste("The following factors were not found in the source_romics_object:",
               paste(missing_factors, collapse = ", ")))
  }

  # Get sample names
  source_samples <- colnames(source_romics_object$metadata)
  target_samples <- colnames(target_romics_object$metadata)

  # Check sample compatibility
  samples_in_source_not_target <- setdiff(source_samples, target_samples)
  samples_in_target_not_source <- setdiff(target_samples, source_samples)

  if(length(samples_in_source_not_target) > 0) {
    stop(paste("The source object contains samples not present in the target object:",
               paste(samples_in_source_not_target, collapse = ", ")))
  }

  # Handle partial transfer case
  if(length(samples_in_target_not_source) > 0) {
    if(!allow_partial) {
      stop(paste("The target object contains samples not present in the source object:",
                 paste(samples_in_target_not_source, collapse = ", "),
                 "\nSet allow_partial = TRUE to transfer factors anyway.",
                 "\nSamples not in source will be labeled as:", missing_label))
    } else {
      message(paste(length(samples_in_target_not_source),
                    "sample(s) in target are not present in source and will be labeled as:",
                    missing_label))
    }
  }

  # Check for existing factors in target and handle based on overwrite parameter
  target_factors <- rownames(target_romics_object$metadata)
  existing_factors <- factors[factors %in% target_factors]
  if(length(existing_factors) > 0 && !overwrite) {
    stop(paste("The following factors already exist in the target_romics_object:",
               paste(existing_factors, collapse = ", "),
               "\nSet overwrite = TRUE to replace them."))
  } else if(length(existing_factors) > 0) {
    warning(paste("The following factors will be overwritten in the target_romics_object:",
                  paste(existing_factors, collapse = ", ")))
  }

  # ── Fast transfer: build all new rows at once then do a single rbind ──────────

  # Extract all requested factors from source as a character matrix (factors x source_samples)
  source_block <- source_romics_object$metadata[factors, , drop = FALSE]

  # Build aligned matrix (factors x target_samples) in one vectorised operation
  # Columns present in source are copied; missing columns get missing_label
  n_factors <- length(factors)
  n_targets <- length(target_samples)

  # Pre-fill with missing_label
  aligned_matrix <- matrix(missing_label,
                           nrow     = n_factors,
                           ncol     = n_targets,
                           dimnames = list(factors, target_samples))

  # Fill in the samples that exist in source using index matching (no loop needed)
  shared_samples <- intersect(target_samples, source_samples)
  if(length(shared_samples) > 0) {
    aligned_matrix[, shared_samples] <- as.matrix(source_block[, shared_samples, drop = FALSE])
  }

  # Remove any factors that already exist in target (overwrite case)
  if(length(existing_factors) > 0) {
    target_romics_object$metadata <- target_romics_object$metadata[
      !rownames(target_romics_object$metadata) %in% existing_factors, , drop = FALSE]
  }

  # Single rbind for all factors at once
  target_romics_object$metadata <- rbind(target_romics_object$metadata, aligned_matrix)

  # Update steps
  target_romics_object <- romicsUpdateSteps(target_romics_object, arguments)

  message(paste(length(factors), "factor(s) transferred successfully."))
  return(target_romics_object)
}


#' romicsChangeLevelName()
#' @description Changes the name(s) of one or more levels in a factor of the romics_object metadata layer.
#' @param romics_object A romics_object created using romicsCreateObject().
#' @param factor_name Character string. The name of the factor whose level names should be changed.
#' @param old_names Character vector. The current name(s) of the level(s) to be changed.
#' @param new_names Character vector. The new name(s) for the level(s). Must be the same length as old_names.
#' @details This function allows renaming of factor levels in the metadata layer of a romics_object.
#' @return A romics_object with updated factor level names in the metadata layer.
#' @author Geremy Clair
#' @export
romicsChangeLevelName <- function(romics_object,
                                  factor_name,
                                  old_names,
                                  new_names) {

  arguments <- as.list(match.call())

  # Input validation
  if(!is.romicsObject(romics_object) || missing(romics_object)) {
    stop("<romics_object> is missing or is not in the appropriate format")
  }

  if(missing(factor_name) || !is.character(factor_name) || length(factor_name) != 1) {
    stop("<factor_name> must be a single character string")
  }

  if(missing(old_names) || !is.character(old_names) || length(old_names) < 1) {
    stop("<old_names> must be a character vector with at least one element")
  }

  if(missing(new_names) || !is.character(new_names) || length(new_names) < 1) {
    stop("<new_names> must be a character vector with at least one element")
  }

  # Check if factor exists
  if(!factor_name %in% rownames(romics_object$metadata)) {
    stop(paste0("Factor '", factor_name, "' not found in metadata.\n",
                "Available factors: ", paste(romicsFactorNames(romics_object), collapse = ", ")))
  }

  # Validate new_names length
  if(length(new_names) != 1 && length(new_names) != length(old_names)) {
    stop(paste0("<new_names> must have length 1 (collapse all) or same length as <old_names> (",
                length(old_names), ")"))
  }

  # Extract current factor values
  factor_values <- as.character(romics_object$metadata[factor_name, ])
  unique_levels <- unique(factor_values)

  # Check if all old_names exist in the factor
  missing_names <- old_names[!old_names %in% unique_levels]
  if(length(missing_names) > 0) {
    stop(paste0("The following level names were not found in factor '", factor_name, "': ",
                paste(missing_names, collapse = ", "), "\n",
                "Available levels: ", paste(unique_levels, collapse = ", ")))
  }

  # Check for conflicts (new names that already exist and aren't being replaced)
  if(length(new_names) > 1) {
    # One-to-one mapping
    for(i in seq_along(new_names)) {
      if(new_names[i] %in% unique_levels && !new_names[i] %in% old_names) {
        warning(paste0("New name '", new_names[i], "' already exists as a level in factor '",
                       factor_name, "'. This may cause confusion."))
      }
    }
  } else {
    # Many-to-one mapping
    if(new_names %in% unique_levels && !new_names %in% old_names) {
      message(paste0("New name '", new_names, "' already exists as a level. ",
                     "The specified old levels will be merged with this existing level."))
    }
  }

  # Perform the renaming using vectorized operations
  new_factor_values <- factor_values

  if(length(new_names) == 1) {
    # Many-to-one: collapse multiple levels to one name
    new_factor_values[factor_values %in% old_names] <- new_names
  } else {
    # One-to-one: map each old name to corresponding new name
    lookup <- setNames(new_names, old_names)
    to_change <- factor_values %in% old_names
    new_factor_values[to_change] <- lookup[factor_values[to_change]]
  }

  # Get new unique levels (maintaining order of first appearance)
  new_unique_levels <- unique(new_factor_values)

  # Update the metadata with new factor values as a proper factor
  romics_object$metadata[factor_name, ] <- factor(new_factor_values, levels = new_unique_levels)

  # Update steps
  romics_object <- romicsUpdateSteps(romics_object, arguments)

  return(romics_object)
}

#' romicsFactorProportionTestPairwise
#' @description Tests whether the proportions of each level of a given factor differ
#' significantly between pairs of conditions. Uses an omnibus chi-square test on the
#' full contingency table for each pairwise comparison, then derives per-level
#' directionality from standardized Pearson residuals.
#' @param romics_object An romics_object created using romicsCreateObject().
#' @param test_factor A character string indicating the factor whose level proportions
#'   will be tested (e.g., "Leiden_clust_names", "cell_type", "AGER_presence").
#' @param condition_factor A character string indicating the condition factor (e.g., "condition").
#' @param method A character string specifying the test method: "chisq" (chi-square test) or
#'   "fisher" (Fisher's exact test via monte-carlo simulation for tables > 2x2). Default: "chisq".
#' @param simulate_p A logical indicating whether to simulate p-values (recommended when expected
#'   cell frequencies are low, i.e. < 5). Default: TRUE.
#' @param B Integer. Number of replicates for p-value simulation. Default: 10000.
#' @param padj_method Correction method for multiple testing across pairwise comparisons. Default: "BH".
#' @param residual_threshold Minimum absolute standardized Pearson residual required to flag a
#'   level as meaningfully over/under-represented. Conventionally 2.0. Default: 2.0.
#' @param min_effect_size Minimum absolute difference in proportions (%) required as a
#'   descriptive filter. Default: 5.
#' @param min_fold_change Minimum fold change in proportions required as a descriptive filter.
#' @author Geremy Clair
#' @export
romicsFactorProportionTestPairwise <- function(romics_object,
                                         test_factor,
                                         condition_factor   = "condition",
                                         method             = "chisq",
                                         simulate_p         = TRUE,
                                         B                  = 10000,
                                         padj_method        = "BH",
                                         residual_threshold = 2.0,
                                         min_effect_size    = 5,
                                         min_fold_change    = 1.5) {

  # ── Input validation ──────────────────────────────────────────────────────────
  if (!is.romicsObject(romics_object) || missing(romics_object)) {
    stop("<romics_object> is missing or is not in the appropriate format.")
  }
  if (missing(test_factor)) {
    stop("<test_factor> is required.")
  }
  if (!method %in% c("chisq", "fisher")) {
    stop("<method> must be 'chisq' or 'fisher'.")
  }

  # ── Extract metadata ──────────────────────────────────────────────────────────
  metadata <- data.frame(t(romics_object$metadata), stringsAsFactors = TRUE)

  if (!test_factor %in% colnames(metadata)) {
    stop(paste0("test_factor '", test_factor, "' not found in metadata."))
  }
  if (!condition_factor %in% colnames(metadata)) {
    stop(paste0("condition_factor '", condition_factor, "' not found in metadata."))
  }

  all_conditions <- levels(droplevels(metadata[[condition_factor]]))
  if (length(all_conditions) < 2) {
    stop("Need at least 2 condition levels for pairwise comparisons.")
  }

  pairwise_combinations <- combn(all_conditions, 2, simplify = FALSE)
  n_pairs               <- length(pairwise_combinations)

  message("=== Pairwise Proportion Tests (Omnibus) ===")
  message(paste("Test factor      :", test_factor))
  message(paste("Condition factor :", condition_factor))
  message(paste("Conditions       :", paste(all_conditions, collapse = ", ")))
  message(paste("Pairwise comparisons:", n_pairs))
  message(paste("Test method      :", method,
                if (simulate_p) paste0("(simulated p, B=", B, ")") else ""))
  message(paste("P-adj method     :", padj_method, "(across", n_pairs, "comparisons)"))
  message(paste("Residual threshold:", residual_threshold))
  message(paste("Effect size filter: |diff| >=", min_effect_size, "%, FC >=", min_fold_change))

  # ── Full contingency table (all conditions, for reference) ────────────────────
  full_contingency <- table(metadata[[test_factor]], metadata[[condition_factor]])
  full_proportions <- prop.table(full_contingency, margin = 2) * 100

  # ── Initialise result containers ──────────────────────────────────────────────
  results <- list()
  results$full_contingency_table <- full_contingency
  results$full_proportions       <- full_proportions
  results$parameters <- list(
    test_factor        = test_factor,
    condition_factor   = condition_factor,
    method             = method,
    simulate_p         = simulate_p,
    B                  = B,
    padj_method        = padj_method,
    residual_threshold = residual_threshold,
    min_effect_size    = min_effect_size,
    min_fold_change    = min_fold_change
  )
  results$pairwise_comparisons <- list()

  omnibus_rows  <- list()
  per_level_all <- list()

  # ── Main loop: one omnibus test per pair ──────────────────────────────────────
  for (pair_idx in seq_along(pairwise_combinations)) {

    pair  <- pairwise_combinations[[pair_idx]]
    cond1 <- pair[1]
    cond2 <- pair[2]
    comparison_name <- paste(cond1, "vs", cond2, sep = "_")

    message(paste("\n--- Processing:", comparison_name, "---"))

    # Subset to this pair
    pair_rows     <- metadata[[condition_factor]] %in% c(cond1, cond2)
    pair_metadata <- metadata[pair_rows, ]
    pair_metadata[[condition_factor]] <- droplevels(pair_metadata[[condition_factor]])

    # Full K x 2 contingency table; enforce column order
    pair_contingency <- table(pair_metadata[[test_factor]],
                              pair_metadata[[condition_factor]])
    pair_contingency <- pair_contingency[, c(cond1, cond2), drop = FALSE]

    # Column-wise proportions (% within each condition)
    pair_proportions <- prop.table(pair_contingency, margin = 2) * 100

    # ── Omnibus test ────────────────────────────────────────────────────────────
    if (method == "chisq") {
      omnibus_test      <- suppressWarnings(
        chisq.test(pair_contingency, simulate.p.value = simulate_p, B = B)
      )
      omnibus_p         <- omnibus_test$p.value
      omnibus_statistic <- as.numeric(omnibus_test$statistic)
      omnibus_df        <- as.numeric(omnibus_test$parameter)
      std_residuals     <- omnibus_test$stdres

    } else {
      omnibus_test      <- fisher.test(pair_contingency,
                                       simulate.p.value = TRUE, B = B)
      omnibus_p         <- omnibus_test$p.value
      omnibus_statistic <- NA
      omnibus_df        <- NA

      # Compute standardized Pearson residuals manually for Fisher
      expected      <- outer(rowSums(pair_contingency),
                             colSums(pair_contingency)) / sum(pair_contingency)
      row_margins   <- rowSums(pair_contingency) / sum(pair_contingency)
      col_margins   <- colSums(pair_contingency) / sum(pair_contingency)
      std_residuals <- (pair_contingency - expected) /
        sqrt(expected * outer(1 - row_margins, 1 - col_margins))
    }

    message(paste("Omnibus p-value:", format.pval(omnibus_p, digits = 3),
                  "| statistic:", round(omnibus_statistic, 3)))

    # Store omnibus summary row
    omnibus_rows[[comparison_name]] <- data.frame(
      Comparison        = comparison_name,
      Cond1             = cond1,
      Cond2             = cond2,
      N_cond1           = sum(pair_contingency[, cond1]),
      N_cond2           = sum(pair_contingency[, cond2]),
      Omnibus_statistic = omnibus_statistic,
      Omnibus_df        = omnibus_df,
      Omnibus_p         = omnibus_p,
      stringsAsFactors  = FALSE
    )

    # ── Per-level descriptive stats + residuals ─────────────────────────────────
    factor_levels  <- rownames(pair_contingency)
    per_level_rows <- list()

    for (level in factor_levels) {

      prop_c1    <- pair_proportions[level, cond1]
      prop_c2    <- pair_proportions[level, cond2]
      count_c1   <- pair_contingency[level, cond1]
      count_c2   <- pair_contingency[level, cond2]

      abs_diff    <- abs(prop_c1 - prop_c2)
      signed_diff <- prop_c1 - prop_c2   # positive = enriched in cond1

      if (prop_c1 == 0 && prop_c2 == 0) {
        fold_change <- 1
      } else if (prop_c2 == 0) {
        fold_change <- Inf
      } else {
        fold_change <- prop_c1 / prop_c2
      }
      fold_change_abs <- max(fold_change, 1 / fold_change, na.rm = TRUE)

      resid_cond1 <- std_residuals[level, cond1]
      resid_cond2 <- std_residuals[level, cond2]

      # Build row with generic placeholders then rename to actual condition names
      row <- data.frame(
        Comparison      = comparison_name,
        Level           = level,
        N_cond1         = count_c1,
        N_cond2         = count_c2,
        Prop_cond1      = prop_c1,
        Prop_cond2      = prop_c2,
        Abs_Diff        = abs_diff,
        Signed_Diff     = signed_diff,
        Fold_Change     = fold_change,
        Fold_Change_Abs = fold_change_abs,
        Residual_cond1  = resid_cond1,
        Residual_cond2  = resid_cond2,
        Omnibus_p       = omnibus_p,
        stringsAsFactors = FALSE
      )

      colnames(row)[colnames(row) == "N_cond1"]        <- paste0("N_",        cond1)
      colnames(row)[colnames(row) == "N_cond2"]        <- paste0("N_",        cond2)
      colnames(row)[colnames(row) == "Prop_cond1"]     <- paste0("Prop_",     cond1)
      colnames(row)[colnames(row) == "Prop_cond2"]     <- paste0("Prop_",     cond2)
      colnames(row)[colnames(row) == "Residual_cond1"] <- paste0("Residual_", cond1)
      colnames(row)[colnames(row) == "Residual_cond2"] <- paste0("Residual_", cond2)

      per_level_rows[[level]] <- row
    }

    # Within-comparison rbind is safe — all rows share identical columns
    comparison_df            <- do.call(rbind, per_level_rows)
    rownames(comparison_df)  <- NULL

    results$pairwise_comparisons[[comparison_name]] <- list(
      contingency_table = pair_contingency,
      proportions       = pair_proportions,
      std_residuals     = std_residuals,
      omnibus_test      = omnibus_test,
      per_level         = comparison_df
    )

    per_level_all[[comparison_name]] <- comparison_df
  }

  # ── Adjust omnibus p-values across all pairwise comparisons ──────────────────
  # omnibus_rows all have identical columns — safe to do.call(rbind)
  omnibus_df_all              <- do.call(rbind, omnibus_rows)
  rownames(omnibus_df_all)    <- NULL
  omnibus_df_all$Omnibus_padj <- p.adjust(omnibus_df_all$Omnibus_p, method = padj_method)
  results$omnibus_results     <- omnibus_df_all

  # ── Combine per-level results across comparisons ──────────────────────────────
  # Comparisons have different condition-named columns -> use bind_rows (fills NAs)
  all_results           <- dplyr::bind_rows(per_level_all)
  rownames(all_results) <- NULL

  # Merge adjusted omnibus p-value
  all_results <- merge(
    all_results,
    omnibus_df_all[, c("Comparison", "Omnibus_padj")],
    by    = "Comparison",
    all.x = TRUE
  )

  # ── Extract per-row cond1 residual (used for direction and residual filter) ───
  # After bind_rows, condition-specific Residual_ columns exist but are NA for
  # rows from other comparisons. Extract the correct residual per row using the
  # cond1 name stored in omnibus_df_all.
  all_results$Residual_cond1_val <- mapply(
    function(comp, lvl) {
      cond1_name <- omnibus_df_all$Cond1[omnibus_df_all$Comparison == comp]
      col_name   <- paste0("Residual_", cond1_name)
      if (length(col_name) == 1 && col_name %in% colnames(all_results)) {
        val <- all_results[[col_name]][all_results$Comparison == comp &
                                         all_results$Level      == lvl]
        if (length(val) == 1) val else NA_real_
      } else {
        NA_real_
      }
    },
    all_results$Comparison,
    all_results$Level
  )
  all_results$Residual_cond1_val <- as.numeric(all_results$Residual_cond1_val)

  # ── Apply significance filters ────────────────────────────────────────────────
  all_results$Passes_Omnibus     <- all_results$Omnibus_padj        <  0.05
  all_results$Passes_Residual    <- abs(all_results$Residual_cond1_val) >= residual_threshold
  all_results$Passes_Effect_Size <- all_results$Abs_Diff             >= min_effect_size
  all_results$Passes_Fold_Change <- all_results$Fold_Change_Abs      >= min_fold_change
  all_results$Significant_Final  <- (all_results$Passes_Omnibus      &
                                       all_results$Passes_Residual     &
                                       all_results$Passes_Effect_Size  &
                                       all_results$Passes_Fold_Change)

  all_results$Direction <- ifelse(
    !all_results$Significant_Final,
    "Not Significant",
    ifelse(
      all_results$Residual_cond1_val > 0,
      paste0("Enriched_in_", sub("_vs_.*", "", all_results$Comparison)),
      paste0("Depleted_in_", sub("_vs_.*", "", all_results$Comparison))
    )
  )

  results$all_pairwise_results <- all_results

  # ── Print summary ─────────────────────────────────────────────────────────────
  message("\n=== Omnibus Test Summary ===")
  message(sprintf("%-25s  %8s  %8s  %s", "Comparison", "p_omnibus", "p_adj", "Sig"))
  message(strrep("-", 55))
  for (i in seq_len(nrow(omnibus_df_all))) {
    message(sprintf("%-25s  %8s  %8s  %s",
                    omnibus_df_all$Comparison[i],
                    format.pval(omnibus_df_all$Omnibus_p[i],    digits = 3),
                    format.pval(omnibus_df_all$Omnibus_padj[i], digits = 3),
                    if (omnibus_df_all$Omnibus_padj[i] < 0.05) "***" else "ns"))
  }

  message("\n=== Per-Level Results (Significant Final) ===")
  sig <- all_results[all_results$Significant_Final, ]
  if (nrow(sig) > 0) {
    sig_sorted <- sig[order(sig$Comparison, sig$Omnibus_padj), ]
    for (i in seq_len(nrow(sig_sorted))) {
      message(sprintf("  [%s] %s | diff=%.1f%% | FC=%.2f | residual=%.2f | padj=%s | %s",
                      sig_sorted$Comparison[i],
                      sig_sorted$Level[i],
                      sig_sorted$Signed_Diff[i],
                      sig_sorted$Fold_Change[i],
                      sig_sorted$Residual_cond1_val[i],
                      format.pval(sig_sorted$Omnibus_padj[i], digits = 3),
                      sig_sorted$Direction[i]))
    }
  } else {
    message("  No levels passed all filters.")
  }

  message(paste0("\nTotal levels tested       : ", nrow(all_results)))
  message(paste0("Passing omnibus p_adj<0.05: ", sum(all_results$Passes_Omnibus,     na.rm = TRUE)))
  message(paste0("Passing residual threshold: ", sum(all_results$Passes_Residual,    na.rm = TRUE)))
  message(paste0("Passing effect size filter: ", sum(all_results$Passes_Effect_Size, na.rm = TRUE)))
  message(paste0("Passing fold change filter: ", sum(all_results$Passes_Fold_Change, na.rm = TRUE)))
  message(paste0("Significant (all filters) : ", sum(all_results$Significant_Final,  na.rm = TRUE)))

  # ── Summary table (long format: one row per level x comparison) ──────────────
  # Extract per-row Prop_cond1 and Prop_cond2 using condition names from omnibus table
  prop_cond1_vals <- mapply(function(comp, lvl) {
    cond1_name <- omnibus_df_all$Cond1[omnibus_df_all$Comparison == comp]
    col_name   <- paste0("Prop_", cond1_name)
    if (length(col_name) == 1 && col_name %in% colnames(all_results)) {
      val <- all_results[[col_name]][all_results$Comparison == comp &
                                       all_results$Level      == lvl]
      if (length(val) == 1) round(val, 3) else NA_real_
    } else NA_real_
  }, all_results$Comparison, all_results$Level)

  prop_cond2_vals <- mapply(function(comp, lvl) {
    cond2_name <- omnibus_df_all$Cond2[omnibus_df_all$Comparison == comp]
    col_name   <- paste0("Prop_", cond2_name)
    if (length(col_name) == 1 && col_name %in% colnames(all_results)) {
      val <- all_results[[col_name]][all_results$Comparison == comp &
                                       all_results$Level      == lvl]
      if (length(val) == 1) round(val, 3) else NA_real_
    } else NA_real_
  }, all_results$Comparison, all_results$Level)

  # Build summary table with labelled proportion columns
  summary_table <- data.frame(
    Comparison        = all_results$Comparison,
    Level             = all_results$Level,
    Prop_cond1        = as.numeric(prop_cond1_vals),
    Prop_cond2        = as.numeric(prop_cond2_vals),
    Signed_Diff       = round(all_results$Signed_Diff,        3),
    Fold_Change       = round(all_results$Fold_Change,        3),
    Residual          = round(all_results$Residual_cond1_val, 3),
    Omnibus_p         = all_results$Omnibus_p,
    Omnibus_padj      = all_results$Omnibus_padj,
    Direction         = all_results$Direction,
    Significant_Final = all_results$Significant_Final,
    stringsAsFactors  = FALSE
  )

  # Label Prop columns with actual condition names: "Prop_CTL (cond1)", etc.
  # Merge cond1/cond2 labels per comparison for readable column naming
  # Since different rows have different cond1/cond2, we label generically but
  # add a Cond1 and Cond2 column so the reader knows which condition is which
  cond_labels <- omnibus_df_all[, c("Comparison", "Cond1", "Cond2")]
  summary_table <- merge(summary_table, cond_labels, by = "Comparison", all.x = TRUE)

  # Reorder columns for readability
  summary_table <- summary_table[, c("Comparison", "Cond1", "Cond2", "Level",
                                     "Prop_cond1", "Prop_cond2",
                                     "Signed_Diff", "Fold_Change", "Residual",
                                     "Omnibus_p", "Omnibus_padj",
                                     "Direction", "Significant_Final")]

  # Sort: significant first, then by comparison and descending |residual|
  summary_table <- summary_table[order(!summary_table$Significant_Final,
                                       summary_table$Comparison,
                                       -abs(summary_table$Residual)), ]
  rownames(summary_table) <- NULL

  # Print summary table
  message("\n=== Summary Table ===")
  print(summary_table, row.names = FALSE)

  results$summary_table <- summary_table

  return(results)
}

#' romicsProportionHeatmap
#' @description Plots a heatmap of signed proportion differences (%) from the
#'   output of romicsProportionTestPairwise(). Rows are factor levels
#'   ordered alphabetically, columns are pairwise comparisons. Fill colour
#'   encodes Signed_Diff (proportion in cond1 minus proportion in cond2).
#'   Each tile displays the rounded signed difference followed by a significance level.
#' @param proportion_test_result The list returned by romicsProportionTestPairwise.
#' @param title A character string for the plot title. Default:
#'   "Proportion Difference Heatmap".
#' @param midpoint Numeric. Centre of the diverging colour scale. Default: 0.
#' @param low_color Colour for negative differences (depleted in cond1).
#'   Default: "blue".
#' @param mid_color Colour for zero difference. Default: "white".
#' @param high_color Colour for positive differences (enriched in cond1).
#'   Default: "red".
#' @param label_size Numeric. Font size of the tile labels. Default: 3.
#' @param text_size Numeric. Base font size for the plot. Default: 11.
#' @param interactive Logical. If TRUE returns a plotly object; if FALSE
#'   returns a ggplot2 object. Default: FALSE.
#' @return A ggplot2 object, or a plotly object if interactive = TRUE.
#' @author Geremy Clair
#' @export
romicsProportionHeatmap <- function(proportion_test_result,
                                    title      = "Proportion Difference Heatmap",
                                    midpoint   = 0,
                                    low_color  = "blue",
                                    mid_color  = "white",
                                    high_color = "red",
                                    label_size = 3,
                                    text_size  = 11,
                                    interactive = FALSE) {

  # ── Input validation ──────────────────────────────────────────────────────────
  if (!is.list(proportion_test_result) ||
      !"summary_table" %in% names(proportion_test_result)) {
    stop("proportion_test_result must be the list returned by romicsProportionTestPairwise().")
  }

  summary_table <- proportion_test_result$summary_table

  required_cols <- c("Comparison", "Level", "Signed_Diff",
                     "Significant_Final", "Omnibus_padj")
  missing_cols  <- setdiff(required_cols, colnames(summary_table))
  if (length(missing_cols) > 0) {
    stop(paste("summary_table is missing columns:", paste(missing_cols, collapse = ", ")))
  }

  # ── Prepare plot data ─────────────────────────────────────────────────────────
  plot_data <- summary_table[, c("Comparison", "Level", "Signed_Diff",
                                 "Significant_Final", "Omnibus_padj",
                                 "Residual", "Prop_cond1", "Prop_cond2",
                                 "Cond1", "Cond2")]

  # Order levels alphabetically (decreasing so A appears at top of y-axis)
  level_order          <- sort(unique(plot_data$Level), decreasing = TRUE)
  plot_data$Level      <- factor(plot_data$Level,      levels = level_order)
  plot_data$Comparison <- factor(plot_data$Comparison)

  # Tile label: rounded signed diff + significance stars in parentheses
  # Stars only shown when the level passes all filters (Significant_Final == TRUE)
  plot_data$tile_label <- paste0(
    round(plot_data$Signed_Diff, 1),
    dplyr::case_when(
      !plot_data$Significant_Final   ~ "",
      plot_data$Omnibus_padj < 0.001 ~ " (***)",
      plot_data$Omnibus_padj < 0.01  ~ " (**)",
      plot_data$Omnibus_padj < 0.05  ~ " (*)",
      TRUE                           ~ ""
    )
  )

  # Tooltip text for plotly
  plot_data$tooltip <- paste0(
    "Comparison: ",  plot_data$Comparison,                    "\n",
    "Level: ",       plot_data$Level,                         "\n",
    "Cond1 (", plot_data$Cond1, "): ",
    round(plot_data$Prop_cond1, 2), "%\n",
    "Cond2 (", plot_data$Cond2, "): ",
    round(plot_data$Prop_cond2, 2), "%\n",
    "Signed diff: ", round(plot_data$Signed_Diff, 2),   "%\n",
    "Residual: ",    round(plot_data$Residual,    2),    "\n",
    "padj: ",        signif(plot_data$Omnibus_padj, 3),  "\n",
    "Significant: ", plot_data$Significant_Final
  )

  # Symmetric colour scale limit
  max_abs <- max(abs(plot_data$Signed_Diff), na.rm = TRUE)

  # ── Build ggplot ──────────────────────────────────────────────────────────────
  p <- ggplot2::ggplot(plot_data,
                       ggplot2::aes(x    = Comparison,
                                    y    = Level,
                                    fill = Signed_Diff,
                                    text = tooltip)) +
    ggplot2::geom_tile(color = "white", linewidth = 1) +
    ggplot2::geom_text(ggplot2::aes(label = tile_label),
                       color = "black",
                       size  = label_size) +
    ggplot2::scale_fill_gradient2(
      low      = low_color,
      mid      = mid_color,
      high     = high_color,
      midpoint = midpoint,
      limits   = c(-max_abs, max_abs),
      name     = "Signed Diff\n(cond1 - cond2, %)"
    ) +
    ggplot2::labs(
      title    = title,
      subtitle = "Red = enriched in cond1 | Blue = enriched in cond2 | (*) p<0.05  (**) p<0.01  (***) p<0.001",
      x        = "Comparison",
      y        = "Level"
    ) +
    ggplot2::theme_minimal(base_size = text_size) +
    ggplot2::theme(
      axis.text.x   = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
      plot.title    = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 9)
    )

  # ── Return ────────────────────────────────────────────────────────────────────
  if (interactive) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      stop("Package 'plotly' is required for interactive = TRUE.")
    }
    return(plotly::ggplotly(p, tooltip = "text"))
  }

  return(p)
}

#' romicsImportFactor()
#' @description Import factors from external objects (character vector, data frame, or transposed data frame)
#' into a romics_object's metadata layer. Supports partial matching with "not_defined" for missing values.
#' @param romics_object A romics_object created using romicsCreateObject().
#' @param data_source A character vector, data frame, or transposed data frame containing factor values.
#'        - Character vector: names must match colnames(romics_object$data)
#'        - Data frame: colnames must match colnames(romics_object$data), rownames specify factor names
#'        - Transposed data frame: rownames must match colnames(romics_object$data), colnames are factor names
#' @param factor_name Character vector of factor names. Used when data_source is a character vector or
#'        when rownames of data frame are numbered/non-informative. If NULL, extracts from rownames. Default: NULL
#' @param allow_partial Logical. If TRUE, allows partial matching of samples with "not_defined" for missing values.
#'        If FALSE, throws error if not all samples are covered. Default: TRUE
#' @param overwrite Logical. If TRUE, allows overwriting existing factors with the same name.
#'        If FALSE, throws error if factor already exists. Default: FALSE
#' @param verbose Logical. Print messages about import progress. Default: TRUE
#' @return A romics_object with new factors added to the metadata layer.
#' @details This function is flexible and handles multiple input formats:
#' - Single factor as named character vector: names are sample IDs, values are factor levels
#' - Multiple factors as data frame: columns are samples, rows are factors
#' - Multiple factors as transposed data frame: rows are samples, columns are factors
#'
#' Partial matching allows importing data with fewer samples than the romics_object. Unmatched
#' samples receive "not_defined" values and are reported in the console.
#' @examples \dontrun{
#' # Import from named character vector
#' factor_vector <- c(sample1 = "groupA", sample2 = "groupB", sample3 = "groupA")
#' romics_obj <- romicsImportFactor(romics_obj, factor_vector, factor_name = "my_factor")
#'
#' # Import multiple factors from data frame (samples as columns)
#' factor_df <- data.frame(sample1 = c("groupA", "typeX"),
#'                         sample2 = c("groupB", "typeY"),
#'                         sample3 = c("groupA", "typeX"),
#'                         row.names = c("group", "type"))
#' romics_obj <- romicsImportFactor(romics_obj, factor_df)
#'
#' # Import from transposed data frame (samples as rows)
#' factor_df_t <- data.frame(group = c("groupA", "groupB", "groupA"),
#'                           type = c("typeX", "typeY", "typeX"),
#'                           row.names = c("sample1", "sample2", "sample3"))
#' romics_obj <- romicsImportFactor(romics_obj, factor_df_t)
#' }
#' @author Geremy Clair
#' @export
romicsImportFactor <- function(romics_object,
                               data_source,
                               factor_name = NULL,
                               allow_partial = TRUE,
                               overwrite = FALSE,
                               verbose = TRUE) {

  arguments <- as.list(match.call())

  # Input validation
  if (!is.romicsObject(romics_object) || missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }

  if (missing(data_source)) {
    stop("data_source is required")
  }

  sample_names <- colnames(romics_object$data)
  n_samples <- length(sample_names)

  # ── Handle character vector input ──────────────────────────────────────────
  if (is.character(data_source) && !is.data.frame(data_source)) {
    if (is.null(names(data_source))) {
      stop("Character vector data_source must have names corresponding to sample names")
    }

    if (is.null(factor_name)) {
      factor_name <- "imported_factor"
    }

    if (length(factor_name) != 1) {
      stop("When data_source is a character vector, factor_name must be a single string")
    }

    # Match samples
    matched_samples <- match(sample_names, names(data_source))
    new_factor <- rep("not_defined", n_samples)
    matched_idx <- !is.na(matched_samples)
    new_factor[matched_idx] <- data_source[matched_samples[matched_idx]]

    n_matched <- sum(!is.na(matched_samples))
    n_missing <- n_samples - n_matched

    if (n_missing > 0) {
      if (!allow_partial) {
        stop("Only ", n_matched, " out of ", n_samples, " samples were found in data_source.\n",
             "Set allow_partial = TRUE to allow partial matching with 'not_defined' for missing samples.")
      }
      if (verbose) {
        message("Imported factor '", factor_name, "': ", n_matched, " matched, ", n_missing, " set to 'not_defined'")
      }
    } else {
      if (verbose) {
        message("Imported factor '", factor_name, "': all ", n_samples, " samples matched")
      }
    }

    # Check if factor already exists
    if (factor_name %in% rownames(romics_object$metadata)) {
      if (!overwrite) {
        stop("Factor '", factor_name, "' already exists in metadata. Set overwrite = TRUE to replace it.")
      }
      if (verbose) {
        message("Overwriting existing factor '", factor_name, "'")
      }
      # Remove the existing factor and add the new one
      romics_object$metadata <- romics_object$metadata[rownames(romics_object$metadata) != factor_name, , drop = FALSE]
    }

    romics_object$metadata <- rbind(romics_object$metadata, factor_name = new_factor)

  # ── Handle data frame input ────────────────────────────────────────────────
  } else if (is.data.frame(data_source)) {
    # Determine if data frame is transposed (rows=samples) or normal (columns=samples)

    # Check if colnames match sample names (normal orientation: columns=samples)
    colnames_match <- all(colnames(data_source) %in% sample_names)
    rownames_match <- all(rownames(data_source) %in% sample_names)

    if (colnames_match && !rownames_match) {
      # Normal orientation: columns are samples, rows are factors
      if (verbose) {
        message("Detected normal orientation: columns=samples, rows=factors")
      }

      # Extract factor names from rownames
      df_factor_names <- rownames(data_source)

      # If factor_name is provided, use it to rename
      if (!is.null(factor_name)) {
        if (length(factor_name) != nrow(data_source)) {
          stop("factor_name length (", length(factor_name), ") must match data_source rows (",
               nrow(data_source), ")")
        }
        df_factor_names <- factor_name
      }

      # Match samples and import each factor
      for (i in seq_along(df_factor_names)) {
        factor_name_i <- df_factor_names[i]
        factor_values_i <- as.character(data_source[i, ])

        # Reorder to match romics_object samples
        matched_samples <- match(sample_names, colnames(data_source))
        new_factor <- rep("not_defined", n_samples)
        matched_idx <- !is.na(matched_samples)
        new_factor[matched_idx] <- factor_values_i[matched_samples[matched_idx]]

        n_matched <- sum(!is.na(matched_samples))
        n_missing <- n_samples - n_matched

        if (n_missing > 0) {
          if (!allow_partial) {
            stop("Only ", n_matched, " out of ", n_samples, " samples were found for factor '",
                 factor_name_i, "'")
          }
          if (verbose) {
            message("  - Factor '", factor_name_i, "': ", n_matched, " matched, ", n_missing, " set to 'not_defined'")
          }
        } else {
          if (verbose) {
            message("  - Factor '", factor_name_i, "': all ", n_samples, " samples matched")
          }
        }

        # Check if factor already exists
        if (factor_name_i %in% rownames(romics_object$metadata)) {
          if (!overwrite) {
            stop("Factor '", factor_name_i, "' already exists in metadata. Set overwrite = TRUE to replace it.")
          }
          if (verbose) {
            message("    Overwriting existing factor '", factor_name_i, "'")
          }
          # Remove the existing factor and add the new one
          romics_object$metadata <- romics_object$metadata[rownames(romics_object$metadata) != factor_name_i, , drop = FALSE]
        }

        romics_object$metadata <- rbind(romics_object$metadata, new_factor)
        rownames(romics_object$metadata)[nrow(romics_object$metadata)] <- factor_name_i
      }

    } else if (rownames_match && !colnames_match) {
      # Transposed orientation: rows are samples, columns are factors
      if (verbose) {
        message("Detected transposed orientation: rows=samples, columns=factors")
      }

      # Extract factor names from colnames
      df_factor_names <- colnames(data_source)

      # If factor_name is provided, use it to rename
      if (!is.null(factor_name)) {
        if (length(factor_name) != ncol(data_source)) {
          stop("factor_name length (", length(factor_name), ") must match data_source columns (",
               ncol(data_source), ")")
        }
        df_factor_names <- factor_name
      }

      # Match samples and import each factor
      for (i in seq_along(df_factor_names)) {
        factor_name_i <- df_factor_names[i]

        # Reorder to match romics_object samples
        matched_samples <- match(sample_names, rownames(data_source))
        new_factor <- rep("not_defined", n_samples)
        matched_idx <- !is.na(matched_samples)
        new_factor[matched_idx] <- as.character(data_source[matched_samples[matched_idx], i])

        n_matched <- sum(!is.na(matched_samples))
        n_missing <- n_samples - n_matched

        if (n_missing > 0) {
          if (!allow_partial) {
            stop("Only ", n_matched, " out of ", n_samples, " samples were found for factor '",
                 factor_name_i, "'")
          }
          if (verbose) {
            message("  - Factor '", factor_name_i, "': ", n_matched, " matched, ", n_missing, " set to 'not_defined'")
          }
        } else {
          if (verbose) {
            message("  - Factor '", factor_name_i, "': all ", n_samples, " samples matched")
          }
        }

        # Check if factor already exists
        if (factor_name_i %in% rownames(romics_object$metadata)) {
          if (!overwrite) {
            stop("Factor '", factor_name_i, "' already exists in metadata. Set overwrite = TRUE to replace it.")
          }
          if (verbose) {
            message("    Overwriting existing factor '", factor_name_i, "'")
          }
          # Remove the existing factor and add the new one
          romics_object$metadata <- romics_object$metadata[rownames(romics_object$metadata) != factor_name_i, , drop = FALSE]
        }

        romics_object$metadata <- rbind(romics_object$metadata, new_factor)
        rownames(romics_object$metadata)[nrow(romics_object$metadata)] <- factor_name_i
      }

    } else {
      stop("Could not determine data frame orientation.\n",
           "Either colnames or rownames must match the romics_object sample names (",
           paste(sample_names[1:min(3, length(sample_names))], collapse = ", "), ", ...)")
    }

  } else {
    stop("data_source must be a character vector or data frame")
  }

  # Update processing steps
  romics_object <- romicsUpdateSteps(romics_object, arguments)

  return(romics_object)
}
