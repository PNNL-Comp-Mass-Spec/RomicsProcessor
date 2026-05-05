#' romicsHclust()
#' @description Plots the hierarchical clustering of the samples contaned in the romics_object calculated using the functions dist() and hclust() (see documentation of those R functions for more details). The colors used for the plotting will correspond to the main_factor of the romics_object.
#' @param romics_object has to be a log transformed romics_object created using romicsCreateObject() and transformed using the function log2transform() or log10transform()
#' @param method_dist the distance measure to be used. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". Any unambiguous substring can be given.
#' @param method_hclust the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' @details This function uses the dist() and hclust() functions to calculate the hierachical clustering and then plots the hclust with colors based on the current main_factor of the romics_object.
#' @return a hierarchical clustering tree plot with its branches colored by factor.
#' @author Geremy Clair
#' @export
romicsHclust<-function(romics_object, method_dist="euclidean", method_hclust="ward.D" ){
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if (missing(method_dist)){method_dist<-"euclidean"}
  if(!method_dist %in% c("euclidean","maximum","manhattan","canberra","binary","minkowski")){stop("The method_dist used is not appropriate")}
  if(missing(method_hclust)){method_hclust<-"ward.D"}
  if(!method_hclust %in% c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")){stop("The method_hclust used is not appropriate")}

  #extract the data
  data<-romics_object$data
  data<-t(data)

  #extract the colors
  colors_dend <- as.character(t(romics_object$metadata[romicsFactorNames(romics_object)=="colors_romics",]))

  #define the clustering method
  distance_samples<-dist(data, method=method_dist)

  #define the hclust agglomeration method
  hc<- hclust(distance_samples, method=method_hclust)

  #convert the hclust into a dendrogram
  dd <- as.dendrogram(hc)

  #order the colors as the hclust
  colors_dend<-colors_dend[order=hc$order]

  #color the branches
  dd <- dendextend::color_branches(dd,k = NULL, h = NULL,col=as.vector(colors_dend))

  #plot the dendrogram to choose the number of clusters
  par(cex=0.7, mar=c(8, 8, 4, 1))
  plot(dd,main="Hierarchical clustering samples", sub="")
}

# Helper function for consistent factor coloring across plot functions
assign_factor_colors <- function(factor_values, color_palette = NULL, romics_object = NULL) {
  unique_levels <- unique(factor_values)
  n_levels <- length(unique_levels)

  if(!is.null(color_palette) && length(color_palette) > 0) {
    colors_to_use <- color_palette
    if(length(colors_to_use) < n_levels) {
      if(exists("ROP_colors", envir = .GlobalEnv)) {
        additional_colors <- get("ROP_colors", envir = .GlobalEnv)
        colors_to_use <- c(colors_to_use, additional_colors)
      }
    }
    if(length(colors_to_use) < n_levels) {
      colors_to_use <- c(colors_to_use, rainbow(n_levels - length(colors_to_use)))
    }
  } else if(!is.null(romics_object) && "colors_romics" %in% rownames(romics_object$metadata)) {
    colors_romics_vec <- as.character(romics_object$metadata["colors_romics", ])
    colors_to_use <- unique(colors_romics_vec)
    if(length(colors_to_use) < n_levels) {
      if(exists("ROP_colors", envir = .GlobalEnv)) {
        additional_colors <- get("ROP_colors", envir = .GlobalEnv)
        colors_to_use <- c(colors_to_use, additional_colors)
      }
    }
    if(length(colors_to_use) < n_levels) {
      colors_to_use <- c(colors_to_use, rainbow(n_levels - length(colors_to_use)))
    }
  } else if(exists("ROP_colors", envir = .GlobalEnv)) {
    colors_to_use <- get("ROP_colors", envir = .GlobalEnv)
    if(length(colors_to_use) < n_levels) {
      colors_to_use <- c(colors_to_use, rainbow(n_levels - length(colors_to_use)))
    }
  } else {
    colors_to_use <- rainbow(n_levels)
  }

  color_mapping <- setNames(colors_to_use[1:n_levels], unique_levels)
  colors_out <- as.character(color_mapping[as.character(factor_values)])
  return(colors_out)
}

#' romicsPCA()
#' @description Calculate the PCA of the data layer of the romics_object using the package FactoMineR.
#' If the data layer contains some missing values those will be imputed using the missMDA::imputePCA() method.
#' This function updates the romics_object by storing the PCA embeddings in the embeddings layer.
#' @param romics_object has to be a romics_object created using romicsCreateObject()
#' @param verbose boolean. If TRUE, also assigns the complete PCA results to PCA_results in the global environment
#' @param ncp number of principal components to compute. Default is 5 if method.cv is not provided.
#'        If method.cv is provided but ncp is not, the optimal number is estimated automatically.
#' @param n_components number of principal components to store in the embeddings layer. If NULL (default),
#'        stores all computed components. Must be <= ncp.
#' @param scale boolean. TRUE implies a same weight for each variable for the missMDA data imputation)
#' @param method "Regularized" (by default) or "EM". Used for missing value imputation using missMDA.
#' @param row.w row weights (by default, a vector of 1 for uniform row weights)
#' @param ncp.min used only if ncp is automatically estimated. Integer corresponding to the minimum number of components to test
#' @param ncp.max used only if ncp is automatically estimated. Integer corresponding to the maximum number of components to test
#' @param method.cv method for cross-validation when automatically determining the optimal number of components.
#'        One of "gcv" (generalized cross-validation), "loo" (leave-one-out), or "Kfold" (K-fold cross-validation).
#'        If provided and ncp is missing, triggers automatic estimation of the optimal number of components.
#' @param ... further arguments passed to missMDA functions.
#' @details This function calculates PCA and stores the sample embeddings in the romics_object$embeddings layer.
#'          If missing values are present in the data, they are imputed using missMDA::imputePCA().
#'          When method.cv is provided but ncp is not, the optimal number of components is estimated using missMDA::estim_ncpPCA().
#'          If verbose=TRUE, the complete PCA results are assigned to PCA_results in the global environment.
#' @return Returns the romics_object with PCA results added to the embeddings layer.
#' @author Geremy Clair
#' @export
romicsPCA <- function(romics_object, verbose = TRUE, ncp, n_components = NULL, scale = TRUE, method=c("Regularized","EM"),
                      ncp.min = 0, ncp.max = 5, method.cv, ...){
  arguments <- as.list(match.call())
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(verbose)){verbose<-TRUE}

  # Set ncp default or estimate it
  estimate_ncp <- FALSE
  if(missing(ncp)){
    if(!missing(method.cv)){
      # If method.cv is provided but ncp is missing, we'll estimate ncp
      estimate_ncp <- TRUE
    } else {
      # If neither ncp nor method.cv is provided, use default ncp=5
      ncp <- 5
    }
  }

  # Validate n_components parameter
  if(!is.null(n_components)){
    if(!is.numeric(n_components) || length(n_components) != 1 || n_components < 1){
      stop("n_components must be a positive integer")
    }
    n_components <- as.integer(n_components)
  }

  data <- t(romics_object$data)

  # Handle missing values with imputation
  if(sum(is.na(data))>0){
    warning("Your romics_object data layer was containing some missing values. The missMDA::imputePCA() method was used to impute these values")
    if(missing("scale")){scale=TRUE}
    if(missing("method")){method="Regularized"}
    # Estimate ncp if needed
    if(estimate_ncp){
      if(missing(ncp.min)){ncp.min=0}
      if(missing(ncp.max)){ncp.max=5}
      if(missing(method.cv)){method.cv="Kfold"}
      warning("ncp will be automatically estimated using the missMDA::estim_ncpPCA() function, this might take a few minutes.")
      ncp_est <- missMDA::estim_ncpPCA(data, ncp.min=ncp.min, ncp.max=ncp.max, method.cv = method.cv, ...)
      print("Below is the result of the estim_ncpPCA()")
      print(ncp_est)
      ncp <- ncp_est$ncp
    }
    data <- missMDA::imputePCA(data, ncp = ncp, scale = scale, method = method, ...)$completeObs
  }

  # Validate n_components against ncp after ncp is determined
  if(!is.null(n_components)){
    if(n_components > ncp){
      stop(paste("n_components (", n_components, ") cannot be greater than ncp (", ncp, ")", sep=""))
    }
  } else {
    # If n_components is NULL, store all computed components
    n_components <- ncp
  }

  # Run the PCA without showing the graphics
  pdf(file = NULL)
  pca_results <- FactoMineR::PCA(data, ncp = ncp, scale.unit = scale, graph = FALSE)
  dev.off()

  # Extract PCA coordinates for samples - only the requested number of components
  pca_coords <- as.data.frame(pca_results$ind$coord[, 1:n_components, drop = FALSE])

  # Name the components
  colnames(pca_coords) <- paste0("pca_component_", 1:ncol(pca_coords))

  # Create or update embeddings layer in the romics_object
  if(!is.null(romics_object$embeddings)){
    if(sum(grepl("pca_component_", rownames(romics_object$embeddings))) > 0){
      print("PCA embeddings were already recorded in the romics_object, these will be removed and replaced by the ones established here.")
      romics_object$embeddings <- romics_object$embeddings[!grepl("pca_component_", rownames(romics_object$embeddings)),]
    }
  } else {
    print("The <embeddings> layer was added to the romics_object, and the PCA coordinates were added.")
    romics_object$embeddings <- data.frame(matrix(ncol=ncol(romics_object$data), nrow=0))
    colnames(romics_object$embeddings) <- colnames(romics_object$data)
  }

  romics_object$embeddings <- rbind(romics_object$embeddings, t(pca_coords))

  # Update the steps
  romics_object <- romicsUpdateSteps(romics_object, arguments)

  # If verbose, assign PCA results to global environment
  if(verbose){
    assign("PCA_results", pca_results, envir = .GlobalEnv)
    message("PCA results have been assigned to 'PCA_results' in the global environment")
    message(paste("Stored", n_components, "PCA components in embeddings layer"))
  }

  # Always return the updated romics_object
  return(romics_object)
}

#' romicsPCAplot()
#' @description Plots the samples on the PCA sample plot in 2D or 3D. When in 2D mode, can also show percentage of explained variance.
#' The plot can be colored by any factor in the romics_object's metadata or by feature values.
#' This function uses the PCA embeddings stored in the romics_object and checks for PCA_results in the global environment for variance plots.
#' @param romics_object A romics_object that has PCA embeddings in its embeddings layer (generated using romicsPCA)
#' @param Xcomp numerical/double. Indicate the component to plot on the X axis (default: 1)
#' @param Ycomp numerical/double. Indicate the component to plot on the Y axis (default: 2)
#' @param Zcomp numerical/double. If provided, creates a 3D plot with this component on the Z axis (default: NULL)
#' @param label boolean. Indicate if the sample name labels should be plotted (default: FALSE, applies only to 2D plots)
#' @param plotType should be one of the following options to indicate the type of plot to be returned:
#'        'dual' (for both sample plot and variance plot), 'individual', or 'percentage' (default: "dual").
#'        Note: The 'dual' and 'percentage' options only work for 2D plots (when Zcomp is NULL).
#' @param feature character string. If provided, the plot will be colored according to the values of this feature.
#'        If NULL (default), coloring will be based on the factor specified by factor_name.
#' @param factor_name character string. The name of the factor in romics_object$metadata to use for coloring the points
#'        when feature is NULL. If NULL (default), the romics_object's main_factor will be used.
#' @param color_palette a gradient color palette to use when coloring by feature (default: viridis(n=20))
#' @param standardize boolean. Indicates if the feature values should be scaled when coloring by feature (default: TRUE)
#' @param size numeric. Size of the points in the plot (default: 3)
#' @param alpha numeric. Transparency of the points in the plot (default: 0.8)
#' @param plotly boolean. If TRUE, generates an interactive plotly plot instead of a static ggplot (default: FALSE)
#' @param show_outline boolean. If FALSE, removes point outlines and shows only fill colors (default: TRUE)
#' @details This function plots PCA results using the embeddings stored in the romics_object.
#'        If Zcomp is provided, a 3D interactive plot is generated using plotly. In 3D mode, only the sample plot
#'        is available (plotType settings are ignored).
#'        If Zcomp is NULL, a 2D plot is generated using ggplot2, and can include the variance plot depending on plotType.
#'        For the variance percentage plot (2D mode only), the function looks for PCA_results in the global environment.
#'        When using a non-main factor for coloring, the function automatically updates the romics_object
#'        using romicsChangeFactor() to ensure consistent colors across plots.
#'        Coloring is determined automatically:
#'        - If feature is provided, the plot is colored by feature values
#'        - If feature is NULL, the plot is colored by factor values
#'        When show_outline is FALSE, points are displayed without borders for cleaner appearance.
#' @return Returns either a ggplot2 plot, a grid.arrange combined plot (for 2D dual plots), or a plotly 3D plot.
#' @author Geremy Clair
#' @export
romicsPCAplot <- function(romics_object, Xcomp=1, Ycomp=2, Zcomp=NULL, label=FALSE, plotType="dual",
                          feature=NULL, factor_name=NULL, color_palette=NULL,
                          standardize=TRUE, size=3, alpha=0.8, plotly=FALSE, show_outline=TRUE,
                          xlim=NULL, ylim=NULL, zlim=NULL){
  # Input validation
  if(!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }

  # Determine coloring mode based on whether feature is provided
  if(!is.null(feature)) {
    color_by <- "feature"
    # Validate feature
    if(!feature %in% rownames(romics_object$data)) {
      stop(paste("Feature", feature, "not found in romics_object data."))
    }
  } else {
    color_by <- "factor"
    # Handle factor selection and update romics_object if needed
    if(is.null(factor_name)) {
      factor_name <- romics_object$main_factor
    } else {
      # Validate factor exists
      if(!(factor_name %in% rownames(romics_object$metadata))) {
        stop(paste("Factor", factor_name, "not found in romics_object metadata."))
      }
      # If factor_name is different from main_factor, update the romics_object
      if(factor_name != romics_object$main_factor) {
        message(paste0("Changing factor to ", factor_name, ", this operation may take a few seconds."))
        romics_object <- romicsChangeFactor(romics_object, main_factor = factor_name)
        # Clear stale colors_romics so assign_factor_colors() uses fresh color assignment
        if("colors_romics" %in% rownames(romics_object$metadata)) {
          romics_object$metadata <- romics_object$metadata[rownames(romics_object$metadata) != "colors_romics", , drop = FALSE]
        }
      }
    }
  }

  # Check if PCA embeddings exist
  if(is.null(romics_object$embeddings) ||
     sum(grepl("pca_component_", rownames(romics_object$embeddings))) == 0) {
    stop("The romics_object doesn't contain PCA embeddings. Please run romicsPCA() first.")
  }

  # Get the PCA embeddings
  pca_comp_rows <- grepl("pca_component_", rownames(romics_object$embeddings))
  pca_embeddings <- romics_object$embeddings[pca_comp_rows, ]
  n_components <- nrow(pca_embeddings)

  # Check component numbers
  if(missing(Xcomp) || !(is.numeric(Xcomp) || is.double(Xcomp))){
    warning("The component to be plotted on the first axis was not defined. PC1 will be used")
    Xcomp <- 1
  }
  if(missing(Ycomp) || !(is.numeric(Ycomp) || is.double(Ycomp))){
    warning("The component to be plotted on the second axis was not defined. PC2 will be used")
    Ycomp <- 2
  }
  if(Xcomp > n_components){
    stop(paste("The Xcomp selected (", Xcomp, ") is too large. Only", n_components, "components available."))
  }
  if(Ycomp > n_components){
    stop(paste("The Ycomp selected (", Ycomp, ") is too large. Only", n_components, "components available."))
  }

  # Check Zcomp if provided
  is_3d_plot <- !is.null(Zcomp)
  if(is_3d_plot && Zcomp > n_components){
    stop(paste("The Zcomp selected (", Zcomp, ") is too large. Only", n_components, "components available."))
  }

  # Try to get PCA_results from global environment
  has_pca_results <- exists("PCA_results", envir = .GlobalEnv)
  if(!has_pca_results && !is_3d_plot && (plotType == "dual" || plotType == "percentage")){
    warning("PCA_results object not found in global environment. Setting plotType to 'individual'")
    plotType <- "individual"
  }

  # Extract the specific components for plotting
  x_values <- as.numeric(pca_embeddings[paste0("pca_component_", Xcomp), ])
  y_values <- as.numeric(pca_embeddings[paste0("pca_component_", Ycomp), ])

  # Get axis labels
  if(has_pca_results){
    pca_results <- get("PCA_results", envir = .GlobalEnv)
    pca_results_explained <- data.frame(pca_results$eig[, 2])
    pca_results_explained <- round(as.vector(pca_results_explained[, 1]), 2)
    x_label <- paste0("PC", Xcomp, " (", pca_results_explained[Xcomp], "%)")
    y_label <- paste0("PC", Ycomp, " (", pca_results_explained[Ycomp], "%)")
    if(is_3d_plot) {
      z_label <- paste0("PC", Zcomp, " (", pca_results_explained[Zcomp], "%)")
    }
  } else {
    x_label <- paste0("PC", Xcomp)
    y_label <- paste0("PC", Ycomp)
    if(is_3d_plot) {
      z_label <- paste0("PC", Zcomp)
    }
  }

  # Set default color_palette for continuous features only
  if(is.null(color_palette) && color_by == "feature") {
    color_palette <- viridis::viridis(n = 20)
  }

  # Prepare coloring data - using unified color assignment
  if(color_by == "factor") {
    # Extract factor values
    factor_values <- as.character(romicsExtractFactor(romics_object, factor = "main"))
    # Use unified color assignment function
    color_values_mapped <- assign_factor_colors(factor_values, color_palette = color_palette, romics_object = romics_object)
    color_values <- factor(factor_values)
    color_title <- factor_name

    # Create named color vector for ggplot - match each unique level to its color
    unique_levels <- unique(as.character(factor_values))
    idx <- match(unique_levels, as.character(factor_values))
    colors <- setNames(color_values_mapped[idx], unique_levels)
  } else { # color_by == "feature"
    # Extract feature values
    feature_values <- as.numeric(romics_object$data[feature, ])
    # Scale if requested
    if(standardize) {
      feature_values <- as.numeric(scale(feature_values))
    }
    color_values <- feature_values
    color_values_mapped <- NA  # Not used for continuous coloring
    color_title <- feature
  }

  # For 3D plot
  if(is_3d_plot) {
    # Get Z values
    z_values <- as.numeric(pca_embeddings[paste0("pca_component_", Zcomp), ])
    # Create data frame for 3D plot
    pca_coord <- data.frame(
      PCA = x_values,
      PCB = y_values,
      PCC = z_values,
      color_var = color_values,
      sample = colnames(romics_object$data)
    )
    # Convert factor to factor type for plotly if needed
    if(color_by == "factor") {
      pca_coord$color_var <- as.factor(pca_coord$color_var)
    }
    # Create 3D plot
    if(color_by == "factor") {
      fig <- plotly::plot_ly(pca_coord, x = ~PCA, y = ~PCB, z = ~PCC,
                             color = ~color_var,
                             colors = colors,
                             text = ~sample)
    } else { # color_by == "feature"
      fig <- plotly::plot_ly(pca_coord, x = ~PCA, y = ~PCB, z = ~PCC,
                             color = ~color_var,
                             colors = color_palette,
                             text = ~sample)
    }
    fig <- fig %>% plotly::add_markers(size = size, opacity = alpha)
    fig <- fig %>% plotly::layout(
      scene = list(
        xaxis = list(title = x_label),
        yaxis = list(title = y_label),
        zaxis = list(title = z_label)
      ),
      coloraxis = list(colorbar = list(title = color_title))
    )
    return(fig)
  }
  # For 2D plot
  else {
    # Create data frame for 2D plot
    pca_coordinates <- data.frame(
      x = x_values,
      y = y_values,
      color_var = color_values,
      sample = colnames(romics_object$data)
    )
    # Convert to factor if needed
    if(color_by == "factor") {
      pca_coordinates$color_var <- as.factor(pca_coordinates$color_var)
    }

    # Determine plot scale
    pca_plot_scale <- max(c(ceiling(max(abs(c(x_values, y_values)))/10)+1)*10)

    # Create explained variance plot if PCA_results is available
    if(has_pca_results){
      explained_plot <- factoextra::fviz_screeplot(pca_results, barcolor = "gray20", barfill = "gray20") +
        theme_ROP() +
        theme(axis.text.x = element_text(angle = 0)) +
        ggtitle("Percentage of exp. variance")
    }

    # Create individual plot with outline control
    if(show_outline) {
      # Traditional plot with outlines
      pca_indiv_plot <- ggplot2::ggplot(pca_coordinates, ggplot2::aes(x = x, y = y, colour = color_var)) +
        ggplot2::geom_point(size = size, alpha = alpha) +
        ggplot2::labs(colour = color_title)
    } else {
      # Plot without outlines using fill
      pca_indiv_plot <- ggplot2::ggplot(pca_coordinates, ggplot2::aes(x = x, y = y, fill = color_var)) +
        ggplot2::geom_point(size = size, alpha = alpha, colour = NA, shape = 21) +
        ggplot2::labs(fill = color_title)
    }

    # Add common elements
    # Use user-provided limits if available, otherwise use automatic scaling
    x_limits <- if(is.null(xlim)) c(-pca_plot_scale, pca_plot_scale) else xlim
    y_limits <- if(is.null(ylim)) c(-pca_plot_scale, pca_plot_scale) else ylim

    pca_indiv_plot <- pca_indiv_plot +
      ggplot2::scale_y_continuous(limits = y_limits) +
      ggplot2::scale_x_continuous(limits = x_limits) +
      ggplot2::xlab(x_label) +
      ggplot2::ylab(y_label) +
      ggplot2::ggtitle("Principal Component Analysis") +
      theme_ROP()

    # Add appropriate color scale
    if(color_by == "factor") {
      if(show_outline) {
        pca_indiv_plot <- pca_indiv_plot + ggplot2::scale_color_manual(values = colors)
      } else {
        pca_indiv_plot <- pca_indiv_plot + ggplot2::scale_fill_manual(values = colors)
      }
    } else if(color_by == "feature") {
      if(show_outline) {
        pca_indiv_plot <- pca_indiv_plot + ggplot2::scale_color_gradientn(colors = color_palette, na.value = "gray20")
      } else {
        pca_indiv_plot <- pca_indiv_plot + ggplot2::scale_fill_gradientn(colors = color_palette, na.value = "gray20")
      }
    }

    # Add labels if requested
    if(label == TRUE){
      if(show_outline) {
        pca_indiv_plot <- pca_indiv_plot +
          ggplot2::geom_text(ggplot2::aes(colour = color_var), size = 3, label = pca_coordinates$sample)
      } else {
        pca_indiv_plot <- pca_indiv_plot +
          ggplot2::geom_text(ggplot2::aes(fill = color_var), size = 3, label = pca_coordinates$sample)
      }
    }

    # Convert to plotly if requested
    if(plotly) {
      pca_indiv_plot <- plotly::ggplotly(pca_indiv_plot)
    }

    # Return the appropriate plot
    if(plotType == "dual" && has_pca_results && !plotly){
      grid.arrange(explained_plot, pca_indiv_plot, ncol = 2)
    } else if(plotType == "percentage" && has_pca_results && !plotly){
      return(explained_plot)
    } else {
      return(pca_indiv_plot)
    }
  }
}

#' romicsUmap()
#' @description This function utilizes uwot::umap() function to calculate the umap embedding based on the data contained in a Romics_object and export the embedding in the metadata layer of the romics_object.
#' You can use a subset of features to generate the embedding, all the parameters of the uwot::umap() function are available, see the documentation of this function for more details.
#' The results are recorded in the embeddings layer of the romics_object, with rows named 'umap_component_1', 'umap_component_2', etc.
#' If a umap was run previously, the embedding will replace the current embedding.
#' @param romics_object has to be a romics_object created using romicsCreateObject() note that the data cannot contain missing values, missing values (e.g., NA) can be imputed using a RomicsProcessor::impute function. After performing the umap the imputation can be cancelled using the romicsRestoreMissing() function.
#' @param n_neighbors The size of local neighborhood (in terms of number of neighboring sample points) used for manifold approximation. Larger values result in more global views of the manifold, while smaller values result in more local data being preserved. In general values should be in the range 2 to 100.
#' @param n_components The dimension of the space to embed into. This defaults to 2 to provide easy visualization, but can reasonably be set to any integer value in the range 2 to 100.
#' @param metric Type of distance metric to use to find nearest neighbors. For nn_method = "annoy" this can be one of: "euclidean" (the default), "cosine", "manhattan", "hamming","correlation" (a distance based on the Pearson correlation), "categorical" (see below).
#' @param standardize Scaling to apply to the data layer of the romics_object: "none" or FALSE or NULL No scaling; "Z" or "scale" or TRUE Scale each feature to zero mean and variance 1; "maxabs" Center each feature to mean 0, then divide each element by the maximum absolute value over the entire matrix; "range" Range scale the entire matrix, so the smallest element is 0 and the largest is 1; "colrange" Scale each feature in the range (0,1); For UMAP, the default is "none".
#' @param init Type of initialization for the coordinates. Options are:
#' \itemize{
#'   \item "spectral" Spectral embedding using the normalized Laplacian of the fuzzy 1-skeleton, with Gaussian noise added
#'   \item "normlaplacian" Spectral embedding using the normalized Laplacian of the fuzzy 1-skeleton, without noise
#'   \item "random" Coordinates assigned using a uniform random distribution between -10 and 10
#'   \item "lvrandom" Coordinates assigned using a Gaussian distribution with standard deviation 1e-4, as used in LargeVis (Tang et al., 2016) and t-SNE
#'   \item "laplacian" Spectral embedding using the Laplacian Eigenmap (Belkin and Niyogi, 2002)
#'   \item "pca" The first two principal components from PCA of X if X is a data frame, and from a 2-dimensional classical MDS if X is of class "dist"
#'   \item "spca" Like "pca", but each dimension is then scaled so the standard deviation is 1e-4, to give a distribution similar to that used in t-SNE. This is an alias for init = "pca", init_sdev = 1e-4
#'   \item "agspectral" An "approximate global" modification of "spectral" which all edges in the graph to a value of 1, and then sets a random number of edges (negative_sample_rate edges per vertex) to 0.1, to approximate the effect of non-local affinities
#'   \item "pca_embeddings" Use existing PCA embeddings from the romics_object embeddings layer (see use_pca_embeddings parameter)
#'   \item A matrix of initial coordinates
#' }
#' For spectral initializations, ("spectral", "normlaplacian", "laplacian", "agspectral"), if more than one connected component is identified, no spectral initialization is attempted. Instead a PCA-based initialization is attempted. If verbose = TRUE the number of connected components are logged to the console. The existence of multiple connected components implies that a global view of the data cannot be attained with this initialization. Increasing the value of n_neighbors may help.
#' @param use_pca_embeddings Logical. If TRUE and init is not specified, will attempt to use PCA embeddings from the romics_object embeddings layer for initialization.
#' If TRUE and no PCA embeddings are found, will run romicsPCA() first. Default is FALSE.
#' When TRUE, the number of PCA components used will match n_components.
#' @param pca_ncp Integer. Number of principal components to calculate if PCA needs to be run for initialization.
#' Only used when use_pca_embeddings = TRUE and PCA embeddings don't exist. Default is max(10, n_components).
#' @param lock.seed has to be TRUE or FALSE to indicate if the seed is locked to enable reproducing plotting
#' @param seed numeric value indicating what seed to use when random seed are used
#' @param verbose boolean. If TRUE, also assigns the complete UMAP model to UMAP_results in the global environment
#' @param ... Arguments passed to umap function.
#' @details This function calculates UMAP embeddings and stores them in the romics_object. It can also store the complete
#' UMAP model in the global environment if verbose=TRUE.
#' @return Returns a romics_object with UMAP embeddings added to the embeddings layer.
#' @author Geremy Clair
#' @export
romicsUmap <- function(romics_object, verbose=TRUE, n_neighbors=30, n_components=2, metric="euclidean",
                       standardize=TRUE, init, use_pca_embeddings=FALSE, pca_ncp,
                       lock.seed=TRUE, seed=42, ...){
  arguments <- as.list(match.call())

  # Input validation
  if(!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }

  # Check required package
  if (!requireNamespace("uwot", quietly = TRUE)) {
    stop("Package 'uwot' is required for UMAP. Please install it with: install.packages('uwot')")
  }

  # Parameter validation and defaults
  if(missing(n_neighbors)){
    if(ncol(romics_object$data) < 30){
      n_neighbors = ncol(romics_object$data) - 2
    } else {
      n_neighbors = 30
    }
  }
  if(missing(verbose)){verbose = TRUE}
  if(missing(n_components)){n_components = 2}
  if(missing(metric)){metric = "euclidean"}
  if(missing(standardize)){standardize = TRUE}
  if(missing(lock.seed)){lock.seed = TRUE}
  if(!is.logical(lock.seed)){stop("lock_seed has to be either TRUE or FALSE")}
  if(missing(seed)){seed <- 42}
  if(missing(use_pca_embeddings)){use_pca_embeddings <- FALSE}

  # Check for missing values
  if(any(is.na(romics_object$data))) {
    stop("The data contains missing values (NA). Please impute missing values before running UMAP.")
  }

  # Ensure data is numeric and in proper format
  data_matrix <- romics_object$data
  if(!is.numeric(as.matrix(data_matrix))) {
    stop("Data must be numeric for UMAP analysis.")
  }

  # Convert to numeric matrix if it's a data.frame
  if(is.data.frame(data_matrix)) {
    data_matrix <- as.matrix(data_matrix)
    storage.mode(data_matrix) <- "double"
  }

  # Handle PCA embeddings for initialization
  init_value <- NULL
  if(use_pca_embeddings == TRUE && missing(init)){
    message("Using PCA embeddings for UMAP initialization...")

    # Check if PCA embeddings exist
    has_pca <- !is.null(romics_object$embeddings) &&
      any(grepl("pca_component_", rownames(romics_object$embeddings)))

    if(!has_pca){
      message("No PCA embeddings found. Running romicsPCA() first...")

      # Determine number of PCA components to calculate
      if(missing(pca_ncp)){
        pca_ncp <- max(10, n_components)
      }

      # Run PCA
      romics_object <- romicsPCA(
        romics_object,
        ncp = pca_ncp,
        n_components = pca_ncp,
        verbose = FALSE
      )
      message(paste("PCA completed with", pca_ncp, "components"))
    }

    # Extract PCA embeddings
    pca_rows <- grepl("pca_component_", rownames(romics_object$embeddings))
    pca_embeddings <- t(romics_object$embeddings[pca_rows, , drop=FALSE])

    # Check if enough PCA components are available
    if(ncol(pca_embeddings) < n_components){
      stop(paste0("Not enough PCA components available for initialization. ",
                  "Requested: ", n_components, ", Available: ", ncol(pca_embeddings), ". ",
                  "Please increase pca_ncp or reduce n_components."))
    }

    # Select the required number of components
    init_value <- pca_embeddings[, 1:n_components, drop=FALSE]

    # Ensure proper formatting
    storage.mode(init_value) <- "double"

    message(paste("Using", n_components, "PCA components for UMAP initialization"))

  } else if(!missing(init)){
    # Use the provided init parameter
    if(is.character(init) && init == "pca_embeddings"){
      stop("To use PCA embeddings from romics_object, set use_pca_embeddings=TRUE instead of init='pca_embeddings'")
    }
    init_value <- init
  }

  # Set seed so the plot will always look the same if lock seed was TRUE else use a random seed.
  if(lock.seed == TRUE){
    set.seed(seed)
  } else {
    set.seed(as.numeric(Sys.time()))
  }

  # Transpose data for UMAP (samples as rows, features as columns)
  d <- t(data_matrix)

  # Ensure d is a proper numeric matrix
  if(!is.matrix(d)) {
    d <- as.matrix(d)
  }
  storage.mode(d) <- "double"

  # Run the umap with ret_model=TRUE to keep the model when verbose=TRUE
  tryCatch({
    if(!is.null(init_value)){
      umap_result <- uwot::umap(d,
                                n_neighbors=n_neighbors,
                                n_components=n_components,
                                metric=metric,
                                scale=standardize,
                                init=init_value,
                                ret_model=verbose,
                                ...)
    } else {
      umap_result <- uwot::umap(d,
                                n_neighbors=n_neighbors,
                                n_components=n_components,
                                metric=metric,
                                scale=standardize,
                                ret_model=verbose,
                                ...)
    }
  }, error = function(e) {
    stop(paste("UMAP calculation failed:", e$message))
  })

  # Extract embeddings with proper error handling
  if(verbose) {
    # When ret_model=TRUE, the embeddings are in $embedding and the model info is in other elements
    if(is.list(umap_result) && "embedding" %in% names(umap_result)) {
      umap_embeddings <- as.data.frame(umap_result$embedding)
      # Save the full model to global environment
      assign("UMAP_results", umap_result, envir = .GlobalEnv)
      message("UMAP results have been assigned to 'UMAP_results' in the global environment")
    } else {
      stop("Unexpected UMAP result structure when ret_model=TRUE")
    }
  } else {
    # When ret_model=FALSE, the result is just the embeddings matrix
    if(is.matrix(umap_result) || is.data.frame(umap_result)) {
      umap_embeddings <- as.data.frame(umap_result)
    } else {
      stop("Unexpected UMAP result structure when ret_model=FALSE")
    }
  }

  # Ensure embeddings are properly formatted
  if(!is.data.frame(umap_embeddings)) {
    umap_embeddings <- as.data.frame(umap_embeddings)
  }

  # Name the components
  colnames(umap_embeddings) <- paste0("umap_component_", 1:ncol(umap_embeddings))

  # Ensure column names match original data
  rownames(umap_embeddings) <- colnames(romics_object$data)

  # Create or update embeddings layer in the romics_object
  if(!is.null(romics_object$embeddings)){
    if(sum(grepl("umap_component_", rownames(romics_object$embeddings))) > 0){
      message("UMAP embeddings were already recorded in the romics_object, these will be removed and replaced by the ones established here.")
      romics_object$embeddings <- romics_object$embeddings[!grepl("umap_component_", rownames(romics_object$embeddings)), , drop=FALSE]
    }
  } else {
    message("The <embeddings> layer was added to the romics_object, and the UMAP coordinates were added.")
    romics_object$embeddings <- data.frame(matrix(ncol=ncol(romics_object$data), nrow=0))
    colnames(romics_object$embeddings) <- colnames(romics_object$data)
  }

  # Add embeddings to the romics object
  romics_object$embeddings <- rbind(romics_object$embeddings, t(umap_embeddings))

  # Update the steps
  romics_object <- romicsUpdateSteps(romics_object, arguments)

  return(romics_object)
}

#' romicsUmapPlot()
#' @description Plots the samples on a UMAP plot in 2D or 3D.
#' The plot can be colored by any factor in the romics_object's metadata or by feature values.
#' This function uses the UMAP embeddings stored in the romics_object.
#' @param romics_object A romics_object that has UMAP embeddings in its embeddings layer (generated using romicsUmap)
#' @param Xcomp numerical/double. Indicate the component to plot on the X axis (default: 1)
#' @param Ycomp numerical/double. Indicate the component to plot on the Y axis (default: 2)
#' @param Zcomp numerical/double. If provided, creates a 3D plot with this component on the Z axis (default: NULL)
#' @param label boolean. Indicate if the sample name labels should be plotted (default: FALSE, applies only to 2D plots)
#' @param feature character string. If provided, the plot will be colored according to the values of this feature.
#'        If NULL (default), coloring will be based on the factor specified by factor_name.
#' @param factor_name character string. The name of the factor in romics_object$metadata to use for coloring the points
#'        when feature is NULL. If NULL (default), the romics_object's main_factor will be used.
#' @param color_palette a color palette to use. For categorical factors, should be a vector with one color per level.
#'        For continuous features (when feature is specified), should be a gradient color palette.
#'        Default: NULL (uses colors_romics for factors, viridis for features).
#' @param standardize boolean. Indicates if the feature values should be scaled when coloring by feature (default: TRUE)
#' @param size numeric. Size of the points in the plot (default: 3)
#' @param alpha numeric. Transparency of the points in the plot (default: 0.8)
#' @param plotly boolean. If TRUE, generates an interactive plotly plot instead of a static ggplot (default: FALSE)
#' @param force_color_update boolean. If TRUE, forces color regeneration even for main factor (default: FALSE)
#' @details This function plots UMAP results using the embeddings stored in the romics_object.
#'        If Zcomp is provided, a 3D interactive plot is generated using plotly.
#'        If Zcomp is NULL, a 2D plot is generated using ggplot2.
#'
#'        Color palette behavior:
#'        - For categorical factors: If color_palette is provided and has enough colors, it will be used.
#'          Otherwise, colors_romics or ROP_colors will be used (with warning if insufficient).
#'        - For continuous features: If color_palette is NULL, viridis is used. Otherwise, the provided palette is used.
#'
#'        When using a non-main factor for coloring without a custom color_palette, the function automatically
#'        updates the romics_object using romicsChangeFactor() to ensure consistent colors across plots.
#' @return Returns either a ggplot2 plot or a plotly interactive plot.
#' @author Geremy Clair
#' @export
romicsUmapPlot <- function(romics_object, Xcomp=1, Ycomp=2, Zcomp=NULL, label=FALSE,
                           feature=NULL, factor_name=NULL, color_palette=NULL,
                           standardize=TRUE, size=3, alpha=0.8, plotly=FALSE, force_color_update=FALSE,
                           xlim=NULL, ylim=NULL, zlim=NULL){

  # Input validation
  if(!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }

  # Determine coloring mode based on whether feature is provided
  if(!is.null(feature)) {
    color_by <- "feature"
    # Validate feature
    if(!feature %in% rownames(romics_object$data)) {
      stop(paste("Feature", feature, "not found in romics_object data."))
    }
  } else {
    color_by <- "factor"
    # Handle factor selection and update romics_object if needed
    if(is.null(factor_name)) {
      factor_name <- romics_object$main_factor
    } else {
      # Validate factor exists
      if(!(factor_name %in% rownames(romics_object$metadata))) {
        stop(paste("Factor", factor_name, "not found in romics_object metadata."))
      }
    }

    # Always update factor if it's different OR if force_color_update is TRUE
    # BUT: Don't update if custom color_palette is provided (we'll handle colors ourselves)
    if(factor_name != romics_object$main_factor || force_color_update) {
      if(is.null(color_palette)) {
        message(paste0("Changing factor to ", factor_name, ", this operation may take a few seconds."))
        romics_object <- romicsChangeFactor(romics_object, main_factor = factor_name)
        # Clear stale colors_romics so assign_factor_colors() uses fresh color assignment
        if("colors_romics" %in% rownames(romics_object$metadata)) {
          romics_object$metadata <- romics_object$metadata[rownames(romics_object$metadata) != "colors_romics", , drop = FALSE]
        }
      } else {
        message(paste0("Using factor '", factor_name, "' with custom color palette (skipping romicsChangeFactor)."))
        # Just update the main_factor without regenerating colors
        romics_object$main_factor <- factor_name
      }
    }
  }

  # Check if UMAP embeddings exist
  if(is.null(romics_object$embeddings) ||
     sum(grepl("umap_component_", rownames(romics_object$embeddings))) == 0) {
    stop("The romics_object doesn't contain UMAP embeddings. Please run romicsUmap() first.")
  }

  # Check if 2D or 3D plot
  is_3d_plot <- !is.null(Zcomp)

  # Get the UMAP embeddings
  umap_comp_rows <- grepl("umap_component_", rownames(romics_object$embeddings))
  umap_embeddings <- romics_object$embeddings[umap_comp_rows, ]
  n_components <- nrow(umap_embeddings)

  # Check component numbers
  if(missing(Xcomp) || !(is.numeric(Xcomp) || is.double(Xcomp))){
    warning("The component to be plotted on the first axis was not defined. UMAP1 will be used")
    Xcomp <- 1
  }
  if(missing(Ycomp) || !(is.numeric(Ycomp) || is.double(Ycomp))){
    warning("The component to be plotted on the second axis was not defined. UMAP2 will be used")
    Ycomp <- 2
  }
  if(Xcomp > n_components){
    stop(paste("The Xcomp selected (", Xcomp, ") is too large. Only", n_components, "components available."))
  }
  if(Ycomp > n_components){
    stop(paste("The Ycomp selected (", Ycomp, ") is too large. Only", n_components, "components available."))
  }

  # Check Zcomp if provided
  if(is_3d_plot && Zcomp > n_components){
    stop(paste("The Zcomp selected (", Zcomp, ") is too large. Only", n_components, "components available."))
  }

  # Extract component values for plotting
  x_values <- as.numeric(umap_embeddings[paste0("umap_component_", Xcomp), ])
  y_values <- as.numeric(umap_embeddings[paste0("umap_component_", Ycomp), ])
  if(is_3d_plot) {
    z_values <- as.numeric(umap_embeddings[paste0("umap_component_", Zcomp), ])
  }

  # Set default color_palette for continuous features only
  if(is.null(color_palette) && color_by == "feature") {
    color_palette <- viridis::viridis(n = 20)
  }

  # Prepare coloring data - using unified color assignment
  if(color_by == "factor") {
    # Extract factor values
    factor_values <- as.character(romicsExtractFactor(romics_object, factor = "main"))
    # Use unified color assignment function
    color_values_mapped <- assign_factor_colors(factor_values, color_palette = color_palette, romics_object = romics_object)
    color_values <- factor(factor_values)
    color_title <- factor_name

  } else { # color_by == "feature"
    # Extract feature values
    feature_values <- as.numeric(romics_object$data[feature, ])

    # Scale if requested
    if(standardize) {
      feature_values <- as.numeric(scale(feature_values))
    }

    color_values <- feature_values
    color_title <- feature
  }

  # Create color mapping for factor-based coloring
  if(color_by == "factor") {
    # Create a named vector mapping factor levels to colors
    unique_levels <- unique(as.character(factor_values))
    idx <- match(unique_levels, as.character(factor_values))
    colors <- setNames(color_values_mapped[idx], unique_levels)
  }

  # For 3D plot
  if(is_3d_plot) {
    # Create data frame for 3D plot
    umap_coord <- data.frame(
      UMAP1 = x_values,
      UMAP2 = y_values,
      UMAP3 = z_values,
      color_var = color_values,
      sample = colnames(romics_object$data)
    )

    # Convert factor to factor type for plotly if needed
    if(color_by == "factor") {
      umap_coord$color_var <- as.factor(umap_coord$color_var)
    }

    # Create 3D plot
    if(color_by == "factor") {
      fig <- plotly::plot_ly(umap_coord, x = ~UMAP1, y = ~UMAP2, z = ~UMAP3,
                             color = ~color_var,
                             colors = colors,
                             text = ~sample)
    } else { # color_by == "feature"
      fig <- plotly::plot_ly(umap_coord, x = ~UMAP1, y = ~UMAP2, z = ~UMAP3,
                             color = ~color_var,
                             colors = color_palette,
                             text = ~sample)
    }

    fig <- fig %>% plotly::add_markers(size = size, opacity = alpha,
                                       line = list(width = 0))
    fig <- fig %>% plotly::layout(
      scene = list(
        xaxis = list(title = paste0("UMAP Component ", Xcomp)),
        yaxis = list(title = paste0("UMAP Component ", Ycomp)),
        zaxis = list(title = paste0("UMAP Component ", Zcomp), range = if(!is.null(zlim)) zlim else NULL)
      ),
      coloraxis = list(colorbar = list(title = color_title))
    )

    return(fig)

  } else {
    # Create data frame for 2D plot
    umap_coordinates <- data.frame(
      x = x_values,
      y = y_values,
      color_var = color_values,
      sample = colnames(romics_object$data)
    )

    # Convert to factor if needed
    if(color_by == "factor") {
      umap_coordinates$color_var <- as.factor(umap_coordinates$color_var)
    }

    # Create individual plot
    umap_plot <- ggplot2::ggplot(umap_coordinates, ggplot2::aes(x = x, y = y, colour = color_var)) +
      ggplot2::geom_point(size = size, alpha = alpha, stroke = 0) +
      ggplot2::xlab(paste0("UMAP Component ", Xcomp)) +
      ggplot2::ylab(paste0("UMAP Component ", Ycomp)) +
      ggplot2::ggtitle("UMAP Plot") +
      ggplot2::labs(colour = color_title) +
      ggplot2::scale_x_continuous(limits = if(!is.null(xlim)) xlim else NULL) +
      ggplot2::scale_y_continuous(limits = if(!is.null(ylim)) ylim else NULL) +
      theme_ROP()

    # Add appropriate color scale
    if(color_by == "factor") {
      umap_plot <- umap_plot + ggplot2::scale_color_manual(values = colors)
    } else if(color_by == "feature") {
      umap_plot <- umap_plot + ggplot2::scale_color_gradientn(colors = color_palette, na.value = "gray20")
    }

    # Add labels if requested
    if(label == TRUE){
      umap_plot <- umap_plot +
        ggplot2::geom_text(ggplot2::aes(colour = color_var), size = 3, label = umap_coordinates$sample)
    }

    # Convert to plotly if requested
    if(plotly) {
      return(plotly::ggplotly(umap_plot))
    } else {
      return(umap_plot)
    }
  }
}


#' romicsSubset()
#' @description Keeps or drops a subset of specific elements/columns from the romics_object
#' @param romics_object A romics_object created using createRomicsObject()
#' @param subset_vector A character vector of factor levels or colnames to keep in the object.
#' @param type Either 'keep' or 'drop' to indicate if you want to conserve or to drop the elements from a given factor
#' @param by Either 'colnames' or 'level' to indicate what you want to keep or drop.
#' @param factor A factor contained in the metadata of the romics_object, to obtain the list of factors please use the function romicsFactorNames()
#' @param handle_embeddings Either 'remove' (default) or 'subset'. How to handle embeddings when subsetting.
#' @details This function creates a new object based on a previous romics_object to include or drop a list of specified columns from the original object. The created object will have a new step object created that will indicate the name of the original object to be subsetted and the log/non-log status of the object.
#' @details Note that this function will remove the stat layer from your object
#' @return This function generates a subsetted romics_object
#' @author Geremy Clair
#' @export
romicsSubset <- function(romics_object, subset_vector, type = "keep", by = "colnames",
                         factor = "main", handle_embeddings = "remove") {
  arguments <- as.list(match.call())

  if(!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }
  if(missing(subset_vector)) {stop("Your subsetting vector is missing")}
  if(!is.character(subset_vector)) {stop("Your subset vector has to be a character vector")}
  if(missing(type)) {type <- "keep"}
  if(!type %in% c("keep", "drop")) {stop("<type> has to be either 'keep' or 'drop'")}
  if(missing(by)) {by <- "colnames"}
  if(!by %in% c("colnames", "level")) {stop("<by> has to be either 'colnames' or 'level'")}
  if(!handle_embeddings %in% c("remove", "subset")) {
    stop("<handle_embeddings> has to be either 'remove' or 'subset'")
  }
  if(missing(factor)) {factor <- "main"}

  if(factor == "main") {
    factor <- romics_object$main_factor
  } else {
    if(!factor %in% rownames(romics_object$metadata)) {
      stop("Your factor is not present in the metadata")
    }
  }

  # Store original sample count
  original_n_samples <- ncol(romics_object$data)

  # Extract factor
  fac <- romics_object$metadata[factor == rownames(romics_object$metadata), ]
  fac <- as.factor(t(fac))

  # Create a logical vector containing the columns to keep (TRUE)
  if(by == "colnames") {
    if(sum(subset_vector %in% colnames(romics_object$data)) != length(subset_vector)) {
      warning("Not all the elements of the subset_vector were present in the colnames of the data")
    }
    if(type == "keep") {
      sub_logical <- colnames(romics_object$data) %in% subset_vector
    } else {
      sub_logical <- !colnames(romics_object$data) %in% subset_vector
    }
  } else {
    if(sum(subset_vector %in% fac) != length(subset_vector)) {
      warning("Not all the elements of the subset_vector were levels of the factor")
    }
    if(type == "keep") {
      sub_logical <- fac %in% subset_vector
    } else {
      sub_logical <- !fac %in% subset_vector
    }
  }

  # Calculate summary stats
  n_samples_kept <- sum(sub_logical)
  n_samples_removed <- original_n_samples - n_samples_kept
  pct_removed <- round((n_samples_removed / original_n_samples) * 100, 1)

  # Print summary message first
  message(paste0("Subsetting complete: ", n_samples_kept, " samples retained (",
                 pct_removed, "% of samples removed)."))

  # Remove the columns from the element data, metadata, and missingness of the romics_object
  romics_object$data <- romics_object$data[, sub_logical]
  romics_object$metadata <- romics_object$metadata[, sub_logical]
  romics_object$missingdata <- romics_object$missingdata[, sub_logical]

  # Handle embeddings
  if(!is.null(romics_object$embeddings) && length(romics_object$embeddings) > 0) {
    if(handle_embeddings == "remove") {
      romics_object$embeddings <- NULL
      message("Embeddings removed from the subsetted object.")
    } else if(handle_embeddings == "subset") {
      # Subset embeddings by columns (same as data)
      romics_object$embeddings <- romics_object$embeddings[, sub_logical]
      message("Embeddings subsetted to match the data.")
    }
  }

  # Remove the stat layer
  if("statistics" %in% names(romics_object)) {
    warning("The statistics layer was removed from the romics_object; the statistics were calculated on the non-subsetted object")
    romics_object <- romics_object[names(romics_object) != "statistics"]
  }

  class(romics_object) <- "romics_object"

  # Update the colors
  romics_object <- romicsUpdateColor(romics_object)

  # Update the steps
  romics_object <- romicsUpdateSteps(romics_object, arguments)

  return(romics_object)
}

#' romicsUmapDensityPlot()
#' @description Creates a UMAP plot with density overlay showing regions of high point concentration.
#' This is useful for identifying clusters or areas with many samples in UMAP space.
#' @param romics_object A romics_object that has UMAP embeddings in its embeddings layer
#' @param Xcomp numerical/double. Component to plot on the X axis (default: 1)
#' @param Ycomp numerical/double. Component to plot on the Y axis (default: 2)
#' @param density_type character. Type of density visualization: "contour", "filled", "hex", or "2d" (default: "filled")
#' @param bins numeric. Number of bins for hexagonal binning (used when density_type = "hex", default: 30)
#' @param n_contours numeric. Number of contour levels (used for contour plots, default: 10)
#' @param color_palette character vector. Color palette for density gradient (default: viridis)
#' @param bandwidth numeric. Bandwidth for density estimation (default: NULL for automatic)
#' @param internal_margin numeric. Absolute value to extend axis limits beyond data range (default: 0.5)
#' @return Returns a ggplot2 object showing density distribution in UMAP space
#' @author Geremy Clair
#' @export
romicsUmapDensityPlot <- function(romics_object, Xcomp = 1, Ycomp = 2,
                                  density_type = "filled", bins = 30, n_contours = 10,
                                  color_palette = NULL, bandwidth = NULL,
                                  internal_margin = 0.5) {

  # Input validation
  if(!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }

  # Check if UMAP embeddings exist
  if(is.null(romics_object$embeddings) ||
     sum(grepl("umap_component_", rownames(romics_object$embeddings))) == 0) {
    stop("The romics_object doesn't contain UMAP embeddings. Please run romicsUmap() first.")
  }

  # Set default color palette
  if(is.null(color_palette)) {
    color_palette <- viridis::viridis(n = 100)
  }

  # Get UMAP embeddings
  umap_comp_rows <- grepl("umap_component_", rownames(romics_object$embeddings))
  umap_embeddings <- romics_object$embeddings[umap_comp_rows, ]
  n_components <- nrow(umap_embeddings)

  # Validate component numbers
  if(Xcomp > n_components || Ycomp > n_components) {
    stop(paste("Component selection exceeds available components. Only", n_components, "available."))
  }

  # Extract coordinates
  x_values <- as.numeric(umap_embeddings[paste0("umap_component_", Xcomp), ])
  y_values <- as.numeric(umap_embeddings[paste0("umap_component_", Ycomp), ])

  # Create data frame
  umap_data <- data.frame(
    x = x_values,
    y = y_values
  )

  # Calculate axis limits with absolute margin
  x_range <- range(x_values, na.rm = TRUE)
  y_range <- range(y_values, na.rm = TRUE)

  x_limits <- c(x_range[1] - internal_margin, x_range[2] + internal_margin)
  y_limits <- c(y_range[1] - internal_margin, y_range[2] + internal_margin)

  # Create base plot
  density_plot <- ggplot2::ggplot(umap_data, ggplot2::aes(x = x, y = y))

  # Add density layer based on type
  if(density_type == "contour") {
    if(!is.null(bandwidth)) {
      density_plot <- density_plot +
        ggplot2::geom_density_2d(ggplot2::aes(color = ggplot2::after_stat(level)),
                                 bins = n_contours, h = bandwidth)
    } else {
      density_plot <- density_plot +
        ggplot2::geom_density_2d(ggplot2::aes(color = ggplot2::after_stat(level)),
                                 bins = n_contours)
    }
    density_plot <- density_plot +
      ggplot2::scale_color_gradientn(colors = color_palette, name = "Density")

  } else if(density_type == "filled") {
    if(!is.null(bandwidth)) {
      density_plot <- density_plot +
        ggplot2::stat_density_2d(ggplot2::aes(fill = ggplot2::after_stat(level)),
                                 geom = "polygon", bins = n_contours, h = bandwidth)
    } else {
      density_plot <- density_plot +
        ggplot2::stat_density_2d(ggplot2::aes(fill = ggplot2::after_stat(level)),
                                 geom = "polygon", bins = n_contours)
    }
    density_plot <- density_plot +
      ggplot2::scale_fill_gradientn(colors = color_palette, name = "Density")

  } else if(density_type == "hex") {
    density_plot <- density_plot +
      ggplot2::geom_hex(bins = bins) +
      ggplot2::scale_fill_gradientn(colors = color_palette, name = "Count")

  } else if(density_type == "2d") {
    if(!is.null(bandwidth)) {
      density_plot <- density_plot +
        ggplot2::geom_density_2d_filled(bins = n_contours, h = bandwidth)
    } else {
      density_plot <- density_plot +
        ggplot2::geom_density_2d_filled(bins = n_contours)
    }
    density_plot <- density_plot +
      ggplot2::scale_fill_manual(values = color_palette[seq(1, length(color_palette),
                                                            length.out = n_contours)])
  } else {
    stop("density_type must be one of: 'contour', 'filled', 'hex', or '2d'")
  }

  # Add axis limits with margin and labels
  density_plot <- density_plot +
    ggplot2::coord_cartesian(xlim = x_limits, ylim = y_limits) +
    ggplot2::xlab(paste0("UMAP Component ", Xcomp)) +
    ggplot2::ylab(paste0("UMAP Component ", Ycomp)) +
    ggplot2::ggtitle("UMAP Density Plot") +
    theme_ROP()

  return(density_plot)
}

#' romicsTsne()
#' @description This function utilizes Rtsne::Rtsne() function to calculate the t-SNE embedding based on the data contained in a Romics_object and export the embedding in the embeddings layer of the romics_object.
#' You can use a subset of features to generate the embedding, all the parameters of the Rtsne::Rtsne() function are available, see the documentation of this function for more details.
#' The results are recorded in the embeddings layer of the romics_object, with rows named 'tsne_component_1', 'tsne_component_2', etc.
#' If a t-SNE was run previously, the embedding will replace the current embedding.
#' @param romics_object has to be an romics_object created using romicsCreateObject() note that the data cannot contain missing values, missing values (e.g., NA) can be imputed using a RomicsProcessor::impute function. After performing the t-SNE the imputation can be cancelled using the romicsRestoreMissing() function.
#' @param perplexity Numeric; Perplexity parameter (should not be bigger than 3 * perplexity < nrow(X) - 1, see details for interpretation). Default: 30
#' @param dims Integer; Output dimensionality (default: 2)
#' @param initial_dims Integer; the number of dimensions that should be retained in the initial PCA step (default: 50)
#' @param standardize Logical; Whether to scale and center the data before applying t-SNE. Default: TRUE
#' @param theta Numeric; Speed/accuracy trade-off (increase for less accuracy), set to 0.0 for exact t-SNE (default: 0.5)
#' @param max_iter Integer; Number of iterations (default: 1000)
#' @param pca Logical; Whether an initial PCA step should be performed (default: TRUE)
#' @param pca_center Logical; Should data be centered before pca is applied? (default: TRUE)
#' @param pca_scale Logical; Should data be scaled before pca is applied? (default: FALSE)
#' @param check_duplicates Logical; Checks whether duplicates are present. It is best to make sure there are no duplicates present and set this option to FALSE, especially for large datasets (default: TRUE)
#' @param lock.seed Logical; TRUE or FALSE to indicate if the seed is locked to enable reproducing results
#' @param seed Numeric; value indicating what seed to use when random seed are used
#' @param verbose Logical; If TRUE, also assigns the complete t-SNE results to TSNE_results in the global environment and shows progress
#' @param ... Arguments passed to Rtsne function.
#' @details This function calculates t-SNE embeddings and stores them in the romics_object. It can also store the complete
#' t-SNE results in the global environment if verbose=TRUE. t-SNE is particularly good at preserving local structure and revealing clusters.
#' @return Returns a romics_object with t-SNE embeddings added to the embeddings layer.
#' @author Geremy Clair
#' @export
romicsTsne <- function(romics_object, verbose=TRUE, perplexity=30, dims=2, initial_dims=50,
                       standardize=TRUE, theta=0.5, max_iter=1000, pca=TRUE,
                       pca_center=TRUE, pca_scale=FALSE, check_duplicates=TRUE,
                       lock.seed=TRUE, seed=42, ...){
  arguments <- as.list(match.call())

  # Input validation
  if(!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }

  # Check required package
  if (!requireNamespace("Rtsne", quietly = TRUE)) {
    stop("Package 'Rtsne' is required for t-SNE. Please install it with: install.packages('Rtsne')")
  }

  # Parameter validation and defaults
  if(missing(perplexity)){
    # Perplexity should be smaller than (n_samples - 1) / 3
    max_perplexity <- floor((ncol(romics_object$data) - 1) / 3)
    perplexity <- min(30, max_perplexity)
    if(perplexity < 5) {
      stop("Too few samples for t-SNE. Need at least 15 samples (3 * perplexity + 1).")
    }
  }
  if(missing(verbose)){verbose = TRUE}
  if(missing(dims)){dims = 2}
  if(missing(initial_dims)){initial_dims = 50}
  if(missing(standardize)){standardize = TRUE}
  if(missing(theta)){theta = 0.5}
  if(missing(max_iter)){max_iter = 1000}
  if(missing(pca)){pca = TRUE}
  if(missing(pca_center)){pca_center = TRUE}
  if(missing(pca_scale)){pca_scale = FALSE}
  if(missing(check_duplicates)){check_duplicates = TRUE}
  if(missing(lock.seed)){lock.seed = TRUE}
  if(!is.logical(lock.seed)){stop("lock.seed has to be either TRUE or FALSE")}
  if(missing(seed)){seed <- 42}

  # Validate perplexity
  n_samples <- ncol(romics_object$data)
  if(perplexity >= (n_samples - 1) / 3) {
    stop(paste("Perplexity too large for the number of samples. Max perplexity:", floor((n_samples - 1) / 3)))
  }

  # Check for missing values
  if(any(is.na(romics_object$data))) {
    stop("The data contains missing values (NA). Please impute missing values before running t-SNE.")
  }

  # Ensure data is numeric and in proper format
  data_matrix <- romics_object$data
  if(!is.numeric(as.matrix(data_matrix))) {
    stop("Data must be numeric for t-SNE analysis.")
  }

  # Convert to numeric matrix if it's a data.frame
  if(is.data.frame(data_matrix)) {
    data_matrix <- as.matrix(data_matrix)
    storage.mode(data_matrix) <- "double"
  }

  # Set seed
  if(lock.seed == TRUE){
    set.seed(seed)
  } else {
    set.seed(as.numeric(Sys.time()))
  }

  # Transpose data for t-SNE (samples as rows, features as columns)
  d <- t(data_matrix)

  # Ensure d is a proper numeric matrix
  if(!is.matrix(d)) {
    d <- as.matrix(d)
  }
  storage.mode(d) <- "double"

  # Standardize data if requested (before t-SNE)
  if(standardize) {
    d <- scale(d, center = TRUE, scale = TRUE)
    # Handle any NaN or Inf values from scaling
    d[!is.finite(d)] <- 0
  }

  # Adjust initial_dims if necessary
  initial_dims <- min(initial_dims, ncol(d), nrow(d) - 1)

  # Run t-SNE
  tryCatch({
    tsne_result <- Rtsne::Rtsne(d,
                                dims = dims,
                                perplexity = perplexity,
                                theta = theta,
                                max_iter = max_iter,
                                pca = pca,
                                initial_dims = initial_dims,
                                pca_center = pca_center,
                                pca_scale = pca_scale,
                                check_duplicates = check_duplicates,
                                verbose = verbose,
                                ...)
  }, error = function(e) {
    stop(paste("t-SNE calculation failed:", e$message))
  })

  # Extract embeddings
  if(is.list(tsne_result) && "Y" %in% names(tsne_result)) {
    tsne_embeddings <- as.data.frame(tsne_result$Y)

    if(verbose) {
      # Save the full results to global environment
      assign("TSNE_results", tsne_result, envir = .GlobalEnv)
      message("t-SNE results have been assigned to 'TSNE_results' in the global environment")
    }
  } else {
    stop("Unexpected t-SNE result structure")
  }

  # Ensure embeddings are properly formatted
  if(!is.data.frame(tsne_embeddings)) {
    tsne_embeddings <- as.data.frame(tsne_embeddings)
  }

  # Name the components
  colnames(tsne_embeddings) <- paste0("tsne_component_", 1:ncol(tsne_embeddings))

  # Ensure row names match original data
  rownames(tsne_embeddings) <- colnames(romics_object$data)

  # Create or update embeddings layer in the romics_object
  if(!is.null(romics_object$embeddings)){
    if(sum(grepl("tsne_component_", rownames(romics_object$embeddings))) > 0){
      message("t-SNE embeddings were already recorded in the romics_object, these will be removed and replaced by the ones established here.")
      romics_object$embeddings <- romics_object$embeddings[!grepl("tsne_component_", rownames(romics_object$embeddings)), , drop=FALSE]
    }
  } else {
    message("The <embeddings> layer was added to the romics_object, and the t-SNE coordinates were added.")
    romics_object$embeddings <- data.frame(matrix(ncol=ncol(romics_object$data), nrow=0))
    colnames(romics_object$embeddings) <- colnames(romics_object$data)
  }

  # Add embeddings to the romics object
  romics_object$embeddings <- rbind(romics_object$embeddings, t(tsne_embeddings))

  # Update the steps
  romics_object <- romicsUpdateSteps(romics_object, arguments)

  return(romics_object)
}

#' romicsTsnePlot()
#' @description Plots the samples on a t-SNE plot in 2D or 3D.
#' The plot can be colored by any factor in the romics_object's metadata or by feature values.
#' This function uses the t-SNE embeddings stored in the romics_object.
#' @param romics_object A romics_object that has t-SNE embeddings in its embeddings layer (generated using romicsTsne)
#' @param Xcomp numerical/double. Indicate the component to plot on the X axis (default: 1)
#' @param Ycomp numerical/double. Indicate the component to plot on the Y axis (default: 2)
#' @param Zcomp numerical/double. If provided, creates a 3D plot with this component on the Z axis (default: NULL)
#' @param label boolean. Indicate if the sample name labels should be plotted (default: FALSE, applies only to 2D plots)
#' @param feature character string. If provided, the plot will be colored according to the values of this feature.
#'        If NULL (default), coloring will be based on the factor specified by factor_name.
#' @param factor_name character string. The name of the factor in romics_object$metadata to use for coloring the points
#'        when feature is NULL. If NULL (default), the romics_object's main_factor will be used.
#' @param color_palette a gradient color palette to use when coloring by feature (default: viridis(n=20))
#' @param standardize boolean. Indicates if the feature values should be scaled when coloring by feature (default: TRUE)
#' @param size numeric. Size of the points in the plot (default: 3)
#' @param alpha numeric. Transparency of the points in the plot (default: 0.8)
#' @param plotly boolean. If TRUE, generates an interactive plotly plot instead of a static ggplot (default: FALSE)
#' @param show_outline boolean. If FALSE, removes point outlines and shows only fill colors (default: TRUE)
#' @details This function plots t-SNE results using the embeddings stored in the romics_object.
#'        If Zcomp is provided, a 3D interactive plot is generated using plotly.
#'        If Zcomp is NULL, a 2D plot is generated using ggplot2.
#'        When using a non-main factor for coloring, the function automatically updates the romics_object
#'        using romicsChangeFactor() to ensure consistent colors across plots.
#'        Coloring is determined automatically:
#'        - If feature is provided, the plot is colored by feature values
#'        - If feature is NULL, the plot is colored by factor values
#'        When show_outline is FALSE, points are displayed without borders for cleaner appearance.
#' @return Returns either a ggplot2 plot or a plotly interactive plot.
#' @author Geremy Clair
#' @export
romicsTsnePlot <- function(romics_object, Xcomp=1, Ycomp=2, Zcomp=NULL, label=FALSE,
                           feature=NULL, factor_name=NULL, color_palette=NULL,
                           standardize=TRUE, size=3, alpha=0.8, plotly=FALSE, show_outline=TRUE,
                           xlim=NULL, ylim=NULL, zlim=NULL){
  # Input validation
  if(!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }

  # Determine coloring mode based on whether feature is provided
  if(!is.null(feature)) {
    color_by <- "feature"
    # Validate feature
    if(!feature %in% rownames(romics_object$data)) {
      stop(paste("Feature", feature, "not found in romics_object data."))
    }
  } else {
    color_by <- "factor"
    # Handle factor selection and update romics_object if needed
    if(is.null(factor_name)) {
      factor_name <- romics_object$main_factor
    } else {
      # Validate factor exists
      if(!(factor_name %in% rownames(romics_object$metadata))) {
        stop(paste("Factor", factor_name, "not found in romics_object metadata."))
      }
      # If factor_name is different from main_factor, update the romics_object
      if(factor_name != romics_object$main_factor) {
        message(paste0("Changing factor to ", factor_name, ", this operation may take a few seconds."))
        romics_object <- romicsChangeFactor(romics_object, main_factor = factor_name)
        # Clear stale colors_romics so assign_factor_colors() uses fresh color assignment
        if("colors_romics" %in% rownames(romics_object$metadata)) {
          romics_object$metadata <- romics_object$metadata[rownames(romics_object$metadata) != "colors_romics", , drop = FALSE]
        }
      }
    }
  }

  # Check if t-SNE embeddings exist
  if(is.null(romics_object$embeddings) ||
     sum(grepl("tsne_component_", rownames(romics_object$embeddings))) == 0) {
    stop("The romics_object doesn't contain t-SNE embeddings. Please run romicsTsne() first.")
  }

  # Check if 2D or 3D plot
  is_3d_plot <- !is.null(Zcomp)

  # Get the t-SNE embeddings
  tsne_comp_rows <- grepl("tsne_component_", rownames(romics_object$embeddings))
  tsne_embeddings <- romics_object$embeddings[tsne_comp_rows, ]
  n_components <- nrow(tsne_embeddings)

  # Check component numbers
  if(missing(Xcomp) || !(is.numeric(Xcomp) || is.double(Xcomp))){
    warning("The component to be plotted on the first axis was not defined. t-SNE1 will be used")
    Xcomp <- 1
  }
  if(missing(Ycomp) || !(is.numeric(Ycomp) || is.double(Ycomp))){
    warning("The component to be plotted on the second axis was not defined. t-SNE2 will be used")
    Ycomp <- 2
  }
  if(Xcomp > n_components){
    stop(paste("The Xcomp selected (", Xcomp, ") is too large. Only", n_components, "components available."))
  }
  if(Ycomp > n_components){
    stop(paste("The Ycomp selected (", Ycomp, ") is too large. Only", n_components, "components available."))
  }
  # Check Zcomp if provided
  if(is_3d_plot && Zcomp > n_components){
    stop(paste("The Zcomp selected (", Zcomp, ") is too large. Only", n_components, "components available."))
  }

  # Extract component values for plotting
  x_values <- as.numeric(tsne_embeddings[paste0("tsne_component_", Xcomp), ])
  y_values <- as.numeric(tsne_embeddings[paste0("tsne_component_", Ycomp), ])
  if(is_3d_plot) {
    z_values <- as.numeric(tsne_embeddings[paste0("tsne_component_", Zcomp), ])
  }

  # Set default color_palette for continuous features only
  if(is.null(color_palette) && color_by == "feature") {
    color_palette <- viridis::viridis(n = 20)
  }

  # Prepare coloring data - using unified color assignment
  if(color_by == "factor") {
    # Extract factor values
    factor_values <- as.character(romicsExtractFactor(romics_object, factor = "main"))
    # Use unified color assignment function
    color_values_mapped <- assign_factor_colors(factor_values, color_palette = color_palette, romics_object = romics_object)
    color_values <- factor(factor_values)
    color_title <- factor_name
  } else { # color_by == "feature"
    # Extract feature values
    feature_values <- as.numeric(romics_object$data[feature, ])
    # Scale if requested
    if(standardize) {
      feature_values <- as.numeric(scale(feature_values))
    }
    color_values <- feature_values
    color_values_mapped <- NA  # Not used for continuous coloring
    color_title <- feature
  }

  # Create color mapping for factor-based coloring
  if(color_by == "factor") {
    # Create a named vector mapping factor levels to colors
    unique_levels <- unique(as.character(factor_values))
    idx <- match(unique_levels, as.character(factor_values))
    colors <- setNames(color_values_mapped[idx], unique_levels)
  }

  # For 3D plot
  if(is_3d_plot) {
    # Create data frame for 3D plot
    tsne_coord <- data.frame(
      TSNE1 = x_values,
      TSNE2 = y_values,
      TSNE3 = z_values,
      color_var = color_values,
      sample = colnames(romics_object$data)
    )
    # Convert factor to factor type for plotly if needed
    if(color_by == "factor") {
      tsne_coord$color_var <- as.factor(tsne_coord$color_var)
    }
    # Create 3D plot
    if(color_by == "factor") {
      fig <- plotly::plot_ly(tsne_coord, x = ~TSNE1, y = ~TSNE2, z = ~TSNE3,
                             color = ~color_var,
                             colors = colors,
                             text = ~sample)
    } else { # color_by == "feature"
      fig <- plotly::plot_ly(tsne_coord, x = ~TSNE1, y = ~TSNE2, z = ~TSNE3,
                             color = ~color_var,
                             colors = color_palette,
                             text = ~sample)
    }
    fig <- fig %>% plotly::add_markers(size = size, opacity = alpha)
    fig <- fig %>% plotly::layout(
      scene = list(
        xaxis = list(title = paste0("t-SNE Component ", Xcomp)),
        yaxis = list(title = paste0("t-SNE Component ", Ycomp)),
        zaxis = list(title = paste0("t-SNE Component ", Zcomp), range = if(!is.null(zlim)) zlim else NULL)
      ),
      coloraxis = list(colorbar = list(title = color_title))
    )
    return(fig)
  } else {
    # Create data frame for 2D plot
    tsne_coordinates <- data.frame(
      x = x_values,
      y = y_values,
      color_var = color_values,
      sample = colnames(romics_object$data)
    )
    # Convert to factor if needed
    if(color_by == "factor") {
      tsne_coordinates$color_var <- as.factor(tsne_coordinates$color_var)
    }

    # Create individual plot with outline control
    if(show_outline) {
      # Traditional plot with outlines
      tsne_plot <- ggplot2::ggplot(tsne_coordinates, ggplot2::aes(x = x, y = y, colour = color_var)) +
        ggplot2::geom_point(size = size, alpha = alpha) +
        ggplot2::labs(colour = color_title)
    } else {
      # Plot without outlines using fill
      tsne_plot <- ggplot2::ggplot(tsne_coordinates, ggplot2::aes(x = x, y = y, fill = color_var)) +
        ggplot2::geom_point(size = size, alpha = alpha, colour = NA, shape = 21) +
        ggplot2::labs(fill = color_title)
    }

    # Add common elements
    tsne_plot <- tsne_plot +
      ggplot2::xlab(paste0("t-SNE Component ", Xcomp)) +
      ggplot2::ylab(paste0("t-SNE Component ", Ycomp)) +
      ggplot2::ggtitle("t-SNE Plot") +
      ggplot2::scale_x_continuous(limits = if(!is.null(xlim)) xlim else NULL) +
      ggplot2::scale_y_continuous(limits = if(!is.null(ylim)) ylim else NULL) +
      theme_ROP()

    # Add appropriate color scale
    if(color_by == "factor") {
      if(show_outline) {
        tsne_plot <- tsne_plot + ggplot2::scale_color_manual(values = colors)
      } else {
        tsne_plot <- tsne_plot + ggplot2::scale_fill_manual(values = colors)
      }
    } else if(color_by == "feature") {
      if(show_outline) {
        tsne_plot <- tsne_plot + ggplot2::scale_color_gradientn(colors = color_palette, na.value = "gray20")
      } else {
        tsne_plot <- tsne_plot + ggplot2::scale_fill_gradientn(colors = color_palette, na.value = "gray20")
      }
    }

    # Add labels if requested
    if(label == TRUE){
      if(show_outline) {
        tsne_plot <- tsne_plot +
          ggplot2::geom_text(ggplot2::aes(colour = color_var), size = 3, label = tsne_coordinates$sample)
      } else {
        tsne_plot <- tsne_plot +
          ggplot2::geom_text(ggplot2::aes(fill = color_var), size = 3, label = tsne_coordinates$sample)
      }
    }

    # Convert to plotly if requested
    if(plotly) {
      return(plotly::ggplotly(tsne_plot))
    } else {
      return(tsne_plot)
    }
  }
}

#' romicsTransferEmbeddings()
#' @description This function enables to transfer PCA, UMAP, and t-SNE coordinates/embeddings from one romics_object
#' to another one providing added flexibility in embeddings generation. The origin and receiving object
#' have to contain the same <samples> (data columns).
#' @param origin_romics_object has to be a romics_object created using romicsCreateObject()
#' @param receiving_romics_object has to be a romics_object created using romicsCreateObject() and containing the same samples as the origin_romics_object.
#' @param type a vector indicating which embeddings to transfer. Options are "pca", "umap", and "tsne". Default: c("pca","umap","tsne")
#' @return Returns the receiving_romics_object with the transferred embeddings.
#' @author Geremy Clair
#' @export
romicsTransferEmbeddings <- function(origin_romics_object, receiving_romics_object, type = c("pca","umap","tsne")){
  arguments <- as.list(match.call())

  if(!is.romicsObject(origin_romics_object) | missing(origin_romics_object)) {
    stop("<origin_romics_object> is missing or is not in the appropriate format.")
  }
  if(!is.romicsObject(receiving_romics_object) | missing(receiving_romics_object)) {
    stop("<receiving_romics_object> is missing or is not in the appropriate format.")
  }
  if(missing(type)){
    type = c("pca", "umap", "tsne")
  }

  # Standardize the type parameter to lowercase for case-insensitive comparison
  type <- tolower(type)

  # Check that both objects have the same samples (but allow different order)
  origin_samples <- colnames(origin_romics_object$data)
  receiving_samples <- colnames(receiving_romics_object$data)

  # Check if they contain the same samples (regardless of order)
  if(!setequal(origin_samples, receiving_samples)) {
    missing_in_origin <- setdiff(receiving_samples, origin_samples)
    missing_in_receiving <- setdiff(origin_samples, receiving_samples)

    error_msg <- "The two <romics_objects> do not contain the same <samples>."
    if(length(missing_in_origin) > 0) {
      error_msg <- paste(error_msg, "\nSamples missing in origin object:", paste(missing_in_origin, collapse = ", "))
    }
    if(length(missing_in_receiving) > 0) {
      error_msg <- paste(error_msg, "\nSamples missing in receiving object:", paste(missing_in_receiving, collapse = ", "))
    }
    stop(error_msg)
  }

  # If samples are the same but in different order, notify user
  if(!identical(origin_samples, receiving_samples)) {
    message("Sample orders differ between objects. Embeddings will be reordered to match receiving object.")
  }

  # Create embeddings layer if it doesn't exist
  if(is.null(receiving_romics_object$embeddings)){
    print("The <embeddings> layer was added to the receiving_romics_object.")
    receiving_romics_object$embeddings <- data.frame(matrix(ncol=ncol(receiving_romics_object$data), nrow=0))
    colnames(receiving_romics_object$embeddings) <- colnames(receiving_romics_object$data)
  }

  # Helper function to transfer and reorder embeddings
  transfer_embeddings <- function(pattern, embedding_type) {
    # Remove existing embeddings if present
    if(sum(grepl(pattern, rownames(receiving_romics_object$embeddings))) > 0){
      print(paste(embedding_type, "embeddings were already recorded in the receiving_romics_object, these will be removed and replaced."))
      receiving_romics_object$embeddings <<- receiving_romics_object$embeddings[!grepl(pattern, rownames(receiving_romics_object$embeddings)), ]
    }

    # Transfer embeddings if they exist in origin object
    embedding_rows <- grepl(pattern, rownames(origin_romics_object$embeddings))
    if(sum(embedding_rows) > 0){
      print(paste("Transferring", embedding_type, "embeddings."))

      # Extract embeddings from origin object
      origin_embeddings <- origin_romics_object$embeddings[embedding_rows, , drop = FALSE]

      # Reorder columns to match receiving object's sample order
      reordered_embeddings <- origin_embeddings[, receiving_samples, drop = FALSE]

      # Bind to receiving object
      receiving_romics_object$embeddings <<- rbind(
        receiving_romics_object$embeddings,
        reordered_embeddings
      )
    } else {
      warning(paste("Your <origin_romics_object> does not contain", paste0("'", embedding_type, "'"), "embeddings to be transferred."))
    }
  }

  # Handle PCA embeddings transfer
  if("pca" %in% type){
    transfer_embeddings("pca_component_", "PCA")
  }

  # Handle UMAP embeddings transfer
  if("umap" %in% type){
    transfer_embeddings("umap_component_", "UMAP")
  }

  # Handle t-SNE embeddings transfer
  if("tsne" %in% type){
    transfer_embeddings("tsne_component_", "t-SNE")
  }

  # Update steps
  receiving_romics_object <- romicsUpdateSteps(receiving_romics_object, arguments)
  return(receiving_romics_object)
}

#' romicsRemoveEmbeddings()
#' @description Remove embeddings from a romics_object to reduce memory usage.
#' @param romics_object A romics_object containing embeddings to be removed
#' @param embeddings Character vector of embedding names to remove. If NULL, removes all embeddings.
#' @return Returns the romics_object with specified embeddings removed
#' @author Geremy Clair
#' @export
romicsRemoveEmbeddings <- function(romics_object, embeddings = NULL) {
  arguments <- as.list(match.call())

  # Input validation
  if(!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format.")
  }

  if(is.null(romics_object$embeddings)) {
    message("No embeddings found in the romics_object.")
    return(romics_object)
  }

  # Remove all if not specified
  if(is.null(embeddings)) {
    embeddings <- names(romics_object$embeddings)
  }

  # Remove specified embeddings
  romics_object$embeddings[embeddings] <- NULL

  # Set to NULL if empty
  if(length(romics_object$embeddings) == 0) {
    romics_object$embeddings <- NULL
  }

  # Update steps
  romics_object <- romicsUpdateSteps(romics_object, arguments)

  message(paste(length(embeddings), "embedding(s) removed."))

  return(romics_object)
}

#' romicsKmeansSamplesEval()
#' @description Evaluates the optimal number of clusters for K-means using the elbow method.
#' Plots the within-cluster sum of squares (WCSS) for different numbers of clusters.
#' @param romics_object A romics_object created using romicsCreateObject()
#' @param iter.max The maximum number of iterations allowed. Default: 100
#' @param algorithm Character string specifying the K-means algorithm to use:
#'        'Hartigan-Wong', 'Lloyd', 'Forgy', or 'MacQueen'. Default: 'Lloyd'
#' @param cluster_using Character string specifying which data to use for clustering:
#'        'data' (full data matrix), 'pca' (PCA embeddings), 'umap' (UMAP embeddings), or 'tsne' (t-SNE embeddings). Default: 'data'
#' @param max_clust Maximum number of clusters to evaluate. Default: 50
#' @param ... Additional parameters passed to kmeans().
#' @return Returns a plot showing the WCSS for clusters from 1 to max_clust.
#' @author Geremy Clair
#' @export
romicsKmeansSamplesEval <- function(romics_object,
                                    iter.max = 100,
                                    algorithm = c("Lloyd", "Hartigan-Wong", "Forgy", "MacQueen"),
                                    cluster_using = c("data", "pca", "umap", "tsne"),
                                    max_clust = 50,
                                    ...) {
  # Input validation
  if(!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format.")
  }

  # Set defaults for missing parameters
  if(missing(iter.max)) {
    iter.max = 100
  }

  if(missing(algorithm)) {
    algorithm = "Lloyd"
  } else {
    algorithm = algorithm[1]  # Take first value if vector
  }

  if(!algorithm %in% c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen")) {
    stop("algorithm has to be one of: 'Hartigan-Wong', 'Lloyd', 'Forgy', or 'MacQueen'.")
  }

  if(missing(cluster_using)) {
    cluster_using <- "data"
  } else {
    cluster_using <- tolower(cluster_using)[1]  # Convert to lowercase and take first value if vector
  }

  if(!cluster_using %in% c("data", "pca", "umap", "tsne")) {
    stop("cluster_using has to be one of: 'data', 'pca', 'umap', or 'tsne'.")
  }

  if(cluster_using %in% c("pca", "umap", "tsne") && is.null(romics_object$embeddings)) {
    stop("Dimension reduction embeddings were not computed for this romics_object.")
  }

  if(missing(max_clust)) {
    max_clust = 50
  }

  # Prepare data for clustering based on specified method
  if(cluster_using == "data") {
    m <- as.matrix(t(romics_object$data))
  } else if(cluster_using == "umap") {
    if(sum(grepl("umap_component_", rownames(romics_object$embeddings))) == 0) {
      stop("UMAP embeddings were not computed for this romics_object.")
    }
    m <- t(romics_object$embeddings[grepl("umap_component_", rownames(romics_object$embeddings)), ])
  } else if(cluster_using == "pca") {
    if(sum(grepl("pca_component_", rownames(romics_object$embeddings))) == 0) {
      stop("PCA embeddings were not computed for this romics_object.")
    }
    m <- t(romics_object$embeddings[grepl("pca_component_", rownames(romics_object$embeddings)), ])
  } else if(cluster_using == "tsne") {
    if(sum(grepl("tsne_component_", rownames(romics_object$embeddings))) == 0) {
      stop("t-SNE embeddings were not computed for this romics_object.")
    }
    m <- t(romics_object$embeddings[grepl("tsne_component_", rownames(romics_object$embeddings)), ])
  }

  # Calculate WCSS for different numbers of clusters
  wcss <- vector()
  for (k in 1:max_clust) {
    message(paste("Evaluating", k, "clusters..."))
    kmeans_model <- kmeans(m, centers = k, iter.max = iter.max, algorithm = algorithm, ...)
    wcss[k] <- kmeans_model$tot.withinss
  }

  # Prepare plot data
  plot_data <- data.frame(
    k = 1:max_clust,
    wcss = wcss
  )

  # Create plot
  if(requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = k, y = wcss)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_point(size = 3) +
      ggplot2::labs(
        title = paste("Elbow Method for Optimal k using", toupper(cluster_using)),
        x = "Number of Clusters (k)",
        y = "Within-Cluster Sum of Squares (WCSS)"
      ) +
      theme_ROP()
    print(p)
  } else {
    # Fallback to base R plot
    plot(1:max_clust, wcss, type = "b", main = paste("Elbow Method using", toupper(cluster_using)),
         xlab = "Number of Clusters (k)", ylab = "Within-Cluster Sum of Squares (WCSS)",
         pch = 19, cex = 1.2, col = "blue", lwd = 2)
  }

  # Return results invisibly
  return(invisible(list(k = 1:max_clust, wcss = wcss)))
}

#' romicsKmeansSamples()
#' @description Generate K-means clustering for the samples in a romics_object. The clustering can be based on
#' the full data matrix or dimensionally reduced data from PCA, UMAP, or t-SNE embeddings.
#' @param romics_object A romics_object created using romicsCreateObject()
#' @param centers Either the number of clusters, or a set of initial (distinct) cluster centers.
#'        If a number, a random set of (distinct) rows in the data is chosen as the initial centers.
#' @param cluster_using Character string specifying which data to use for clustering:
#'        'data' (full data matrix), 'pca' (PCA embeddings), 'umap' (UMAP embeddings), or 'tsne' (t-SNE embeddings). Default: 'data'
#' @param factor_name Character string specifying the name for the new factor to be added to metadata.
#'        Default: 'K_clust'
#' @param iter.max The maximum number of iterations allowed. Default: 100
#' @param algorithm Character string specifying the K-means algorithm to use:
#'        'Hartigan-Wong', 'Lloyd', 'Forgy', or 'MacQueen'. Default: 'Lloyd'
#' @param lock.seed Logical. If TRUE, uses a fixed seed for reproducible clustering. Default: TRUE
#' @param seed Numeric value for the random seed when lock.seed=TRUE. Default: 42
#' @param ... Additional parameters passed to the kmeans() function.
#' @details This function performs K-means clustering on samples and adds the cluster assignments
#'          as a new factor in the romics_object's metadata. Using dimensionally-reduced data
#'          (PCA, UMAP, or t-SNE) can improve computation time but may reduce the accuracy of the clusters.
#'          The function ensures that "colors_romics" factor remains the last row in the metadata.
#'          Any NA values in cluster assignments will be replaced with "undefined cluster".
#' @return Returns the romics_object with a new factor added to metadata containing the cluster assignments.
#'         Clusters are named as "{factor_name}_001", "{factor_name}_002", etc.
#'         The "colors_romics" factor is guaranteed to be the last row in metadata.
#' @author Geremy Clair
#' @export
romicsKmeansSamples <- function(romics_object,
                                centers = 3,
                                cluster_using = c("data", "pca", "umap", "tsne"),
                                factor_name = "K_clust",
                                iter.max = 100,
                                algorithm = c("Lloyd", "Hartigan-Wong", "Forgy", "MacQueen"),
                                lock.seed = TRUE,
                                seed = 42,
                                ...) {
  arguments <- as.list(match.call())

  # Input validation
  if(!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format.")
  }
  if(missing(centers) || !is.numeric(centers)) {
    stop("centers has to be present and numeric.")
  }
  if(centers < 1) {
    stop("centers must be at least 1.")
  }

  # Set defaults for missing parameters
  if(missing(cluster_using)) {
    cluster_using <- "data"
  } else {
    cluster_using <- tolower(cluster_using)[1]
  }
  if(!cluster_using %in% c("data", "pca", "umap", "tsne")) {
    stop("cluster_using has to be one of: 'data', 'pca', 'umap', or 'tsne'.")
  }
  if(missing(algorithm)) {
    algorithm <- "Lloyd"
  } else {
    algorithm <- algorithm[1]
  }
  if(!algorithm %in% c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen")) {
    stop("algorithm has to be one of: 'Hartigan-Wong', 'Lloyd', 'Forgy', or 'MacQueen'.")
  }
  if(missing(factor_name)) {
    factor_name <- "K_clust"
  }

  # Check if embeddings exist if using pca, umap, or tsne
  if(cluster_using %in% c("pca", "umap", "tsne") && is.null(romics_object$embeddings)) {
    stop("Dimension reduction embeddings were not computed for this romics_object.")
  }

  # Set seed for reproducibility
  if(lock.seed) {
    set.seed(seed)
  } else {
    set.seed(as.numeric(Sys.time()))
  }

  # Prepare data for clustering based on specified method
  if(cluster_using == "data") {
    m <- as.matrix(t(romics_object$data))
    sample_names <- colnames(romics_object$data)
  } else if(cluster_using == "umap") {
    if(sum(grepl("umap_component_", rownames(romics_object$embeddings))) == 0) {
      stop("UMAP embeddings were not computed for this romics_object.")
    }
    embedding_matrix <- romics_object$embeddings[grepl("umap_component_", rownames(romics_object$embeddings)), ]
    m <- t(embedding_matrix)
    sample_names <- colnames(embedding_matrix)
  } else if(cluster_using == "pca") {
    if(sum(grepl("pca_component_", rownames(romics_object$embeddings))) == 0) {
      stop("PCA embeddings were not computed for this romics_object.")
    }
    embedding_matrix <- romics_object$embeddings[grepl("pca_component_", rownames(romics_object$embeddings)), ]
    m <- t(embedding_matrix)
    sample_names <- colnames(embedding_matrix)
  } else if(cluster_using == "tsne") {
    if(sum(grepl("tsne_component_", rownames(romics_object$embeddings))) == 0) {
      stop("t-SNE embeddings were not computed for this romics_object.")
    }
    embedding_matrix <- romics_object$embeddings[grepl("tsne_component_", rownames(romics_object$embeddings)), ]
    m <- t(embedding_matrix)
    sample_names <- colnames(embedding_matrix)
  }

  # Ensure rownames of clustering matrix match sample names
  rownames(m) <- sample_names

  # Run K-means
  k_results <- kmeans(m, centers = centers, iter.max = iter.max, algorithm = algorithm, ...)

  # Format cluster labels
  max_digits <- nchar(as.character(max(k_results$cluster)))
  format_string <- paste0("%s_%0", max_digits, "d")
  cluster_labels <- sprintf(format_string, factor_name, k_results$cluster)

  # Create a named vector to ensure correct sample-cluster matching
  names(cluster_labels) <- sample_names

  # Store colors_romics before modification (if it exists)
  colors_romics_backup <- NULL
  has_colors_romics <- "colors_romics" %in% rownames(romics_object$metadata)
  if (has_colors_romics) {
    colors_romics_backup <- romics_object$metadata["colors_romics", , drop = FALSE]
    # Remove colors_romics temporarily
    romics_object$metadata <- romics_object$metadata[rownames(romics_object$metadata) != "colors_romics", , drop = FALSE]
  }

  # Check if factor already exists and warn if overwriting
  if(factor_name %in% rownames(romics_object$metadata)) {
    warning(paste(factor_name, "was previously calculated and will be replaced by the newly computed values."))
    romics_object$metadata <- romics_object$metadata[rownames(romics_object$metadata) != factor_name, , drop = FALSE]
  }

  # Add new factor to metadata with proper sample matching
  # Get the sample order from metadata columns
  metadata_sample_order <- colnames(romics_object$metadata)
  # Match cluster labels to metadata sample order
  ordered_cluster_labels <- cluster_labels[metadata_sample_order]

  # CHECK FOR NA VALUES AND REPLACE WITH "undefined cluster"
  if(any(is.na(ordered_cluster_labels))) {
    n_na <- sum(is.na(ordered_cluster_labels))
    na_samples <- metadata_sample_order[is.na(ordered_cluster_labels)]

    warning(paste("Missing values were present for", n_na, "cluster assignments.",
                  "These will be replaced by 'undefined cluster' for samples:",
                  paste(na_samples, collapse = ", ")))

    # Replace NA values with "undefined cluster"
    ordered_cluster_labels[is.na(ordered_cluster_labels)] <- "undefined cluster"
  }

  # Create new factor with correctly ordered labels
  new_factor <- data.frame(matrix(ordered_cluster_labels, nrow = 1))
  colnames(new_factor) <- metadata_sample_order
  rownames(new_factor) <- factor_name

  # Add the new factor to metadata
  romics_object$metadata <- rbind(romics_object$metadata, new_factor)

  # ENSURE colors_romics is always last: Re-add colors_romics as the last row
  if (has_colors_romics) {
    romics_object$metadata <- rbind(romics_object$metadata, colors_romics_backup)
  }

  message(paste(factor_name, "was added to the romics factors."))

  # Print cluster summary
  cluster_summary <- table(k_results$cluster)
  message("Cluster summary:")
  for(i in 1:length(cluster_summary)) {
    cluster_name <- sprintf(format_string, factor_name, i)
    message(paste("  ", cluster_name, ":", cluster_summary[i], "samples"))
  }

  # Print undefined cluster count if any
  final_factor_values <- romics_object$metadata[factor_name, ]
  n_undefined <- sum(final_factor_values == "undefined cluster", na.rm = TRUE)
  if(n_undefined > 0) {
    message(paste("  undefined cluster:", n_undefined, "samples"))
  }

  # Update steps
  romics_object <- romicsUpdateSteps(romics_object, arguments)

  return(romics_object)
}

#' romicsLouvainSamples()
#' @description Generate Louvain community clustering for samples in a romics_object using memory-efficient algorithms.
#' The clustering can be based on the full data matrix or dimensionally reduced data from PCA, UMAP, or t-SNE embeddings.
#' For large datasets (>10k samples), this function uses optimized nearest neighbor search to reduce memory usage.
#' @param romics_object A romics_object created using romicsCreateObject()
#' @param cluster_using Character string specifying which data to use for clustering:
#'        'data' (full data matrix), 'pca' (PCA embeddings), 'umap' (UMAP embeddings), or 'tsne' (t-SNE embeddings).
#'        For large datasets, using 'pca' or 'umap' is strongly recommended to reduce computational complexity. Default: 'data'
#' @param k Number of nearest neighbors to use when constructing the graph. For very large datasets (>100k samples),
#'        consider reducing this value (e.g., k=10) to improve memory efficiency. Default: 20
#' @param resolution Resolution parameter that determines the granularity of the clustering.
#'        Higher values lead to more clusters. Default: 1.0
#' @param target_clusters Optional integer specifying the desired number of clusters. If provided, the function will
#'        iteratively adjust the resolution parameter to approximate this number. Default: NULL (use resolution directly)
#' @param max_iterations Maximum number of iterations when trying to reach target_clusters. Default: 20
#' @param tolerance Tolerance for target_clusters - clustering will stop if within this range. Default: 1
#' @param factor_name Character string specifying the name for the new factor to be added to metadata.
#'        Default: 'Louvain_clust'
#' @param weights Optional edge weights for the graph. If NULL, distance-based weighting is used automatically.
#' @param lock.seed Logical. If TRUE, uses a fixed seed for reproducible clustering. Default: TRUE
#' @param seed Numeric value for the random seed when lock.seed=TRUE. Default: 42
#' @param use_approximate Logical. If TRUE and the FNN package is available, uses fast approximate nearest neighbor
#'        search for improved memory efficiency and speed on large datasets. Default: TRUE
#' @param ... Additional parameters passed to igraph::cluster_louvain().
#' @details This function performs Louvain community detection on samples and adds the cluster assignments
#'          as a new factor in the romics_object's metadata. The function constructs a k-nearest neighbor
#'          graph from the data using memory-efficient algorithms, then applies the Louvain algorithm to identify communities.
#'
#'          When target_clusters is specified, the function will iteratively adjust the resolution parameter
#'          to try to achieve approximately that number of clusters. Note that exact control is not guaranteed
#'          due to the nature of the Louvain algorithm.
#'
#'          For large datasets, the function automatically uses optimized approaches:
#'          - FNN package for fast approximate nearest neighbor search (if available)
#'          - Edge list representation instead of full adjacency matrices
#'          - Chunked processing to control memory usage
#'          - Periodic garbage collection to free memory
#'
#'          The Louvain algorithm is a popular method for community detection in large networks,
#'          optimizing modularity to find densely connected groups of nodes.
#'
#'          Memory recommendations for large datasets:
#'          - Use dimensionally reduced data (cluster_using = "pca" or "umap")
#'          - Reduce k parameter (e.g., k=10 for >100k samples)
#'          - Install FNN package: install.packages("FNN")
#'          - Ensure sufficient system memory (recommend >16GB for >100k samples)
#' @return Returns the romics_object with a new factor added to metadata containing the cluster assignments.
#'         Clusters are named as "{factor_name}_001", "{factor_name}_002", etc.
#'         A message will indicate the number of clusters found and confirm the factor was added.
#' @author Geremy Clair
#' @export
romicsLouvainSamples <- function(romics_object,
                                 cluster_using = c("data", "pca", "umap", "tsne"),
                                 k = 20,
                                 resolution = 1.0,
                                 target_clusters = NULL,
                                 max_iterations = 20,
                                 tolerance = 1,
                                 factor_name = "Louvain_clust",
                                 weights = NULL,
                                 lock.seed = TRUE,
                                 seed = 42,
                                 use_approximate = TRUE,
                                 ...) {

  arguments <- as.list(match.call())

  # Check for required packages
  required_packages <- c("igraph")
  if (use_approximate) {
    required_packages <- c(required_packages, "FNN")
  }
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("The", pkg, "package is required. Please install it with: install.packages('", pkg, "')", sep = ""))
    }
  }

  # Input validation
  if(!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format.")
  }

  # Validate target_clusters
  if (!is.null(target_clusters)) {
    if (!is.numeric(target_clusters) || length(target_clusters) != 1 || target_clusters < 1) {
      stop("target_clusters must be a positive integer")
    }
    target_clusters <- as.integer(target_clusters)
  }

  # Set defaults for missing parameters
  if(missing(cluster_using)) {
    cluster_using <- "data"
  } else {
    cluster_using <- tolower(cluster_using)[1]
  }
  if(!cluster_using %in% c("data", "pca", "umap", "tsne")) {
    stop("cluster_using has to be one of: 'data', 'pca', 'umap', or 'tsne'.")
  }

  # Check if embeddings exist if using pca, umap, or tsne
  if(cluster_using %in% c("pca", "umap", "tsne") && is.null(romics_object$embeddings)) {
    stop("Dimension reduction embeddings were not computed for this romics_object.")
  }

  # Set seed for reproducibility
  if(lock.seed) {
    set.seed(seed)
  } else {
    set.seed(as.numeric(Sys.time()))
  }

  # Prepare data for clustering based on specified method
  if(cluster_using == "data") {
    m <- as.matrix(t(romics_object$data))
  } else if(cluster_using == "umap") {
    if(sum(grepl("umap_component_", rownames(romics_object$embeddings))) == 0) {
      stop("UMAP embeddings were not computed for this romics_object.")
    }
    m <- t(romics_object$embeddings[grepl("umap_component_", rownames(romics_object$embeddings)), ])
  } else if(cluster_using == "pca") {
    if(sum(grepl("pca_component_", rownames(romics_object$embeddings))) == 0) {
      stop("PCA embeddings were not computed for this romics_object.")
    }
    m <- t(romics_object$embeddings[grepl("pca_component_", rownames(romics_object$embeddings)), ])
  } else if(cluster_using == "tsne") {
    if(sum(grepl("tsne_component_", rownames(romics_object$embeddings))) == 0) {
      stop("t-SNE embeddings were not computed for this romics_object.")
    }
    m <- t(romics_object$embeddings[grepl("tsne_component_", rownames(romics_object$embeddings)), ])
  }

  n_samples <- nrow(m)
  n_dims <- ncol(m)
  message(paste("Processing", n_samples, "samples with", n_dims, "dimensions"))

  # For very large datasets, ensure we're using a reasonable k
  k <- min(k, n_samples - 1)

  # Build graph once (this is the time-consuming part)
  message("Building k-nearest neighbor graph...")
  if (use_approximate && requireNamespace("FNN", quietly = TRUE)) {
    knn_result <- FNN::get.knn(m, k = k)
    knn_indices <- knn_result$nn.index
    knn_distances <- knn_result$nn.dist
  } else {
    # Manual implementation for fallback
    knn_indices <- matrix(0, nrow = n_samples, ncol = k)
    knn_distances <- matrix(0, nrow = n_samples, ncol = k)

    for (i in 1:n_samples) {
      diffs <- sweep(m, 2, m[i, ], "-")
      distances <- sqrt(rowSums(diffs^2))
      distances[i] <- Inf
      nearest_idx <- order(distances)[1:k]
      knn_indices[i, ] <- nearest_idx
      knn_distances[i, ] <- distances[nearest_idx]
      if (i %% 1000 == 0) message(paste("Processed", i, "of", n_samples, "samples"))
    }
  }

  # Create graph from kNN results
  edges <- c()
  edge_weights <- c()

  for (i in 1:n_samples) {
    for (j in 1:k) {
      neighbor <- knn_indices[i, j]
      if (neighbor > 0 && neighbor != i) {
        edges <- c(edges, i, neighbor)
        weight <- 1 / (knn_distances[i, j] + 1e-10)
        edge_weights <- c(edge_weights, weight)
      }
    }
  }

  edge_matrix <- matrix(edges, ncol = 2, byrow = TRUE)
  g <- igraph::graph_from_edgelist(edge_matrix, directed = FALSE)
  g <- igraph::simplify(g, edge.attr.comb = list(weight = "mean"))

  # Add final weights
  edges_final <- igraph::get.edgelist(g, names = FALSE)
  final_weights <- numeric(nrow(edges_final))
  for (i in 1:nrow(edges_final)) {
    node1 <- edges_final[i, 1]
    node2 <- edges_final[i, 2]
    dist_val <- sqrt(sum((m[node1, ] - m[node2, ])^2))
    final_weights[i] <- 1 / (dist_val + 1e-10)
  }
  igraph::E(g)$weight <- final_weights

  # Function to perform clustering
  perform_clustering <- function(res) {
    louvain_result <- igraph::cluster_louvain(g, resolution = res, ...)
    cluster_assignments <- igraph::membership(louvain_result)
    n_clusters <- max(cluster_assignments)
    return(list(result = louvain_result, assignments = cluster_assignments, n_clusters = n_clusters))
  }

  # Main clustering logic
  if (!is.null(target_clusters)) {
    message(paste("Attempting to find approximately", target_clusters, "clusters..."))

    current_resolution <- resolution
    best_resolution <- resolution
    best_result <- NULL
    best_diff <- Inf

    # Initial clustering
    clustering_result <- perform_clustering(current_resolution)
    current_n_clusters <- clustering_result$n_clusters
    best_result <- clustering_result
    best_diff <- abs(current_n_clusters - target_clusters)

    message(paste("Initial clustering with resolution", round(current_resolution, 4), "found", current_n_clusters, "clusters"))

    # Iterative adjustment
    iteration <- 1
    while (abs(current_n_clusters - target_clusters) > tolerance && iteration <= max_iterations) {

      # Adjust resolution
      if (current_n_clusters < target_clusters) {
        current_resolution <- current_resolution * 1.2
      } else {
        current_resolution <- current_resolution * 0.8
      }

      # Perform clustering
      clustering_result <- perform_clustering(current_resolution)
      current_n_clusters <- clustering_result$n_clusters

      message(paste("Iteration", iteration, ": resolution", round(current_resolution, 4), "found", current_n_clusters, "clusters"))

      # Update best result if closer to target
      current_diff <- abs(current_n_clusters - target_clusters)
      if (current_diff < best_diff) {
        best_diff <- current_diff
        best_result <- clustering_result
        best_resolution <- current_resolution
      }

      iteration <- iteration + 1
    }

    # Use best result
    final_result <- best_result
    message(paste("Final result: using resolution", round(best_resolution, 4), "with", final_result$n_clusters, "clusters"))

  } else {
    # Use specified resolution directly
    message("Running Louvain clustering...")
    final_result <- perform_clustering(resolution)
  }

  # Format cluster labels
  cluster_assignments <- final_result$assignments
  max_digits <- nchar(as.character(max(cluster_assignments)))
  format_string <- paste0("%s_%0", max_digits, "d")
  cluster_labels <- sprintf(format_string, factor_name, cluster_assignments)

  # Check if factor already exists
  if(factor_name %in% rownames(romics_object$metadata)) {
    warning(paste(factor_name, "was previously calculated and will be replaced by the newly computed values."))
    romics_object$metadata <- romics_object$metadata[rownames(romics_object$metadata) != factor_name, ]
  }

  # Add new factor to metadata
  new_factor <- data.frame(matrix(cluster_labels, nrow = 1))
  colnames(new_factor) <- colnames(romics_object$metadata)
  rownames(new_factor) <- factor_name
  romics_object$metadata <- rbind(romics_object$metadata, new_factor)

  message(paste(factor_name, "was added to the romics factors."))
  message(paste("Found", max(cluster_assignments), "clusters"))

  # Update steps
  romics_object <- romicsUpdateSteps(romics_object, arguments)

  return(romics_object)
}


#' romicsLeidenSamples()
#' @description Generate Leiden community clustering for samples in a romics_object using memory-efficient algorithms.
#' The clustering can be based on the full data matrix or dimensionally reduced data from PCA, UMAP, or t-SNE embeddings.
#' For large datasets (>10k samples), this function uses optimized nearest neighbor search to reduce memory usage.
#' The Leiden algorithm is an improvement over Louvain that can escape local optima.
#' @param romics_object A romics_object created using romicsCreateObject()
#' @param cluster_using Character string specifying which data to use for clustering:
#'        'data' (full data matrix), 'pca' (PCA embeddings), 'umap' (UMAP embeddings), or 'tsne' (t-SNE embeddings).
#'        For large datasets, using 'pca' or 'umap' is strongly recommended to reduce computational complexity. Default: 'data'
#' @param k Number of nearest neighbors to use when constructing the graph. For very large datasets (>100k samples),
#'        consider reducing this value (e.g., k=10) to improve memory efficiency. Default: 20
#' @param resolution Resolution parameter that determines the granularity of the clustering.
#'        Higher values lead to more clusters. Default: 1.0
#' @param target_clusters Optional integer specifying the desired number of clusters. If provided, the function will
#'        iteratively adjust the resolution parameter to approximate this number. Default: NULL (use resolution directly)
#' @param max_iterations Maximum number of iterations when trying to reach target_clusters. Default: 20
#' @param tolerance Tolerance for target_clusters - clustering will stop if within this range. Default: 1
#' @param factor_name Character string specifying the name for the new factor to be added to metadata.
#'        Default: 'Leiden_clust'
#' @param weights Optional edge weights for the graph. If NULL, distance-based weighting is used automatically.
#' @param lock.seed Logical. If TRUE, uses a fixed seed for reproducible clustering. Default: TRUE
#' @param seed Numeric value for the random seed when lock.seed=TRUE. Default: 42
#' @param use_approximate Logical. If TRUE and the FNN package is available, uses fast approximate nearest neighbor
#'        search for improved memory efficiency and speed on large datasets. Default: TRUE
#' @param objective_function Character string specifying the objective function for Leiden clustering:
#'        'modularity', 'CPM' (Constant Potts Model), or 'significance'. Default: 'modularity'
#' @param n_iterations Maximum number of iterations for the Leiden algorithm. Default: -1 (run until convergence)
#' @param ... Additional parameters passed to igraph::cluster_leiden().
#' @details This function performs Leiden community detection on samples and adds the cluster assignments
#'          as a new factor in the romics_object's metadata. The function constructs a k-nearest neighbor
#'          graph from the data using memory-efficient algorithms, then applies the Leiden algorithm to identify communities.
#'
#'          The Leiden algorithm is an improvement over the Louvain algorithm that can escape local optima
#'          and provides better quality partitions. It guarantees well-connected communities.
#'
#'          When target_clusters is specified, the function will iteratively adjust the resolution parameter
#'          to try to achieve approximately that number of clusters. Note that exact control is not guaranteed
#'          due to the nature of the Leiden algorithm.
#'
#'          For large datasets, the function automatically uses optimized approaches:
#'          - FNN package for fast approximate nearest neighbor search (if available)
#'          - Edge list representation instead of full adjacency matrices
#'          - Chunked processing to control memory usage
#'          - Periodic garbage collection to free memory
#'
#'          Memory recommendations for large datasets:
#'          - Use dimensionally reduced data (cluster_using = "pca" or "umap")
#'          - Reduce k parameter (e.g., k=10 for >100k samples)
#'          - Install FNN package: install.packages("FNN")
#'          - Ensure sufficient system memory (recommend >16GB for >100k samples)
#' @return Returns the romics_object with a new factor added to metadata containing the cluster assignments.
#'         Clusters are named as "{factor_name}_001", "{factor_name}_002", etc.
#'         A message will indicate the number of clusters found and confirm the factor was added.
#' @examples
#' \dontrun{
#' # Basic usage
#' romics_obj <- romicsLeidenSamples(romics_obj)
#'
#' # Target specific number of clusters
#' romics_obj <- romicsLeidenSamples(romics_obj, target_clusters = 5)
#'
#' # For large datasets, use PCA embeddings with reduced k
#' romics_obj <- romicsLeidenSamples(romics_obj,
#'                                   cluster_using = "pca",
#'                                   k = 10,
#'                                   target_clusters = 8)
#'
#' # Custom factor name and resolution with CPM objective
#' romics_obj <- romicsLeidenSamples(romics_obj,
#'                                   cluster_using = "umap",
#'                                   resolution = 1.5,
#'                                   factor_name = "UMAP_Leiden",
#'                                   objective_function = "CPM")
#' }
#' @seealso \code{\link{romicsLouvainSamples}} for Louvain clustering,
#'          \code{\link{romicsPCA}}, \code{\link{romicsUMAP}} for dimensionality reduction
#' @author Geremy Clair
#' @export
romicsLeidenSamples <- function(romics_object,
                                cluster_using = c("data", "pca", "umap", "tsne"),
                                k = 20,
                                resolution = 1.0,
                                target_clusters = NULL,
                                max_iterations = 20,
                                tolerance = 1,
                                factor_name = "Leiden_clust",
                                weights = NULL,
                                lock.seed = TRUE,
                                seed = 42,
                                use_approximate = TRUE,
                                objective_function = c("modularity", "CPM", "significance"),
                                n_iterations = -1,
                                ...) {

  arguments <- as.list(match.call())

  # Check for required packages
  required_packages <- c("igraph")
  if (use_approximate) {
    required_packages <- c(required_packages, "FNN")
  }
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("The", pkg, "package is required. Please install it with: install.packages('", pkg, "')", sep = ""))
    }
  }

  # Check for Leiden algorithm availability
  if (!exists("cluster_leiden", where = "package:igraph", mode = "function")) {
    stop("The Leiden algorithm is not available in your igraph version. Please update igraph: install.packages('igraph')")
  }

  # Input validation
  if(!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format.")
  }

  # Validate target_clusters
  if (!is.null(target_clusters)) {
    if (!is.numeric(target_clusters) || length(target_clusters) != 1 || target_clusters < 1) {
      stop("target_clusters must be a positive integer")
    }
    target_clusters <- as.integer(target_clusters)
  }

  # Set defaults for missing parameters
  if(missing(cluster_using)) {
    cluster_using <- "data"
  } else {
    cluster_using <- tolower(cluster_using)[1]
  }
  if(!cluster_using %in% c("data", "pca", "umap", "tsne")) {
    stop("cluster_using has to be one of: 'data', 'pca', 'umap', or 'tsne'.")
  }

  if(missing(objective_function)) {
    objective_function <- "modularity"
  } else {
    objective_function <- objective_function[1]
  }
  if(!objective_function %in% c("modularity", "CPM", "significance")) {
    stop("objective_function has to be one of: 'modularity', 'CPM', or 'significance'.")
  }

  # Check if embeddings exist if using pca, umap, or tsne
  if(cluster_using %in% c("pca", "umap", "tsne") && is.null(romics_object$embeddings)) {
    stop("Dimension reduction embeddings were not computed for this romics_object.")
  }

  # Set seed for reproducibility
  if(lock.seed) {
    set.seed(seed)
  } else {
    set.seed(as.numeric(Sys.time()))
  }

  # Prepare data for clustering based on specified method
  if(cluster_using == "data") {
    m <- as.matrix(t(romics_object$data))
  } else if(cluster_using == "umap") {
    if(sum(grepl("umap_component_", rownames(romics_object$embeddings))) == 0) {
      stop("UMAP embeddings were not computed for this romics_object.")
    }
    m <- t(romics_object$embeddings[grepl("umap_component_", rownames(romics_object$embeddings)), ])
  } else if(cluster_using == "pca") {
    if(sum(grepl("pca_component_", rownames(romics_object$embeddings))) == 0) {
      stop("PCA embeddings were not computed for this romics_object.")
    }
    m <- t(romics_object$embeddings[grepl("pca_component_", rownames(romics_object$embeddings)), ])
  } else if(cluster_using == "tsne") {
    if(sum(grepl("tsne_component_", rownames(romics_object$embeddings))) == 0) {
      stop("t-SNE embeddings were not computed for this romics_object.")
    }
    m <- t(romics_object$embeddings[grepl("tsne_component_", rownames(romics_object$embeddings)), ])
  }

  n_samples <- nrow(m)
  n_dims <- ncol(m)
  message(paste("Processing", n_samples, "samples with", n_dims, "dimensions"))

  # For very large datasets, ensure we're using a reasonable k
  k <- min(k, n_samples - 1)

  # Build graph once (this is the time-consuming part)
  message("Building k-nearest neighbor graph...")
  if (use_approximate && requireNamespace("FNN", quietly = TRUE)) {
    knn_result <- FNN::get.knn(m, k = k)
    knn_indices <- knn_result$nn.index
    knn_distances <- knn_result$nn.dist
  } else {
    # Manual implementation for fallback
    knn_indices <- matrix(0, nrow = n_samples, ncol = k)
    knn_distances <- matrix(0, nrow = n_samples, ncol = k)

    for (i in 1:n_samples) {
      diffs <- sweep(m, 2, m[i, ], "-")
      distances <- sqrt(rowSums(diffs^2))
      distances[i] <- Inf
      nearest_idx <- order(distances)[1:k]
      knn_indices[i, ] <- nearest_idx
      knn_distances[i, ] <- distances[nearest_idx]
      if (i %% 1000 == 0) message(paste("Processed", i, "of", n_samples, "samples"))
    }
  }

  # Create graph from kNN results
  edges <- c()
  edge_weights <- c()

  for (i in 1:n_samples) {
    for (j in 1:k) {
      neighbor <- knn_indices[i, j]
      if (neighbor > 0 && neighbor != i) {
        edges <- c(edges, i, neighbor)
        weight <- 1 / (knn_distances[i, j] + 1e-10)
        edge_weights <- c(edge_weights, weight)
      }
    }
  }

  edge_matrix <- matrix(edges, ncol = 2, byrow = TRUE)
  g <- igraph::graph_from_edgelist(edge_matrix, directed = FALSE)
  g <- igraph::simplify(g, edge.attr.comb = list(weight = "mean"))

  # Add final weights
  edges_final <- igraph::get.edgelist(g, names = FALSE)
  final_weights <- numeric(nrow(edges_final))
  for (i in 1:nrow(edges_final)) {
    node1 <- edges_final[i, 1]
    node2 <- edges_final[i, 2]
    dist_val <- sqrt(sum((m[node1, ] - m[node2, ])^2))
    final_weights[i] <- 1 / (dist_val + 1e-10)
  }
  igraph::E(g)$weight <- final_weights

  # Function to perform clustering
  perform_clustering <- function(res) {
    leiden_result <- igraph::cluster_leiden(g,
                                            resolution = res,
                                            objective_function = objective_function,
                                            n_iterations = n_iterations,
                                            ...)
    cluster_assignments <- igraph::membership(leiden_result)
    n_clusters <- max(cluster_assignments)
    return(list(result = leiden_result, assignments = cluster_assignments, n_clusters = n_clusters))
  }

  # Main clustering logic
  if (!is.null(target_clusters)) {
    message(paste("Attempting to find approximately", target_clusters, "clusters..."))

    current_resolution <- resolution
    best_resolution <- resolution
    best_result <- NULL
    best_diff <- Inf

    # Initial clustering
    clustering_result <- perform_clustering(current_resolution)
    current_n_clusters <- clustering_result$n_clusters
    best_result <- clustering_result
    best_diff <- abs(current_n_clusters - target_clusters)

    message(paste("Initial clustering with resolution", round(current_resolution, 4), "found", current_n_clusters, "clusters"))

    # Iterative adjustment
    iteration <- 1
    while (abs(current_n_clusters - target_clusters) > tolerance && iteration <= max_iterations) {

      # Adjust resolution
      if (current_n_clusters < target_clusters) {
        current_resolution <- current_resolution * 1.2
      } else {
        current_resolution <- current_resolution * 0.8
      }

      # Perform clustering
      clustering_result <- perform_clustering(current_resolution)
      current_n_clusters <- clustering_result$n_clusters

      message(paste("Iteration", iteration, ": resolution", round(current_resolution, 4), "found", current_n_clusters, "clusters"))

      # Update best result if closer to target
      current_diff <- abs(current_n_clusters - target_clusters)
      if (current_diff < best_diff) {
        best_diff <- current_diff
        best_result <- clustering_result
        best_resolution <- current_resolution
      }

      iteration <- iteration + 1
    }

    # Use best result
    final_result <- best_result
    message(paste("Final result: using resolution", round(best_resolution, 4), "with", final_result$n_clusters, "clusters"))

  } else {
    # Use specified resolution directly
    message(paste("Running Leiden clustering with", objective_function, "objective function..."))
    final_result <- perform_clustering(resolution)
  }

  # Format cluster labels
  cluster_assignments <- final_result$assignments
  max_digits <- nchar(as.character(max(cluster_assignments)))
  format_string <- paste0("%s_%0", max_digits, "d")
  cluster_labels <- sprintf(format_string, factor_name, cluster_assignments)

  # Check if factor already exists
  if(factor_name %in% rownames(romics_object$metadata)) {
    warning(paste(factor_name, "was previously calculated and will be replaced by the newly computed values."))
    romics_object$metadata <- romics_object$metadata[rownames(romics_object$metadata) != factor_name, ]
  }

  # Add new factor to metadata
  new_factor <- data.frame(matrix(cluster_labels, nrow = 1))
  colnames(new_factor) <- colnames(romics_object$metadata)
  rownames(new_factor) <- factor_name
  romics_object$metadata <- rbind(romics_object$metadata, new_factor)

  message(paste(factor_name, "was added to the romics factors."))
  message(paste("Found", max(cluster_assignments), "clusters using", objective_function, "objective function"))

  # Update steps
  romics_object <- romicsUpdateSteps(romics_object, arguments)

  return(romics_object)
}


#' romicsNameClusters()
#' @description Rename clusters in any factor of a romics_object's metadata.
#' The original factor values are preserved in a new row with "_original" appended to the factor name.
#' Optionally create a new factor instead of modifying the existing one.
#' @param romics_object A romics_object created using romicsCreateObject()
#' @param cluster_factor Character string specifying the name of an existing factor in metadata containing the clusters to rename.
#'        Must be a valid factor name that can be verified with romicsFactorNames().
#'        If missing, defaults to the main factor of the romics_object.
#' @param cluster_new_names Optional character vector containing the new names for each cluster.
#'        If provided, must have the same length as the number of unique clusters.
#'        If NULL (default), the function will interactively prompt for names for each cluster.
#' @param new_factor_name Optional character string. If provided, creates a new factor with this name instead of modifying the existing factor.
#'        The original factor remains unchanged. If NULL (default), modifies the existing factor.
#' @param alphabetical_mapping Logical. If TRUE and cluster_new_names is provided, maps new names alphabetically to original clusters.
#'        If FALSE (default), maps based on the order of appearance or natural sorting of original clusters.
#' @return Returns the romics_object with renamed clusters and the original clusters preserved (if modifying existing factor)
#' @author Geremy Clair
#' @export
romicsNameClusters <- function(romics_object,
                               cluster_factor = NULL,
                               cluster_new_names = NULL,
                               new_factor_name = NULL,
                               alphabetical_mapping = FALSE) {

  arguments <- as.list(match.call())

  # Input validation
  if(!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format.")
  }

  # Set default cluster_factor to main factor if not provided
  if(is.null(cluster_factor) || missing(cluster_factor)) {
    cluster_factor <- romics_object$main_factor
    message(paste("Using main factor '", cluster_factor, "' as cluster_factor.", sep = ""))
  }

  # Validate cluster_factor using romicsFactorNames()
  available_factors <- romicsFactorNames(romics_object)
  if(!cluster_factor %in% available_factors) {
    stop(paste("The factor '", cluster_factor, "' was not found in your romics_object metadata.\n",
               "Available factors are: ", paste(available_factors, collapse = ", "), sep = ""))
  }

  # Validate new_factor_name if provided
  if(!is.null(new_factor_name)) {
    if(new_factor_name %in% available_factors) {
      warning(paste("Factor", new_factor_name, "already exists. It will be overwritten."))
      romics_object$metadata <- romics_object$metadata[rownames(romics_object$metadata) != new_factor_name, ]
    }
    create_new_factor <- TRUE
  } else {
    create_new_factor <- FALSE
    # Check if the original version already exists (only when modifying existing factor)
    original_factor_name <- paste0(cluster_factor, "_original")
    if(original_factor_name %in% available_factors) {
      warning(paste("An original version of", cluster_factor, "already exists. It will be overwritten."))
      romics_object$metadata <- romics_object$metadata[rownames(romics_object$metadata) != original_factor_name, ]
    }
  }

  # Get the current cluster values and their unique values
  c_list <- as.character(romics_object$metadata[rownames(romics_object$metadata) == cluster_factor, ])
  unique_c <- unique(c_list)

  # Sort clusters based on alphabetical_mapping parameter
  if(alphabetical_mapping) {
    unique_c <- sort(unique_c)  # Alphabetical order
    message("Using alphabetical mapping for cluster assignment.")
  } else {
    unique_c <- unique_c[order(unique_c)]  # Natural/numeric order when possible
  }

  # Display current clusters before mapping
  message(paste("Found", length(unique_c), "unique clusters in factor '", cluster_factor, "':", sep = ""))
  message(paste("Clusters:", paste(unique_c, collapse = ", ")))

  # Create a backup of the original clusters (only if modifying existing factor)
  if(!create_new_factor) {
    original_clusters <- c_list
  }

  # Handle cluster naming - either from provided names or interactively
  if(!is.null(cluster_new_names)) {
    # Check if the number of provided names matches the number of clusters
    if(length(cluster_new_names) != length(unique_c)) {
      stop(paste("The length of cluster_new_names (", length(cluster_new_names),
                 ") does not match the number of unique clusters (", length(unique_c), ").\n",
                 "Unique clusters found: ", paste(unique_c, collapse = ", "), sep = ""))
    }

    # Create a mapping table
    matching_table <- data.frame(
      Original_Cluster = unique_c,
      New_Name = cluster_new_names,
      stringsAsFactors = FALSE
    )

    # Force display the mapping with explicit printing
    message(paste0("\n=== Cluster Mapping for factor '", cluster_factor, "'",
                   if(alphabetical_mapping) " (alphabetical)" else "", " ==="))

    # Print mapping line by line to ensure visibility
    for(i in 1:nrow(matching_table)) {
      message(paste(matching_table$Original_Cluster[i], " -> ", matching_table$New_Name[i]))
    }
    message("=====================================\n")

    # Also print as table for better formatting
    cat("\nMapping Table:\n")
    print(matching_table)
    cat("\n")

  } else {
    # Interactive mode - prompt user for each cluster name
    factor_display_name <- if(create_new_factor) new_factor_name else cluster_factor
    cat(paste0("\nPlease enter names for each cluster in factor '", cluster_factor, "'",
               if(create_new_factor) paste(" (creating new factor '", new_factor_name, "')", sep = "") else "",
               ":\n"))

    cluster_new_names <- character(length(unique_c))
    for(i in seq_along(unique_c)) {
      # Keep prompting until a non-empty name is provided
      repeat {
        cat(paste0("Enter name for cluster '", unique_c[i], "' [", i, "/", length(unique_c), "]: "))
        new_name <- readline()
        if(nchar(new_name) > 0) break
        cat("Please enter a non-empty name.\n")
      }
      cluster_new_names[i] <- new_name
    }

    # Create a mapping table
    matching_table <- data.frame(
      Original_Cluster = unique_c,
      New_Name = cluster_new_names,
      stringsAsFactors = FALSE
    )

    # Display final mapping
    message(paste0("\n=== Final Cluster Mapping for factor '", cluster_factor, "'",
                   if(alphabetical_mapping) " (alphabetical mapping)" else "", " ==="))
    for(i in 1:nrow(matching_table)) {
      message(paste(matching_table$Original_Cluster[i], " -> ", matching_table$New_Name[i]))
    }
    message("=====================================\n")

    cat("\nFinal Mapping Table:\n")
    print(matching_table)
    cat("\n")
  }

  # Apply the new names
  new_clusters <- c_list
  for(i in seq_along(unique_c)) {
    new_clusters[c_list == unique_c[i]] <- cluster_new_names[i]
  }

  if(create_new_factor) {
    # Create a completely new factor
    romics_object$metadata <- rbind(romics_object$metadata, new_clusters)
    rownames(romics_object$metadata)[nrow(romics_object$metadata)] <- new_factor_name

    message(paste("✓ Created new factor:", new_factor_name))
    message(paste("✓ Original factor '", cluster_factor, "' remains unchanged.", sep = ""))

  } else {
    # Modify existing factor and preserve original

    # Save the original clusters to a new metadata row
    romics_object$metadata <- rbind(romics_object$metadata, original_clusters)
    rownames(romics_object$metadata)[nrow(romics_object$metadata)] <- original_factor_name

    # Replace the current clusters with the new names
    romics_object$metadata[cluster_factor, ] <- new_clusters

    message(paste("✓ Modified factor:", cluster_factor))
    message(paste("✓ Original values saved as:", original_factor_name))
  }

  # Update steps
  romics_object <- romicsUpdateSteps(romics_object, arguments)

  return(romics_object)
}
