#' RomicsProcessor: An R package for analyzing omics data
#'
#' @description
#' RomicsProcessor provides a comprehensive toolkit for analyzing omics data (proteomics, metabolomics, etc.)
#' with support for data import, processing, statistical analysis, and visualization.
#'
#' @keywords internal
#' @import ggplot2
#' @import igraph
#' @import factoextra
#' @import dendextend
#' @import ggalluvial
#' @import plotly
#' @import ComplexHeatmap
#' @import ggsignif
#' @import ggtree
"_PACKAGE"

#' @noRd
.onAttach <- function(libname, pkgname) {
  require("ggplot2", quietly = TRUE)
  require("igraph", quietly = TRUE)
  require("factoextra", quietly = TRUE)
  require("dendextend", quietly = TRUE)
  require("ggalluvial", quietly = TRUE)
  require("plotly", quietly = TRUE)
  require("ComplexHeatmap", quietly = TRUE)
  require("ggsignif", quietly = TRUE)
  require("ggtree", quietly = TRUE)
}
