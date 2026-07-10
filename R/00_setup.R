#' Install RomicsProcessor Dependencies
#' @description Installs all required dependencies for RomicsProcessor from CRAN and Bioconductor
#' @param bioconductor logical. Whether to install Bioconductor dependencies. Default: TRUE
#' @details This function should be run once after installing RomicsProcessor to set up all dependencies.
#' @return Invisible NULL
#' @examples
#' \dontrun{
#' installRomicsProcessorDependencies()
#' }
#' @export
installRomicsProcessorDependencies <- function(bioconductor = TRUE) {
  cran_dependencies <- c(
    "uuid", "ggplot2", "plotly", "dplyr", "tidyr", "reshape2", "scales",
    "ggrepel", "ggalluvial", "cowplot", "ggdendro", "dendextend", "RColorBrewer",
    "viridis", "doParallel", "lme4", "lmerTest", "emmeans", "FactoMineR",
    "factoextra", "missMDA", "uwot", "Rtsne", "igraph", "gplots", "circlize",
    "gridExtra", "patchwork", "FNN", "matrixStats", "ggsignif", "data.table"
  )

  bioc_dependencies <- c("ComplexHeatmap", "ggtree")

  message("Installing CRAN dependencies...")
  install.packages(cran_dependencies, dependencies = TRUE)

  if (bioconductor) {
    message("Installing Bioconductor dependencies...")
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(bioc_dependencies, ask = FALSE)
  }

  message("RomicsProcessor dependencies installation complete!")
  invisible(NULL)
}
