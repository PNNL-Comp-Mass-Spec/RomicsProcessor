.onAttach <- function(libname, pkgname) {
  cran_dependencies <- c(
    "uuid",
    "ggplot2",
    "plotly",
    "dplyr",
    "tidyr",
    "reshape2",
    "scales",
    "ggrepel",
    "ggalluvial",
    "cowplot",
    "ggdendro",
    "dendextend",
    "RColorBrewer",
    "viridis",
    "parallel",
    "doParallel",
    "lme4",
    "lmerTest",
    "emmeans",
    "FactoMineR",
    "factoextra",
    "missMDA",
    "uwot",
    "Rtsne",
    "igraph",
    "gplots",
    "circlize",
    "gridExtra",
    "grid",
    "patchwork",
    "FNN",
    "matrixStats",
    "ggsignif",
    "data.table"
  )

  bioc_dependencies <- c(
    "ComplexHeatmap",
    "ggtree"
  )

  missing_cran <- setdiff(cran_dependencies, rownames(installed.packages()))
  missing_bioc <- setdiff(bioc_dependencies, rownames(installed.packages()))

  if (length(missing_cran) > 0) {
    packageStartupMessage("Installing missing CRAN dependencies for RomicsProcessor...")
    suppressWarnings(suppressMessages(
      install.packages(missing_cran, repos = getOption("repos"),
                      quiet = TRUE, dependencies = TRUE)
    ))
  }

  if (length(missing_bioc) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      packageStartupMessage("Installing BiocManager for Bioconductor dependencies...")
      suppressWarnings(suppressMessages(
        install.packages("BiocManager", quiet = TRUE)
      ))
    }
    packageStartupMessage("Installing missing Bioconductor dependencies for RomicsProcessor...")
    suppressWarnings(suppressMessages(
      BiocManager::install(missing_bioc, ask = FALSE, quiet = TRUE)
    ))
  }

  all_dependencies <- c(cran_dependencies, bioc_dependencies)
  for (pkg in all_dependencies) {
    try(library(pkg, character.only = TRUE), silent = TRUE)
  }
}
