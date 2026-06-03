# RomicsProcessor

RomicsProcessor is an R package for analyzing omics data (proteomics, metabolomics, etc.).
The package provides a structured R object (`romics_object`) to store data, metadata, and processing history,
enabling reproducible and FAIR-compatible data analysis. RomicsProcessor also supports creating reusable
analytical pipelines from previously processed objects, allowing for rapid development and method reuse.

**Version 1.17.0** adds comprehensive format conversion capabilities, including seamless interoperability with Seurat objects and SpatialExperiment formats, alongside enhanced spatial analysis and embeddings support.

## Key Features

- **Multi-Omics Support**: Handle bulk, single-cell, and spatial omics datasets (proteomics, metabolomics, genomics, etc.)
- **Reproducible Analysis**: Full processing history tracking with audit trail ensures complete reproducibility
- **FAIR Compliant**: Encapsulated romics_objects ensure Findability, Accessibility, Interoperability, and Reusability
- **Format Interoperability**: Seamless conversion with Seurat objects and SpatialExperiment for ecosystem integration
- **Large-Scale Spatial Data**: Optimized for hundreds of thousands of data points (multiplexed immunofluorescence, MSI)
- **Spatial Map Generation**: Interactive visualization of spatial coordinates with feature overlays, clustering patterns, and factor-based coloring
- **Comprehensive Analysis**: Data import/transformation, normalization, dimensionality reduction, clustering, statistical analysis
- **Advanced Visualization**: Interactive Shiny web interfaces, publication-quality plots, and interactive spatial maps
- **Batch Correction**: Combat/SVA methods for multi-batch omics integration

## Recent Improvements (v1.17.0)

- 🔄 **Format Conversion**: Convert between Seurat objects and romics_objects bidirectionally (SeuratToRomics, romicsToSeurat)
- 🗺️ **Spatial Integration**: Export romics_objects to SpatialExperiment format for BANKSY spatial analysis
- 📍 **Embedding Import**: Automatic PCA, UMAP, and other dimensional reductions preserved during conversions
- 🚀 **Performance**: Optimized metadata handling and vectorized matrix operations
- 📝 **Documentation**: Enhanced roxygen documentation with improved formatting and clarity

## Installation

### Prerequisites

```R
install.packages(“devtools”)
```

### Install from GitHub

```R
devtools::install_github(“PNNL-Comp-Mass-Spec/RomicsProcessor”)
```

### Optional Bioconductor Dependencies

For enhanced functionality, install these optional packages:

```R
# For batch correction (Combat/SVA method)
BiocManager::install(“sva”)

# For advanced omics utilities
install.packages(“pmartR”)

# For ComplexHeatmap and advanced visualization
BiocManager::install(“ComplexHeatmap”)
```

## System Requirements

- **R**: ≥ 4.0.0
- **Memory**: For large datasets (>100k samples), recommend ≥16GB RAM
  - Use dimensionally reduced data (PCA/UMAP) with reduced k parameter for memory efficiency

## Example of use

The folder /Example contain an 
example romics_objeft

## Cite the code

To cite the package please use the following DOI:
[![DOI](https://zenodo.org/badge/206400976.svg)](https://zenodo.org/badge/latestdoi/206400976)

## Authors

**Geremy Clair** (main contact) - geremy.clair@pnnl.gov \
Harsh Bhotika - harsh.bhotika@pnnl.gov \
Brittney Gorman - brittney.gorman@pnnl.gov

Written for the Department of Energy (PNNL, Richland, WA) \
Website: https://omics.pnl.gov/ or https://panomics.pnnl.gov/ \
E-mail: proteomics@pnnl.gov

## License

RomicsProcessor is licensed under the 2-Clause BSD License; 
you may not use this file except in compliance with the License.  You may obtain 
a copy of the License at https://opensource.org/licenses/BSD-2-Clause

Copyright 2019 Battelle Memorial Institute
