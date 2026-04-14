# RomicsProcessor

RomicsProcessor is an R package for analyzing omics data (proteomics, metabolomics, etc.).
The package provides a structured R object (`romics_object`) to store data, metadata, and processing history,
enabling reproducible and FAIR-compatible data analysis. RomicsProcessor also supports creating reusable
analytical pipelines from previously processed objects, allowing for rapid development and method reuse.

**Version 1.11.0** includes major performance optimizations for large datasets, cleaned dependency management,
and new naming convention standardization (camelCase for all exported functions).

## Key Features

- **Reproducible Analysis**: Full processing history tracking with audit trail
- **FAIR Compliant**: Structured data format supports FAIR principles
- **Reusable Pipelines**: Capture and apply processing workflows to new datasets
- **Optimized for Large Data**: Fast k-nearest neighbor clustering (20-100x speedup on large datasets)
- **Comprehensive Analysis**: Statistical tests, dimensionality reduction, visualization, batch correction support

## Recent Improvements (v1.11.0)

- ✨ **KNN Optimization**: Louvain and Leiden clustering now 20-100x faster on large datasets
- 🧹 **Cleaned Dependencies**: Removed 14 unused packages, added 8 missing ones; Imports reduced from 50 to 34
- 📝 **Naming Convention**: All functions now use camelCase (e.g., `is.romicsObject()`, `plotMap()`)
- 🐛 **Bug Fixes**: Fixed parameter mismatches in spatial mapping and single feature plots

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
[R markdown file](https://github.com/PNNL-Comp-Mass-Spec/RomicsProcessor/blob/master/Example/Bacillus_cereus_media_experiment.Rmd)
.

This example consist in an analysis of the proteome of B. cereus grown in different media. 
The Raw LC-MS/MS data is publicly available, on [MassIVE](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=6d0ca42ca79244a49d66a80fd741ba28), FTP at the following address (ftp://massive.ucsd.edu/MSV000085696/).

It was analyzed using MaxQuant set with the 
[attached parameters](https://github.com/PNNL-Comp-Mass-Spec/RomicsProcessor/blob/master/Example/parameters.txt)
. 

The MaxQuant output 
[proteinGroups.txt](https://github.com/PNNL-Comp-Mass-Spec/RomicsProcessor/blob/master/Example/proteinGroups.txt) 
Can be used as data and the 
[metadata](https://github.com/PNNL-Comp-Mass-Spec/RomicsProcessor/blob/master/Example/metadata.csv) 
is also provided.


The different examples dataset are accessible directly inside the package 
[data folder](https://github.com/PNNL-Comp-Mass-Spec/RomicsProcessor/tree/master/data).

## Cite the code

To cite the package please use the following DOI:
[![DOI](https://zenodo.org/badge/206400976.svg)](https://zenodo.org/badge/latestdoi/206400976)

## Contacts

Written by @GeremyClair for the Department of Energy (PNNL, Richland, WA) \
E-mail: geremy.clair@pnnl.gov or proteomics@pnnl.gov \
Website: https://omics.pnl.gov/ or https://panomics.pnnl.gov/

## License

RomicsProcessor is licensed under the 2-Clause BSD License; 
you may not use this file except in compliance with the License.  You may obtain 
a copy of the License at https://opensource.org/licenses/BSD-2-Clause

Copyright 2019 Battelle Memorial Institute
