# RomicsProcessor

RomicsProcessor is an R package for analyzing bulk, single-cell, and spatial omics datasets.
The package provides a structured R object (`romics_object`) to store data, metadata, and complete processing history,
enabling reproducible and FAIR-compatible data analysis. RomicsProcessor enables complete analytical traceability through
UUID-based object tracking and automatic step logging, supporting reusable analytical pipelines and reliable research reproducibility.



## Installation

### Prerequisites

```R
install.packages(“devtools”)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ggtree")
BiocManager::install("ComplexHeatmap")
```

### Install from GitHub

```R
devtools::install_github("PNNL-Comp-Mass-Spec/RomicsProcessor")
```

### Package Installation On-Demand

RomicsProcessor uses a smart dependency system: when you call functions that require optional packages, the package will detect if they're installed and prompt you to install them if needed.

**Examples of functions that prompt for package installation:**

```R
# Batch correction with reference regions or ComBat/SVA methods
# → Will prompt to install: sva
romics_corrected <- romicsBatchCorrection(romics_obj, batch_factor=”batch”, 
                                          method=”ComBat”)

# Format conversion to/from Seurat objects
# → Will prompt to install: Seurat
romics_from_seurat <- SeuratToRomics(seurat_obj)

# Export to SpatialExperiment for BANKSY analysis
# → Will prompt to install: SpatialExperiment
spatial_exp <- romicsToSpatialExperiment(romics_spatial_obj)

# Advanced statistical utilities (pmartR integration)
# → Will prompt to install: pmartR
pmartR_obj <- romicsToPmartR(romics_obj)

# Advanced heatmap visualization
# → Will prompt to install: ComplexHeatmap
heatmap <- romicsComplexHeatmap(romics_obj, features=top_features)
```

When you encounter a package requirement, RomicsProcessor will provide clear installation instructions. No need to install packages until you actually use the functionality!

## System Requirements

- **R**: ≥ 4.5.1 (tested with R 4.5.1, recommended for case study reproducibility)
- **Memory**: For large datasets (>100k samples), recommend ≥16GB RAM
  - Use dimensionally reduced data (PCA/UMAP) with reduced k parameter for memory efficiency

## Getting Started - Case Study Examples

RomicsProcessor includes a comprehensive case study demonstrating the complete workflow:

```R
library(RomicsProcessor)

# Load the processed case study dataset
data(romics_proteins)

# Explore the processed data
head(romics_proteins$data)  # 1,247 proteins × 12 samples
romicsSteps(romics_proteins)  # View complete processing history

# Load the reproducible pipeline
data(romics_pipeline)

# Apply the pipeline to a new dataset
new_result <- applyRomicsPipeline(new_unprocessed_object, romics_pipeline)

# Check object relationships
relationship <- checkRelationRomicsObjects(romics_proteins, another_object)
print(relationship$summary)  # Understand how objects are related

# Generate your own pipeline
my_pipeline <- createRomicsPipeline(my_processed_object)
```

### Case Study: Bacillus cereus Proteomics
The `romics_proteins` dataset contains:
- **Organism**: Spore-forming foodborne pathogen *Bacillus cereus*
- **Conditions**: 4 growth media (soil extract, zucchini puree, Luria-Bertani, AOAC)
- **Samples**: 12 total (3 replicates × 4 conditions)
- **Features**: 1,247 proteins
- **Analysis**: Quality filtering, log2-transformation, median centering, statistical analysis
- **Key Findings**: Media-dependent proteome signatures, sporulation upregulation, toxin production patterns
- **Raw Data**: Available at MassIVE (MSV000085696) for full reproducibility

## Cite the code

To cite the package, please use the following DOI:
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

**Version 2.0.0** - A Major Milestone for Analytical Traceability and FAIR Science

This release completes RomicsProcessor's vision for FAIR-compliant omics analysis with:
- **Complete Object Encapsulation**: Full data, metadata, and processing history within each romics_object
- **Analytical Traceability**: UUID-based relationship tracking, branching detection, and complete lineage recording
- **Embedding Provenance**: Full source tracking for transferred embeddings with divergence documentation
- **Combined Object Lineage**: Complete traceability for combined objects with source UUID and step recording
- **Spatial Omics Support**: Full integration with spatial analysis frameworks and coordinates
- **Robust Reproducibility**: Validated pipeline generation and application with identical result reproduction

## Key Features

### Core Capabilities
- **Multi-Omics Support**: Handle bulk, single-cell, and spatial omics datasets (proteomics, metabolomics, genomics, etc.)
- **Complete Object Encapsulation**: Romics_objects contain data, metadata, statistics, processing history, embeddings, and dependencies
- **Reproducible Analysis**: Full processing history tracking with audit trail ensures complete reproducibility
- **FAIR Compliant**: Findable (UUID), Accessible (R package), Interoperable (format conversion), Reusable (pipelines)
- **Comprehensive Analysis**: Data import/transformation, normalization, dimensionality reduction, clustering, statistical analysis

### Analytical Traceability (v2.0 - NEW)
- **Object Relationship Detection**: `checkRelationRomicsObjects()` identifies branching, linear evolution, or unrelated objects
- **Embedding Provenance**: Full source tracking with divergence documentation for transferred embeddings
- **Combined Object Lineage**: Complete traceability recording which objects were combined and their divergent steps
- **UUID-Based Tracking**: Unique identifiers enable complete analytical chain reconstruction

### Format Interoperability
- **Seurat Integration**: Seamless bidirectional conversion with Seurat objects (SeuratToRomics, romicsToSeurat)
- **SpatialExperiment Support**: Export for BANKSY spatial analysis and integration with Bioconductor ecosystem
- **Embedding Preservation**: Automatic PCA, UMAP, and other dimensional reductions preserved during conversions

### Spatial Omics
- **Large-Scale Spatial Data**: Optimized for hundreds of thousands of data points (multiplexed immunofluorescence, MSI)
- **Spatial Map Generation**: Interactive visualization of spatial coordinates with feature overlays and clustering patterns
- **Spatial Integration**: Format conversion and coordinate preservation for spatial analysis frameworks

### Advanced Features
- **Batch Correction**: ComBat/SVA methods with optional reference region guidance
- **Custom Feature Labeling**: User-specified feature lists for volcano plot annotations
- **Advanced Visualization**: Interactive Shiny web interfaces, publication-quality plots, spatial maps
- **Reproducible Pipelines**: Generate and apply analysis pipelines for method reuse and validation

## v2.0.0 Highlights - Complete Analytical Traceability

### Analytical Traceability System
- 🔗 **UUID-Based Tracking**: Every romics_object has a unique identifier for complete lineage tracing
- 🌿 **Relationship Detection**: Identify if objects branched from common ancestor, are related through steps, or unrelated
- 📊 **Branching Point Detection**: Find exact point where analytical pipelines diverged
- 🔄 **Embedding Provenance**: Full source tracking when embeddings are transferred between objects
- 📦 **Combined Object Lineage**: Document which objects were combined and what divergent steps exist

### Enhanced Traceability Features
- ✅ `checkRelationRomicsObjects()` - Analyze relationships between any two romics_objects
- ✅ `romicsTransferEmbeddings()` - Records origin object UUID, name, and divergent steps
- ✅ `combineRomicsObjects()` - Tracks all source objects and their contributions
- ✅ `applyRomicsPipeline()` - Validate reproducibility with identical result reproduction

### Case Study Validation
The included case study demonstrates end-to-end reproducibility:
- Analyzed *Bacillus cereus* proteome across 4 growth media (12 samples)
- Applied analysis pipeline to unprocessed and reset objects
- Verified identical data, statistics, and embeddings
- Public raw data available (MassIVE MSV000085696)

Copyright 2019 Battelle Memorial Institute
