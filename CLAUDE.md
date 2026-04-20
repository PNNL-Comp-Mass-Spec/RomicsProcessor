# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

RomicsProcessor is an R package for analyzing omics data (proteomics, metabolomics, etc.). The package provides:
- A core data structure (`romics_object`) that stores data, metadata, processing history, and visualization parameters
- Data import functions for common proteomics tools (MaxQuant, FragPipe, DiaNN, Scils)
- Data processing functions (normalization, transformation, missing value imputation)
- Statistical analysis and visualization functions
- Support for creating and applying reusable processing pipelines

## Development Setup

This is an R package project managed with RStudio and devtools. The project file is `RomicsProcessor.Rproj`.

### Key Dependencies
- **Core**: data.table, tidyr, dplyr
- **Visualization**: ggplot2, plotly, ComplexHeatmap, ggforce, ggrepel, ggtree
- **Statistics**: lme4, lmerTest, edgeR, FactoMineR
- **Dimensionality Reduction**: umap, uwot, Rtsne
- **Required**: R >= 4.0.0

Optional packages that enable specific features:
- `pmartR`: Additional omics utilities
- `sva`: Batch correction (SVA method)
- `leiden`: Community detection for network analysis

## Building & Documentation

### Documentation Generation
The package uses roxygen2 for documentation. All exported functions are documented with `#' ` roxygen comments in the R files.

To regenerate documentation after changing function signatures or adding @export tags:
```R
roxygen2::roxygenise()
```
This generates/updates the `NAMESPACE` file and `.Rd` files in the `man/` directory.

### Building the Package
To build the package in R:
```R
devtools::load_all()        # Load package for interactive testing
devtools::check()           # Run full package checks (includes documentation)
devtools::build()           # Build a source package tar.gz
devtools::install()         # Install locally
```

For command line builds:
```bash
R CMD build .               # Creates tarball
R CMD check RomicsProcessor_*.tar.gz  # Full validation
```

## Code Organization

The R code is organized by functional area in the `R/` directory (numbered files):

- **00_Romics_theme_plots.R**: Plotting theme and color palette definitions
- **01_OmicsData_Import.R**: Data import functions (extractMaxQuant, extractFragPipeData, extractDIANN, extractScils)
- **02_Romics_Base_Functions.R**: Core functions (createRomicsObject, is.romics_object, resetRomicsObject)
- **03_CV_calculations.R**: Coefficient of variation calculations
- **04_Manage_Missing.R**: Missing value imputation and handling
- **05_Transformation_normalization.R**: Data transformations (log, quantile norm) and normalization
- **06_Distribution_Plots.R**: Plotting functions for data distributions (violin, boxplot, histogram, ridge, etc.)
- **07_Grouping.R**: Sample grouping, PCA, clustering, and related plots
- **08_Heatmap_and_multifeature_plots.R**: Complex heatmaps and multi-feature visualizations
- **09_Statistics.R**: Statistical tests and analysis functions
- **10_Feature_Plots.R**: Individual feature/protein visualizations
- **11_sva_Batch_Correction.R**: Batch correction using SVA
- **12_PmartConvert.R**: Conversion to/from pmartR objects
- **13_Manage_Outliers.R**: Outlier detection and removal
- **14_Volcano_plot.R**: Volcano plot generation
- **15_Trend_analysis.R**: Trend analysis over time/conditions
- **16_Correlation_plot.R**: Correlation analysis and visualization
- **17_featureOperations.R**: Feature/protein level operations
- **18_Combine_objects.R**: Combining multiple romics_objects
- **19_metadata_plots_and_operations.R**: Metadata-related functions
- **20_spatial_maps.R**: Spatial mapping visualizations
- **21_new_scils_import.R**: Enhanced Scils format import
- **Example_data.R**: Example dataset generation

## Core Data Structure: romics_object

A `romics_object` is a list containing:
- `data`: Quantitative measurements (rows=features/proteins, columns=samples)
- `metadata`: Sample metadata (rows=metadata fields, columns=samples)
- `missingdata`: Boolean matrix tracking missing values
- `original_data` / `original_metadata`: Unmodified versions for comparison
- `main_factor`: Primary experimental factor for grouping
- `colors_romics`: Color assignments for visualization
- `custom_colors`: User-defined colors
- `steps`: Processing history (audit trail of applied functions)
- `omics_type`: Type of data (e.g., "proteomics", "metabolomics")
- `omics_information`: Processing metadata
- `uuid`: Unique object identifier
- `IDs`: Alternative feature/protein identifiers

### Typical Workflow
1. **Import data**: Use `extractMaxQuant()`, `extractFragPipeData()`, etc.
2. **Create romics_object**: `createRomicsObject(data, metadata, main_factor="...")`
3. **Process**: Apply functions like normalization, imputation, transformation
4. **Analyze**: Statistical tests, dimensionality reduction, visualization
5. **Create pipeline**: `createRomicsPipeline()` to capture the analysis workflow
6. **Reuse**: Apply the pipeline to new datasets with `applyRomicsPipeline()`

## Important Patterns & Conventions

### Function Naming
- Import functions: `extract*` (e.g., extractMaxQuant)
- Distribution plots: `distrib*` (e.g., distribBoxplot)
- Correlation functions: `*Correlation*` (e.g., FeatureCorrelation)
- Plot generation: Often end with "plot" or are called directly

### Argument Patterns
- Many functions accept a `romics_object` as the first argument
- Optional: `main_factor` to specify which metadata field to use for grouping
- Many plotting functions accept style parameters (colors, themes, etc.)

### Processing History
All functions that modify a romics_object add entries to the `$steps` list to maintain an audit trail. This is critical for reproducibility.

## Common Tasks

### Adding a new data import function
1. Create extraction function in `01_OmicsData_Import.R` following the pattern of existing extractors
2. Return data and metadata dataframes suitable for `createRomicsObject()`
3. Document with roxygen comments
4. Optionally add corresponding ID extraction function (e.g., `extractMaxQuantIDs`)

### Adding a new statistical test
1. Add function to `09_Statistics.R`
2. Accept romics_object and relevant parameters
3. Return a romics_object with results added to the processing history
4. Document effect on the object structure

### Adding a visualization function
1. Add to appropriate R file based on plot type
2. Accept romics_object and visualization parameters
3. Return ggplot object (or plotly for interactive plots)
4. Ensure consistency with existing color schemes and themes

### Modifying the romics_object structure
This is complex and should maintain backward compatibility:
1. Update `createRomicsObject()` to handle new fields with defaults
2. Update `resetRomicsObject()` if applicable
3. Consider impact on existing pipeline deserialization
4. Document changes thoroughly

## Testing & Validation

Example datasets are included as `.rda` files in the `data/` directory:
- `RomicsProcessor_Example_romics_data.rda`
- `RomicsProcessor_Example_romics_metadata.rda`
- `RomicsProcessor_Example_unprocessed_romics_object.rda`
- `RomicsProcessor_Example_processed_romics_object.rda`
- `RomicsProcessor_Example_romics_pipeline.rda`

These can be used for:
```R
load("data/RomicsProcessor_Example_romics_data.rda")
# Now ROP_data is available for testing
```

## Git Workflow

The repository uses a single master branch. When making changes:
1. Test locally with `devtools::load_all()` and manual testing
2. Run `devtools::check()` before committing
3. Update documentation if function signatures change
4. Include processing history entries when modifying romics_objects
5. Commit with clear messages describing what was changed and why

Recent work has focused on:
- FragPipe and DiaNN import functions
- Quantile normalization
- Bug fixes in data processing
- New visualization and analysis features
