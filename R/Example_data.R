# Case Study 1 Data: Bacillus cereus Proteome Analysis
#
# This file documents the example datasets from the RomicsProcessor case study:
# - romics_proteins: Processed proteomics dataset with complete analysis history
# - romics_pipeline: Reproducible pipeline extracted from the analysis
#
# These datasets are stored as .rda files in the data/ folder and are automatically
# loaded via lazy loading when the package is attached.
#
# CASE STUDY OVERVIEW:
# ====================
# This case study demonstrates RomicsProcessor's end-to-end workflow for reproducible
# omics analysis and biological discovery. We analyzed the proteome of the spore-forming
# foodborne pathogen Bacillus cereus grown in four different media conditions:
# - Soil extract (mimicking environmental conditions)
# - Zucchini puree (food matrix)
# - Luria-Bertani medium (laboratory standard)
# - AOAC limited media (gut-environment mimic)
#
# WORKFLOW:
# =========
# 1. Raw mass spectrometry data was pre-processed with FragPipe
# 2. Data imported into RomicsProcessor with preprocessing settings logged in
#    the omics_information layer
# 3. Quality filtering, log2-transformation, and median centering applied
# 4. All steps, parameters, and timestamps automatically captured in the steps layer
# 5. Statistical analysis identified media-specific proteome signatures
#
# KEY FINDINGS:
# =============
# - Soil extract media increased abundance of sporulation-related proteins (GO 0030435)
# - Gut-mimicking media (AOAC) promoted production of toxins (Hemolysin BL, enterotoxin Nhe)
# - High within-group correlation (Pearson > 0.95) and clear between-group separation
# - Demonstrates that RomicsProcessor can derive meaningful insights into bacterial
#   pathogenesis, with findings corroborating previously documented responses
#
# REPRODUCIBILITY VALIDATION:
# ============================
# The analysis was validated using RomicsProcessor's reproducibility features:
# 1. A duplicate object was created at import time for reproducibility testing
# 2. A pipeline was generated after analysis completion using createRomicsPipeline()
# 3. The pipeline was successfully applied to:
#    - An unprocessed object (romics_unprocessed) created at import time
#    - A reset object (romics_restored_original) recreated from raw data
# 4. In both cases, identical data values and statistical tables were produced
# 5. This provides computational demonstration of end-to-end reproducibility
#
# DATA ACCESS:
# ============
# Raw mass spectrometry files are publicly available at MassIVE:
# - Repository ID: MSV000085696
# - FTP: ftp://massive.ucsd.edu/MSV000085696/
# - This ensures full research transparency and reproducibility
#
# ROMICS_PROTEINS - Processed Proteomics Object:
# ===============================================
# This romics_object contains the complete processed dataset from the case study:
# - 1,247 unique proteins (features)
# - 12 samples (4 media types × 3 replicates)
# - Contains all processed data and complete analysis history
# - Ready for immediate use in downstream analyses
#
# Structure:
#   data: Normalized log2-transformed abundance values (proteins × samples)
#   metadata: Sample groupings by media type and experimental metadata
#   statistics: Statistical test results (ANOVA/t-tests across media)
#   embeddings: PCA and UMAP dimensional reductions
#   steps: Complete processing history with timestamps and parameters
#   omics_information: FragPipe preprocessing configuration logged during import
#   custom_colors: Color scheme for media types
#   uuid: Unique identifier for this specific processed dataset
#   dependencies: Packages and versions used in analysis
#
# ROMICS_PIPELINE - Reproducible Analysis Pipeline:
# ==================================================
# A pipeline object generated from romics_proteins that can be applied to:
# - Unprocessed versions of the same dataset
# - New datasets with identical experimental design
# - Alternative data processing branches for comparative analysis
#
# The pipeline captures the complete processing history including:
# - Quality filtering parameters and thresholds
# - Log2-transformation settings
# - Median centering normalization
# - Statistical test configurations
# - Dimensional reduction algorithms
# - Timestamps for all processing steps
#
# USAGE EXAMPLES:
# ===============
# Load the processed case study data:
#   data(romics_proteins)
#
# Load the reproducible pipeline:
#   data(romics_pipeline)
#
# Apply pipeline to a new dataset:
#   new_result <- applyRomicsPipeline(new_unprocessed_object, romics_pipeline)
#
# Verify reproducibility:
#   identical(new_result$data, expected_processed_data)  # Should be TRUE
#
# Generate your own pipeline from a processed object:
#   my_pipeline <- createRomicsPipeline(my_processed_object)
#
# REFERENCES:
# ===========
# This case study demonstrates the principles and capabilities outlined in the
# RomicsProcessor manuscript. Complete methods are provided in the Supplementary
# Methods section of the accompanying publication.
