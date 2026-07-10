# RomicsProcessor News

## Version 2.0.0 (2026-07-08)

### Major Release - Complete Object Encapsulation & Spatial Omics Support

RomicsProcessor v2.0 represents a milestone release with complete capabilities for encapsulating omics datasets and their entire processing history within the romics_object structure, enabling complete analytical traceability for FAIR-compliant science. This version also introduces comprehensive support for spatial omics datasets.

#### Complete Object Encapsulation & History Tracking
- **Analytical Traceability** (`checkRelationRomicsObjects`): Determine relationships between romics_objects via UUID and processing history comparison
  - Detects object lineage: identical objects, linear evolution, true branching, or unrelated origins
  - Identifies branching points and divergent processing steps across analytical branches
  - Enables complete analytical traceability for FAIR science (Findable, Accessible, Interoperable, Reusable)
  - romics_object now fully encapsulates dataset and complete processing history for reproducible analysis

- **Enhanced Embedding Transfer Tracking** (`romicsTransferEmbeddings`): Complete provenance tracking for transferred embeddings
  - Automatically captures origin object's UUID and name
  - Records all divergent steps from source object with "external_object_" prefix
  - Documents which specific transformations generated the embeddings
  - Stores origin metadata in steps layer for complete analytical chain
  - Enables full traceability of embedding derivation across branching analyses

- **Enhanced Object Combining Tracking** (`combineRomicsObjects`): Complete traceability for combined objects
  - Records all source object names and UUIDs in steps layer
  - Documents divergent processing steps from each combined object with "external_object_" prefix
  - Adds position indicators to track which object contributed what (position=1, position=2, etc.)
  - Stores origin entries for all source objects enabling complete lineage traceability
  - Enables users to understand how combined data was derived from individual sources

#### Spatial Omics Support
- **Enhanced Spatial Integration**: Complete workflow support for spatial omics data analysis
  - Format conversion with SpatialExperiment objects for spatial analysis frameworks
  - Embedding preservation for spatial coordinates and dimensional reductions
  - Interactive spatial visualization with improved mapping capabilities
  - Support for multiple spatial omics platforms (imaging mass spectrometry, spatial transcriptomics, etc.)

#### Enhanced Processing Capabilities
- **Reference Region-Based Batch Correction**: Extended `romicsBatchCorrection()` with reference region support
  - New `method = "reference"` option for reference-guided batch correction
  - `reference_filter` parameter to identify reference regions
  - `use_ref_params_only` option to learn parameters exclusively from reference region

- **Custom Feature Labeling in Volcano Plots**: Enhanced `romicsVolcano()` and `romicsVolcanoByCluster()`
  - New `label_features_list` parameter for user-specified feature labeling

- **Factor Frequency Barplot** (`romicsFactorFrequencyBarplot`): Streamlined visualization for categorical distributions

### Backward Compatibility
- ✅ All new features are fully backward compatible
- ✅ Existing code continues to work without modifications
- ✅ Default parameter values preserve original behavior

### General Improvements
- Comprehensive test suite with 41+ tests demonstrating analytical traceability and branching scenarios
- Many bugs from previous versions have been fixed
- Enhanced documentation and roxygen updates
- Repository cleanup and improved configuration

---