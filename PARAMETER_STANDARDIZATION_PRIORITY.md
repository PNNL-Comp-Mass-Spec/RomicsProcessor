# RomicsProcessor Parameter Standardization Priority List

## Overview
This document prioritizes inconsistencies found in the parameter and function naming analysis. Priorities are based on:
- **Impact**: Number of functions affected, user confusion risk
- **Effort**: Complexity and number of files to modify
- **Breaking Change**: Whether this breaks existing user code
- **Dependencies**: Whether other fixes depend on this one

---

## TIER 1: CRITICAL (Fix First)
**Rationale**: These create the most user confusion and affect the most code

### 1.1 - CRITICAL: Standardize `plot_type` parameter naming
**Status**: BLOCKING other improvements  
**Functions Affected**: 5+ plotting functions
- `FeatureCorrelationHclustPlot()` - uses `plot_type = "dendrogram"`
- `SampleCorrelationHclustPlot()` - uses `plot_type = "dendrogram"`
- `romicsPCAplot()` - uses `plotType = "dual"` ❌ INCONSISTENT
- `multipleFeatureComparisonPlot()` - uses `plot_type = "scatter"`
- `bubblePlotFeatures()` - likely uses similar pattern

**Current State**:
```R
# Inconsistent naming
plot_type = "dendrogram"      # Functions 1 & 2
plotType = "dual"             # Function 3 (camelCase violation!)
```

**Issue**: `plotType` violates the `snake_case` convention used in 4 other similar functions

**Recommendation**: 
- Change `plotType` → `plot_type` in `romicsPCAplot()`, `romicsUmapPlot()`, `romicsTsnePlot()`
- Verify consistency across: `UmapDensityPlot`, `bubblePlotFeatures`, `multipleFeatureComparisonPlot`

**Breaking Change**: YES - users calling `plotType=` will need to update to `plot_type=`

**Effort**: LOW - search & replace + parameter rename in 3-5 functions

**Implementation Order**: 1st (highest priority)

---

### 1.2 - CRITICAL: Eliminate generic `type` parameter
**Status**: CAUSES SIGNIFICANT CONFUSION  
**Functions Affected**: 7 functions with 5+ different semantic meanings

**Current State**:
```R
# SAME parameter name, DIFFERENT meanings
singleFeaturePlot(type = "jb")                    # Plot type (boxplot/jitter)
multipleFeaturePlot(type = "scatter")             # Plot type (scatter/other)
romicsTrend(type = "both")                        # Analysis type (both/separate)
romicsSubset(type = "keep")                       # Filter mode (keep/remove)
romicsTransferEmbeddings(type = c("pca","umap")) # Method selection
romicsOutlierEval(type = "both")                  # Evaluation type
romicsFilterFeature(type = "both")                # Filter type
```

**Issue**: Users cannot guess what `type=` means without reading documentation. Different functions expect different values.

**Recommendation - Rename to specific parameter names**:
| Current | Should Be | Function |
|---------|-----------|----------|
| `type = "jb"` | `plot_type = "jb"` | `singleFeaturePlot` |
| `type = "scatter"` | `plot_type = "scatter"` | `multipleFeaturePlot` |
| `type = "both"` | `analysis_type = "both"` | `romicsTrend` |
| `type = "keep"` | `filter_mode = "keep"` | `romicsSubset` |
| `type = c(...)` | `methods = c(...)` | `romicsTransferEmbeddings` |
| `type = "both"` | `eval_type = "both"` | `romicsOutlierEval` |
| `type = "both"` | `filter_type = "both"` | `romicsFilterFeature` |

**Breaking Change**: YES - all 7 functions change their parameter

**Effort**: MEDIUM - careful renaming in 7 functions, update documentation

**Dependencies**: Overlaps with 1.1 (some functions affected by both)

**Implementation Order**: 2nd (depends on 1.1 for overlap resolution)

---

### 1.3 - CRITICAL: Standardize axis limit parameters
**Status**: INCONSISTENT NAMING  
**Functions Affected**: 2-3 functions + future plotting functions

**Current State**:
```R
singleFeaturePlot(limits = NULL)           # y-axis limits only
romicsVolcano(xlim = NULL, ylim = NULL)    # Both axes
# New additions should use xlim/ylim pattern
```

**Issue**: Users expect `xlim`/`ylim` (standard R/ggplot2 convention) but `singleFeaturePlot` uses non-standard `limits`

**Recommendation**:
- Update `singleFeaturePlot()`: Change `limits` → `ylim` (or add both `xlim` and `ylim`)
- Document: All future plotting functions should use `xlim`, `ylim`, `zlim`
- Update documentation to mention R standard convention

**Breaking Change**: YES - `singleFeaturePlot(limits=...)` becomes `ylim=...`

**Effort**: LOW - rename 1 parameter in 1 function, add documentation

**Implementation Order**: 3rd

---

## TIER 2: HIGH PRIORITY
**Rationale**: Important for code consistency, affects many functions

### 2.1 - HIGH: Standardize boolean parameter naming convention
**Status**: 22 VARIATIONS EXIST  
**Functions Affected**: 20+ functions

**Current Boolean Parameters**:
```R
label = FALSE/TRUE           # 3 functions with DIFFERENT defaults
verbose = TRUE               # 4 functions (consistent)
show_quantiles = TRUE        # 2 functions
percent = TRUE               # 2 functions
padj = TRUE                  # Statistical functions
log_factor = FALSE           # 3 functions
cont.rm = TRUE               # Data cleaning (dot notation!)
site.rm = TRUE               # Data cleaning (dot notation!)
rev.rm = TRUE                # Data cleaning (dot notation!)
scale = TRUE/FALSE           # 5 functions (different defaults)
```

**Issues**:
1. **Prefix inconsistency**: No clear pattern (`label` vs `show_quantiles` vs `verbose`)
2. **Default inconsistency**: `label=TRUE` in one function, `label=FALSE` in another
3. **Dot notation**: `cont.rm`, `site.rm` violate camelCase pattern
4. **Unclear intent**: `scale` could mean many things

**Recommendation - Establish pattern**:
- **For display/output flags**: Use `show_` prefix → `show_labels`, `show_quantiles`, `show_verbose`
- **For filter/removal flags**: Use descriptive names → `remove_contaminants`, `remove_site_only` (not `cont.rm`)
- **For transformation flags**: Use `apply_` prefix → `apply_scaling`, `apply_log_transform`
- **Standardize defaults**: Document why each default (usually `TRUE` for important transformations)

**Breaking Change**: YES - parameter names change in 20+ functions

**Effort**: MEDIUM-HIGH - systematic rename across package, update documentation

**Dependencies**: Should align with 1.2 (some overlap in functions)

**Implementation Order**: 4th

---

### 2.2 - HIGH: Standardize statistical type parameters
**Status**: 4 DIFFERENT PARAMETER NAMES  
**Functions Affected**: 15+ statistical functions

**Current State**:
```R
stat_type = "auto"           # romicsVolcano, romicsVolcanoByCluster
p_type = "p"                 # Pvalue type selector
corr_type = "pearson"        # Correlation type
quantification_type = "LFQ"  # Data extraction type
```

**Issue**: No consistent naming pattern; makes it hard for users to learn

**Recommendation - Use consistent `_type` suffix**:
- `stat_type` (already good) ✓ - statistical test type
- `p_value_type` or `p_adjustment_type` (rename `p_type`) - for p-value selection
- `correlation_type` (rename `corr_type`) - more explicit
- `quantification_type` (keep as is) ✓ - already explicit

**Breaking Change**: YES - 2 functions change parameter names

**Effort**: LOW - rename 2-3 parameters

**Implementation Order**: 5th

---

### 2.3 - HIGH: Consolidate factor parameter naming
**Status**: INCONSISTENT ACROSS 45+ FUNCTIONS  
**Functions Affected**: 45+ functions

**Current State**:
```R
# Two naming patterns in use
factor = "main"              # ~30 functions (STANDARD)
main_factor = "main"         # ~15 functions (NON-STANDARD)
main_factor = "none"         # Some functions with different default
```

**Issue**: Users get confused whether to use `factor=` or `main_factor=`

**Recommendation**:
- Establish `factor = "main"` as the standard
- Phase out `main_factor` in new documentation
- Update old functions to use `factor` where redundant
- Keep `main_factor` in `createRomicsObject()` and core functions (part of object structure)

**Breaking Change**: YES - but can be gradual with deprecation warnings

**Effort**: MEDIUM-HIGH - 15-20 function updates

**Implementation Order**: 6th (lower priority than Tier 1)

---

## TIER 3: MEDIUM PRIORITY
**Rationale**: Nice to fix, but lower impact

### 3.1 - MEDIUM: Remove duplicate `extractScils` function definition
**Status**: DUPLICATE DEFINITION  
**Functions Affected**: 1 function (defined twice)

**Current State**:
```R
extractScils()  # Defined in: 21_new_scils_import.R
extractScils()  # Also defined in: 01_OmicsData_Import.R
```

**Issue**: Two definitions create confusion; may cause namespace conflicts

**Recommendation**:
- Determine which is the "correct" version (check dates, functionality)
- Keep one definition, add deprecation warning to old version
- Update imports/documentation

**Breaking Change**: NO (internal change)

**Effort**: LOW - merge two files, test

**Implementation Order**: 7th (low priority, internal cleanup)

---

### 3.2 - MEDIUM: Add missing `Zcomp` parameter to `romicsUmapDensityPlot`
**Status**: INCONSISTENT WITH SIMILAR FUNCTIONS  
**Functions Affected**: 1 function

**Current State**:
```R
# These have Zcomp:
romicsPCAplot(Xcomp=1, Ycomp=2, Zcomp=NULL, ...)
romicsUmapPlot(Xcomp=1, Ycomp=2, Zcomp=NULL, ...)
romicsTsnePlot(Xcomp=1, Ycomp=2, Zcomp=NULL, ...)

# This is MISSING Zcomp:
romicsUmapDensityPlot(Xcomp=1, Ycomp=2, ...)  # ❌ NO Zcomp
```

**Issue**: Inconsistent API between similar functions; user might expect 3D support

**Recommendation**:
- Add `Zcomp = NULL` parameter to `romicsUmapDensityPlot()`
- Clarify if 3D density plots are supported or if parameter is placeholder

**Breaking Change**: NO (additive change)

**Effort**: LOW - add parameter + documentation

**Implementation Order**: 8th (low priority, additive change)

---

## TIER 4: POLISH (Nice to Have)
**Rationale**: Code quality improvements, not critical for functionality

### 4.1 - POLISH: Standardize parameter spacing convention
**Status**: COSMETIC  
**Functions Affected**: All functions (inconsistently)

**Current State**: Some functions use `param=value`, others use `param = value`

**Recommendation**: Adopt one style (suggest `param = value` with spaces for readability)

**Effort**: LOW but tedious - run R code formatter

---

### 4.2 - POLISH: Document or remove `_old` functions
**Status**: OLD FUNCTION VERSIONS  
**Functions Affected**: 3 functions ending in `_old`
- `romicsTtest_old()`
- `romicsWilcoxTest_old()`
- `romicsGlmBinomial_old()`

**Recommendation**: Add deprecation notices or remove entirely

**Effort**: LOW - add warnings or delete

---

## IMPLEMENTATION ROADMAP

### Phase 1: Critical Fixes (Weeks 1-2)
1. **1.1** - Standardize `plot_type` (rename `plotType`)
2. **1.2** - Eliminate generic `type` parameter (rename to specific names)
3. **1.3** - Standardize axis limits (`limits` → `ylim`)

**Expected breaking changes**: 7-10 functions  
**User impact**: MEDIUM (easy to fix in user code)

---

### Phase 2: High Priority (Weeks 3-4)
4. **2.1** - Standardize boolean naming (systematic rename)
5. **2.2** - Standardize statistical type parameters
6. **2.3** - Consolidate factor parameter naming

**Expected breaking changes**: 35-40 functions  
**User impact**: MEDIUM-HIGH (affects many function calls)

---

### Phase 3: Medium Priority (Week 5)
7. **3.1** - Remove duplicate `extractScils`
8. **3.2** - Add `Zcomp` to `romicsUmapDensityPlot`

**Expected breaking changes**: 0 (internal change + additive)  
**User impact**: LOW

---

### Phase 4: Polish (Week 6)
9. **4.1** - Standardize parameter spacing
10. **4.2** - Document/remove `_old` functions

**Breaking changes**: 0  
**User impact**: NONE

---

## SUMMARY TABLE

| Priority | Fix | Functions | Effort | Breaking | Phase |
|----------|-----|-----------|--------|----------|-------|
| 1.1 | `plotType` → `plot_type` | 5 | LOW | YES | 1 |
| 1.2 | Remove generic `type` | 7 | MEDIUM | YES | 1 |
| 1.3 | `limits` → `ylim` | 2 | LOW | YES | 1 |
| 2.1 | Boolean naming pattern | 20+ | MEDIUM-HIGH | YES | 2 |
| 2.2 | Statistical type params | 2 | LOW | YES | 2 |
| 2.3 | Factor parameter naming | 15-20 | MEDIUM-HIGH | YES | 2 |
| 3.1 | Remove `extractScils` duplicate | 1 | LOW | NO | 3 |
| 3.2 | Add `Zcomp` parameter | 1 | LOW | NO | 3 |
| 4.1 | Parameter spacing | all | LOW | NO | 4 |
| 4.2 | Document `_old` functions | 3 | LOW | NO | 4 |

---

## MIGRATION GUIDE FOR USERS (Draft)

When updating to the next version of RomicsProcessor, use this mapping:

### Phase 1 Changes
```R
# OLD → NEW
plotType = "value"           → plot_type = "value"
type = "jb"                  → plot_type = "jb" (singleFeaturePlot)
type = "both"                → analysis_type = "both" (romicsTrend)
singleFeaturePlot(limits=c(0,100))  → singleFeaturePlot(ylim=c(0,100))
```

### Phase 2 Changes
```R
# OLD → NEW
verbose = TRUE               → show_verbose = TRUE
label = FALSE                → show_labels = FALSE
cont.rm = TRUE               → remove_contaminants = TRUE
stat_type = "auto"           → stat_type = "auto" (NO CHANGE)
main_factor = "cond"         → factor = "cond"
```

---

## RECOMMENDATION FOR DISCUSSION

**Which should we tackle first?**

**Option A: Quick Win (2-3 weeks)**
- Focus only on **Tier 1** (1.1, 1.2, 1.3)
- Get critical inconsistencies fixed
- Plan Tier 2 separately

**Option B: Comprehensive Fix (5-6 weeks)**
- Do all of Tier 1 + Tier 2 together
- One breaking version with all changes at once
- Cleaner for users (one migration, not multiple)

**Option C: Incremental (8+ weeks)**
- Phase 1 → 2 → 3 → 4
- Manage breaking changes gradually
- More user-friendly rollout

**Recommendation**: Option B - do Tier 1 + Tier 2 together as one major version bump
