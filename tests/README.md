# RomicsProcessor Tests

This directory contains unit tests for the RomicsProcessor package using the `testthat` framework.

## Test Organization

- `testthat.R` - Main test runner configuration
- `testthat/test_romics_object.R` - Tests for core romics_object functionality
- `testthat/test_plotting.R` - Tests for plotting functions

## Running Tests

### Run all tests
```r
devtools::test()
```

### Run specific test file
```r
devtools::test_file("tests/testthat/test_romics_object.R")
```

### Run with coverage report
```r
devtools::test_coverage()
```

## Test Categories

### Core Object Tests (`test_romics_object.R`)
- `is.romicsObject()` - Validates romics object recognition
- `createRomicsObject()` - Tests object creation and structure
- `romicsExtractFactor()` - Tests factor extraction
- `romicsFactorNames()` - Tests factor name retrieval

### Plotting Tests (`test_plotting.R`)
- `romicsFactorFrequencyBarplot()` - Tests frequency barplot generation
- Custom color handling
- Error handling for missing factors

## Adding New Tests

When adding new tests:

1. Create a new file `test_*.R` in `testthat/` directory
2. Use descriptive test names following `test_that("description", { ... })`
3. Include both success cases and error cases
4. Use `expect_*` functions from testthat (e.g., `expect_true()`, `expect_equal()`, `expect_error()`)
5. Run tests locally before committing

## Test Coverage Goals

- Core functions: 70%+ coverage
- Public API functions: 60%+ coverage
- Utility functions: 50%+ coverage
