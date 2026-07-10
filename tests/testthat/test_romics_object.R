# Tests for romics_object creation and basic functions
# Note: Full integration tests with createRomicsObject are complex due to data structure requirements
# These tests verify that core functionality and utilities are available

test_that("ROP_colors is available and valid", {
  # Check that the package has the default color palette
  expect_true(exists("ROP_colors"))
  expect_true(is.character(ROP_colors))
  expect_true(length(ROP_colors) > 0)
})

test_that("createRomicsObject function exists", {
  # Check that the main function is exported
  expect_true(exists("createRomicsObject"))
  expect_true(is.function(createRomicsObject))
})

test_that("romicsExtractFactor function exists", {
  # Check that factor extraction function is available
  expect_true(exists("romicsExtractFactor"))
  expect_true(is.function(romicsExtractFactor))
})

test_that("romicsFactorNames function exists", {
  # Check that factor names function is available
  expect_true(exists("romicsFactorNames"))
  expect_true(is.function(romicsFactorNames))
})
