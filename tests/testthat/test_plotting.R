# Tests for plotting functions

test_that("theme_ROP exists and is a function", {
  # Check that the package has the ROP theme
  expect_true(exists("theme_ROP"))
  expect_true(is.function(theme_ROP))
})

test_that("romicsFactorFrequencyBarplot function is exported", {
  # Check that the function exists in the package namespace
  expect_true(exists("romicsFactorFrequencyBarplot"))
  expect_true(is.function(romicsFactorFrequencyBarplot))
})

test_that("ggplot2 package is available", {
  # Verify dependencies
  expect_true(requireNamespace("ggplot2", quietly = TRUE))
})
