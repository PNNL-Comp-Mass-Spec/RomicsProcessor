# Tests for checkRelationRomicsObjects() function
# Demonstrates analytical traceability through UUID and processing history

test_that("checkRelationRomicsObjects function exists", {
  # Verify the function is available
  expect_true(exists("checkRelationRomicsObjects"))
  expect_true(is.function(checkRelationRomicsObjects))
})

test_that("checkRelationRomicsObjects validates inputs", {
  # Test that function validates romics objects
  not_romics <- list(data = matrix(1:100), metadata = data.frame())

  # Should error on invalid input
  expect_error(
    checkRelationRomicsObjects(not_romics, not_romics)
  )
})

test_that("checkRelationRomicsObjects returns proper structure", {
  # Helper to create minimal valid romics-like object
  make_mock_romics <- function(uuid, steps) {
    obj <- list(
      uuid = uuid,
      steps = steps,
      data = data.frame(matrix(rnorm(100), nrow = 10)),
      metadata = data.frame(matrix(1:100, nrow = 10)),
      statistics = data.frame(),
      missingdata = data.frame(),
      original_data = data.frame(matrix(rnorm(100), nrow = 10)),
      original_metadata = data.frame(matrix(1:100, nrow = 10)),
      custom_colors = c("#FF0000"),
      main_factor = c("A", "B"),
      omics_type = "test",
      omics_information = "test",
      dependencies = list(),
      steps = steps
    )
    class(obj) <- "romics_object"
    return(obj)
  }

  mock_romics1 <- make_mock_romics("test-uuid-001", c("romics_object", "step1", "step2"))
  mock_romics2 <- make_mock_romics("test-uuid-001", c("romics_object", "step1", "step2", "step3"))

  # Run comparison without verbose output
  result <- checkRelationRomicsObjects(mock_romics1, mock_romics2, verbose = FALSE)

  # Check result structure
  expect_true(is.list(result))
  expect_true("related" %in% names(result))
  expect_true("uuid1" %in% names(result))
  expect_true("uuid2" %in% names(result))
  expect_true("common_steps" %in% names(result))
  expect_true("branching_point" %in% names(result))
  expect_true("divergent_steps_obj1" %in% names(result))
  expect_true("divergent_steps_obj2" %in% names(result))
  expect_true("summary" %in% names(result))
})

test_that("checkRelationRomicsObjects detects identical objects", {
  # Helper function
  make_mock_romics <- function(uuid, steps) {
    obj <- list(
      uuid = uuid,
      steps = steps,
      data = data.frame(matrix(rnorm(100), nrow = 10)),
      metadata = data.frame(matrix(1:100, nrow = 10)),
      statistics = data.frame(),
      missingdata = data.frame(),
      original_data = data.frame(matrix(rnorm(100), nrow = 10)),
      original_metadata = data.frame(matrix(1:100, nrow = 10)),
      custom_colors = c("#FF0000"),
      main_factor = c("A", "B"),
      omics_type = "test",
      omics_information = "test",      dependencies = list()
    )
    class(obj) <- "romics_object"
    return(obj)
  }

  mock_romics <- make_mock_romics("same-uuid", c("romics_object", "step1", "step2"))
  result <- checkRelationRomicsObjects(mock_romics, mock_romics, verbose = FALSE)

  # Should detect as related
  expect_true(result$related)
  # Should have identical UUIDs
  expect_equal(result$uuid1, result$uuid2)
  # Should have no divergent steps
  expect_equal(length(result$divergent_steps_obj1), 0)
  expect_equal(length(result$divergent_steps_obj2), 0)
})

test_that("checkRelationRomicsObjects detects linear evolution (predecessor)", {
  # Helper function
  make_mock_romics <- function(uuid, steps) {
    obj <- list(
      uuid = uuid,
      steps = steps,
      data = data.frame(matrix(rnorm(100), nrow = 10)),
      metadata = data.frame(matrix(1:100, nrow = 10)),
      statistics = data.frame(),
      missingdata = data.frame(),
      original_data = data.frame(matrix(rnorm(100), nrow = 10)),
      original_metadata = data.frame(matrix(1:100, nrow = 10)),
      custom_colors = c("#FF0000"),
      main_factor = c("A", "B"),
      omics_type = "test",
      omics_information = "test",      dependencies = list()
    )
    class(obj) <- "romics_object"
    return(obj)
  }

  # Object 1: original
  mock_romics1 <- make_mock_romics("same-uuid", c("romics_object", "step1", "step2"))

  # Object 2: derivative with additional steps
  mock_romics2 <- make_mock_romics("same-uuid", c("romics_object", "step1", "step2", "step3", "step4"))

  result <- checkRelationRomicsObjects(mock_romics1, mock_romics2, verbose = FALSE)

  # Should detect as related
  expect_true(result$related)
  # Should have common initial steps
  expect_equal(length(result$common_steps), 3)
  # Object 2 should have divergent steps
  expect_equal(length(result$divergent_steps_obj2), 2)
  # Object 1 should have no divergent steps (it's a predecessor)
  expect_equal(length(result$divergent_steps_obj1), 0)
})

test_that("checkRelationRomicsObjects detects branching (two derivatives)", {
  # Helper function
  make_mock_romics <- function(uuid, steps) {
    obj <- list(
      uuid = uuid,
      steps = steps,
      data = data.frame(matrix(rnorm(100), nrow = 10)),
      metadata = data.frame(matrix(1:100, nrow = 10)),
      statistics = data.frame(),
      missingdata = data.frame(),
      original_data = data.frame(matrix(rnorm(100), nrow = 10)),
      original_metadata = data.frame(matrix(1:100, nrow = 10)),
      custom_colors = c("#FF0000"),
      main_factor = c("A", "B"),
      omics_type = "test",
      omics_information = "test",      dependencies = list()
    )
    class(obj) <- "romics_object"
    return(obj)
  }

  # Common ancestor
  common_steps <- c("romics_object", "step1", "step2")

  # Branch 1: Log2 transform then Quantile norm
  mock_romics1 <- make_mock_romics("same-uuid", c(common_steps, "romicsLog2transform", "romicsQuantileNorm"))

  # Branch 2: Log2 transform then Median norm
  mock_romics2 <- make_mock_romics("same-uuid", c(common_steps, "romicsLog2transform", "romicsMedianNorm"))

  result <- checkRelationRomicsObjects(mock_romics1, mock_romics2, verbose = FALSE)

  # Should detect as related
  expect_true(result$related)
  # Should have common steps up to log2 transform
  expect_equal(length(result$common_steps), 4)
  # Both should have divergent steps
  expect_equal(length(result$divergent_steps_obj1), 1)
  expect_equal(length(result$divergent_steps_obj2), 1)
  # Divergent steps should be different (Quantile vs Median norm)
  expect_false(identical(result$divergent_steps_obj1, result$divergent_steps_obj2))
})

test_that("checkRelationRomicsObjects detects unrelated objects", {
  # Helper function
  make_mock_romics <- function(uuid, steps) {
    obj <- list(
      uuid = uuid,
      steps = steps,
      data = data.frame(matrix(rnorm(100), nrow = 10)),
      metadata = data.frame(matrix(1:100, nrow = 10)),
      statistics = data.frame(),
      missingdata = data.frame(),
      original_data = data.frame(matrix(rnorm(100), nrow = 10)),
      original_metadata = data.frame(matrix(1:100, nrow = 10)),
      custom_colors = c("#FF0000"),
      main_factor = c("A", "B"),
      omics_type = "test",
      omics_information = "test",      dependencies = list()
    )
    class(obj) <- "romics_object"
    return(obj)
  }

  # Object 1
  mock_romics1 <- make_mock_romics("uuid-exp-001", c("romics_object", "step1"))

  # Object 2 with different UUID
  mock_romics2 <- make_mock_romics("uuid-exp-002", c("romics_object", "step1"))

  result <- checkRelationRomicsObjects(mock_romics1, mock_romics2, verbose = FALSE)

  # Should NOT be related
  expect_false(result$related)
  # Summary should indicate different origins
  expect_true(grepl("NOT related", result$summary, ignore.case = TRUE))
})
