# Tests for combineRomicsObjects() with enhanced external object tracking
# Demonstrates complete analytical traceability for combined objects

# Helper function to create valid mock romics objects
make_mock_romics_for_combine <- function(uuid, name, steps, n_samples = 5, n_features = 10) {
  obj <- list(
    uuid = uuid,
    steps = steps,
    data = data.frame(matrix(rnorm(n_features * n_samples), nrow = n_features)),
    metadata = data.frame(matrix(1:(2 * n_samples), nrow = 2)),
    statistics = data.frame(),
    missingdata = data.frame(matrix(rnorm(n_features * n_samples), nrow = n_features)),
    original_data = data.frame(matrix(rnorm(n_features * n_samples), nrow = n_features)),
    original_metadata = data.frame(matrix(1:(2 * n_samples), nrow = 2)),
    custom_colors = c("#FF0000"),
    main_factor = "condition",
    omics_type = "test",
    omics_information = "test",
    dependencies = list()
  )

  colnames(obj$data) <- paste0("sample_", 1:n_samples)
  colnames(obj$metadata) <- paste0("sample_", 1:n_samples)
  colnames(obj$missingdata) <- paste0("sample_", 1:n_samples)
  colnames(obj$original_data) <- paste0("sample_", 1:n_samples)
  colnames(obj$original_metadata) <- paste0("sample_", 1:n_samples)
  rownames(obj$metadata) <- c("condition", "colors_romics")
  rownames(obj$original_metadata) <- c("condition", "colors_romics")

  class(obj) <- "romics_object"
  return(obj)
}

test_that("combineRomicsObjects records source object names", {
  # Create two related objects with same UUID and log status
  common_steps <- c("romics_object", "date|Jan_08_2025|createRomicsObject", "fun|createRomicsObject(...)",
                   "date|Jan_08_2025|log2transform", "fun|log2transform()")

  obj_1 <- make_mock_romics_for_combine("uuid-001", "obj_1", common_steps)
  obj_2 <- make_mock_romics_for_combine("uuid-001", "obj_2", c(common_steps, "date|Jan_08_2025|quantileNormSample", "fun|quantileNormSample()"))

  # Combine them
  combined <- combineRomicsObjects(obj_1, obj_2)

  # Check that source names are recorded
  external_source_entries <- combined$steps[grepl("^external_object_source\\|", combined$steps)]
  expect_length(external_source_entries, 1)
  expect_true(grepl("obj_2", external_source_entries[1]))
})

test_that("combineRomicsObjects records all source UUIDs", {
  # Create objects with consistent log status
  base_steps <- c("romics_object", "date|Jan_08|create", "fun|create()", "date|Jan_08|log2", "fun|log2transform()")
  obj_1 <- make_mock_romics_for_combine("uuid-001", "obj_1", base_steps)
  obj_2 <- make_mock_romics_for_combine("uuid-001", "obj_2", c(base_steps, "date|Jan_08|quant", "fun|quantileNormSample()"))
  obj_3 <- make_mock_romics_for_combine("uuid-001", "obj_3", c(base_steps, "date|Jan_08|med", "fun|medianCenterSample()"))

  # Combine three objects
  combined <- combineRomicsObjects(obj_1, obj_2, obj_3)

  # Check origin entries
  origin_entries <- combined$steps[grepl("^origin\\|", combined$steps)]
  expect_length(origin_entries, 3)
  expect_true(any(grepl("uuid-001", origin_entries)))
})

test_that("combineRomicsObjects records divergent steps from sources", {
  # Create objects with consistent log status
  base_steps <- c("romics_object", "date|Jan_08_2025|createRomicsObject", "fun|createRomicsObject(...)",
                  "date|Jan_08_2025|log2transform", "fun|log2transform()")
  obj_1 <- make_mock_romics_for_combine("uuid-001", "obj_1", base_steps)

  # Divergent object with additional quantileNormSample
  obj_2 <- make_mock_romics_for_combine(
    "uuid-001", "obj_2",
    c(base_steps, "date|Jan_08_2025|quantileNormSample", "fun|quantileNormSample()")
  )

  # Combine them
  combined <- combineRomicsObjects(obj_1, obj_2)

  # Check for external_object_fun entries
  external_fun_entries <- combined$steps[grepl("^fun\\|external_object_fun\\|", combined$steps)]
  # Should have at least the quantileNormSample entry
  expect_gt(length(external_fun_entries), 0)
  expect_true(any(grepl("quantileNormSample", external_fun_entries)))
})

test_that("combineRomicsObjects records combine operation", {
  # Create objects with consistent log status
  base_steps <- c("romics_object", "date|Jan|create", "fun|create()", "date|Jan|log2", "fun|log2transform()")
  obj_1 <- make_mock_romics_for_combine("uuid-001", "obj_1", base_steps)
  obj_2 <- make_mock_romics_for_combine("uuid-001", "obj_2", base_steps)

  initial_step_count <- length(obj_1$steps)

  # Combine
  combined <- combineRomicsObjects(obj_1, obj_2)

  # Should have more steps now
  expect_gt(length(combined$steps), initial_step_count)

  # Should record combine operation
  combine_steps <- combined$steps[grepl("combineRomicsObjects", combined$steps)]
  expect_length(combine_steps, 2)  # date and fun entries
})

test_that("combineRomicsObjects handles unrelated objects", {
  # Create unrelated objects (different UUIDs but same log status)
  base_steps <- c("romics_object", "date|Jan|create", "fun|create()", "date|Jan|log2", "fun|log2transform()")
  obj_1 <- make_mock_romics_for_combine("uuid-001", "obj_1", base_steps)
  obj_2 <- make_mock_romics_for_combine("uuid-002", "obj_2", base_steps)

  # Combine them
  combined <- combineRomicsObjects(obj_1, obj_2)

  # Should still record source object
  external_source_entries <- combined$steps[grepl("^external_object_source\\|", combined$steps)]
  expect_length(external_source_entries, 1)

  # Should record origin entries
  origin_entries <- combined$steps[grepl("^origin\\|", combined$steps)]
  expect_length(origin_entries, 2)
})

test_that("combineRomicsObjects handles three or more objects", {
  # Create three objects with consistent log status
  base_steps <- c("romics_object", "date|Jan|create", "fun|create()", "date|Jan|log2", "fun|log2transform()")
  obj_1 <- make_mock_romics_for_combine("uuid-001", "obj_1", base_steps)
  obj_2 <- make_mock_romics_for_combine("uuid-001", "obj_2", c(base_steps, "date|Jan|quant", "fun|quantileNormSample()"))
  obj_3 <- make_mock_romics_for_combine("uuid-001", "obj_3", c(base_steps, "date|Jan|med", "fun|medianCenterSample()"))

  # Combine all three
  combined <- combineRomicsObjects(obj_1, obj_2, obj_3)

  # Should have entries for both external objects
  external_source_entries <- combined$steps[grepl("^external_object_source\\|", combined$steps)]
  expect_length(external_source_entries, 2)

  # Should have 3 origin entries
  origin_entries <- combined$steps[grepl("^origin\\|", combined$steps)]
  expect_length(origin_entries, 3)
})

test_that("combineRomicsObjects records position index for external sources", {
  # Create objects with consistent log status
  base_steps <- c("romics_object", "date|Jan|create", "fun|create()", "date|Jan|log2", "fun|log2transform()")
  obj_1 <- make_mock_romics_for_combine("uuid-001", "obj_1", base_steps)
  obj_2 <- make_mock_romics_for_combine("uuid-001", "obj_2", base_steps)
  obj_3 <- make_mock_romics_for_combine("uuid-001", "obj_3", base_steps)

  # Combine
  combined <- combineRomicsObjects(obj_1, obj_2, obj_3)

  # Check position indices
  external_source_entries <- combined$steps[grepl("^external_object_source\\|", combined$steps)]
  expect_true(any(grepl("position=2", external_source_entries)))
  expect_true(any(grepl("position=3", external_source_entries)))
})

test_that("combineRomicsObjects preserves data integrity", {
  # Create objects with known data
  n_features <- 5
  n_samples <- 3
  data_1 <- data.frame(matrix(1:15, nrow = 5))
  data_2 <- data.frame(matrix(16:30, nrow = 5))

  obj_1 <- make_mock_romics_for_combine("uuid-001", "obj_1", c("romics_object", "date|Jan|create", "fun|create()"), n_samples = 3, n_features = 5)
  obj_1$data <- data_1
  colnames(obj_1$data) <- c("sample_1", "sample_2", "sample_3")

  obj_2 <- make_mock_romics_for_combine("uuid-001", "obj_2", c("romics_object", "date|Jan|create", "fun|create()"), n_samples = 3, n_features = 5)
  obj_2$data <- data_2
  colnames(obj_2$data) <- c("sample_1", "sample_2", "sample_3")

  # Combine
  combined <- combineRomicsObjects(obj_1, obj_2)

  # Should have more features (5 + 5)
  expect_equal(nrow(combined$data), 10)
  expect_equal(ncol(combined$data), 3)
})

test_that("combineRomicsObjects works with single object", {
  # Single object (edge case) with log status
  base_steps <- c("romics_object", "date|Jan|create", "fun|create()", "date|Jan|log2", "fun|log2transform()")
  obj_1 <- make_mock_romics_for_combine("uuid-001", "obj_1", base_steps)

  # Combine single object
  combined <- combineRomicsObjects(obj_1)

  # Should have origin entry
  origin_entries <- combined$steps[grepl("^origin\\|", combined$steps)]
  expect_length(origin_entries, 1)

  # Should have combine operation recorded
  combine_steps <- combined$steps[grepl("combineRomicsObjects", combined$steps)]
  expect_length(combine_steps, 2)
})
