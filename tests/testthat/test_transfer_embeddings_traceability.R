# Tests for romicsTransferEmbeddings() with enhanced step tracking
# Demonstrates complete analytical traceability for embedding transfers

# Helper function to create valid mock romics objects with embeddings
make_mock_romics_with_embeddings <- function(uuid, steps, embedding_types = c("pca")) {
  n_samples <- 5
  n_features <- 10

  obj <- list(
    uuid = uuid,
    steps = steps,
    data = data.frame(matrix(rnorm(n_features * n_samples), nrow = n_features)),
    metadata = data.frame(matrix(1:(2 * n_samples), nrow = 2)),
    statistics = data.frame(),
    missingdata = data.frame(),
    original_data = data.frame(matrix(rnorm(n_features * n_samples), nrow = n_features)),
    original_metadata = data.frame(matrix(1:(2 * n_samples), nrow = 2)),
    custom_colors = c("#FF0000"),
    main_factor = c("A", "B", "A", "B", "A"),
    omics_type = "test",
    omics_information = "test",
    dependencies = list(),
    embeddings = data.frame(matrix(ncol = n_samples, nrow = 0))
  )

  colnames(obj$data) <- paste0("sample_", 1:n_samples)
  colnames(obj$metadata) <- paste0("sample_", 1:n_samples)
  colnames(obj$embeddings) <- paste0("sample_", 1:n_samples)
  rownames(obj$metadata) <- c("factor", "colors_romics")

  # Add embeddings based on types requested
  if("pca" %in% embedding_types) {
    pca_components <- data.frame(matrix(rnorm(3 * n_samples), nrow = 3))
    rownames(pca_components) <- c("pca_component_1", "pca_component_2", "pca_component_3")
    colnames(pca_components) <- paste0("sample_", 1:n_samples)
    obj$embeddings <- rbind(obj$embeddings, pca_components)
  }

  if("umap" %in% embedding_types) {
    umap_components <- data.frame(matrix(rnorm(2 * n_samples), nrow = 2))
    rownames(umap_components) <- c("umap_component_1", "umap_component_2")
    colnames(umap_components) <- paste0("sample_", 1:n_samples)
    obj$embeddings <- rbind(obj$embeddings, umap_components)
  }

  class(obj) <- "romics_object"
  return(obj)
}

test_that("romicsTransferEmbeddings records origin metadata", {
  # Create base object with PCA embeddings
  obj_a <- make_mock_romics_with_embeddings(
    "uuid-a-001",
    c("romics_object", "date|Jan_08_2025|romicsPCA", "fun|romicsPCA(ncp=5)"),
    embedding_types = c("pca")
  )

  # Create two divergent branches
  obj_b1 <- obj_a
  obj_b1$uuid <- "uuid-a-001"  # Same UUID - branching from obj_a
  obj_b1$steps <- c(obj_a$steps, "date|Jan_08_2025|log2transform", "fun|log2transform()")

  obj_b2 <- obj_a
  obj_b2$uuid <- "uuid-a-001"  # Same UUID - branching from obj_a
  obj_b2$steps <- c(obj_a$steps, "date|Jan_08_2025|quantileNormSample", "fun|quantileNormSample()")

  # Transfer embeddings from b1 to b2
  obj_b2_transferred <- romicsTransferEmbeddings(obj_b1, obj_b2, type = "pca")

  # Check that origin metadata was added
  origin_steps <- obj_b2_transferred$steps[grepl("^origin\\|", obj_b2_transferred$steps)]
  expect_length(origin_steps, 1)
  expect_true(grepl("uuid-a-001", origin_steps[1]))
})

test_that("romicsTransferEmbeddings records divergent steps from origin", {
  # Create base object
  obj_a <- make_mock_romics_with_embeddings(
    "uuid-a-001",
    c("romics_object", "date|Jan_08_2025_14_00_00|romicsPCA", "fun|romicsPCA(ncp=5)"),
    embedding_types = c("pca")
  )

  # Create divergent branch with log transform
  obj_b1 <- obj_a
  obj_b1$uuid <- "uuid-a-001"
  obj_b1$steps <- c(
    obj_a$steps,
    "date|Jan_08_2025_14_10_00|log2transform",
    "fun|log2transform()"
  )

  # Create receiving object without divergence
  obj_b2 <- obj_a
  obj_b2$uuid <- "uuid-a-001"
  obj_b2$steps <- c(obj_a$steps)

  # Transfer embeddings
  obj_b2_transferred <- romicsTransferEmbeddings(obj_b1, obj_b2, type = "pca")

  # Check for external_object steps
  external_steps <- obj_b2_transferred$steps[grepl("^fun\\|external_object_fun\\|", obj_b2_transferred$steps)]
  expect_length(external_steps, 1)
  expect_true(grepl("log2transform", external_steps[1]))
})

test_that("romicsTransferEmbeddings records transfer operation", {
  # Create objects
  obj_a <- make_mock_romics_with_embeddings(
    "uuid-a-001",
    c("romics_object", "date|Jan_08_2025|romicsPCA", "fun|romicsPCA(ncp=5)"),
    embedding_types = c("pca")
  )

  obj_b1 <- obj_a
  obj_b1$uuid <- "uuid-a-001"
  obj_b1$steps <- c(obj_a$steps, "date|Jan_08_2025|log2transform", "fun|log2transform()")

  obj_b2 <- obj_a
  obj_b2$uuid <- "uuid-a-001"
  obj_b2$steps <- c(obj_a$steps)

  initial_step_count <- length(obj_b2$steps)

  # Transfer embeddings
  obj_b2_transferred <- romicsTransferEmbeddings(obj_b1, obj_b2, type = "pca")

  # Should have added: external step(s), transfer date, transfer fun, origin metadata
  expect_gt(length(obj_b2_transferred$steps), initial_step_count)

  # Should have romicsTransferEmbeddings recorded
  transfer_steps <- obj_b2_transferred$steps[grepl("romicsTransferEmbeddings", obj_b2_transferred$steps)]
  expect_length(transfer_steps, 2)  # date and fun entries
})

test_that("romicsTransferEmbeddings works with unrelated objects", {
  # Create two completely unrelated objects with different UUIDs
  obj_a <- make_mock_romics_with_embeddings(
    "uuid-a-001",
    c("romics_object", "date|Jan_08_2025|romicsPCA", "fun|romicsPCA(ncp=5)"),
    embedding_types = c("pca")
  )

  obj_b <- make_mock_romics_with_embeddings(
    "uuid-b-002",  # Different UUID
    c("romics_object", "date|Jan_08_2025|createRomicsObject", "fun|createRomicsObject()"),
    embedding_types = c()
  )

  # Transfer should still work, but with no external_object steps
  obj_b_transferred <- romicsTransferEmbeddings(obj_a, obj_b, type = "pca")

  # Check that origin metadata is still recorded
  origin_steps <- obj_b_transferred$steps[grepl("^origin\\|", obj_b_transferred$steps)]
  expect_length(origin_steps, 1)
  expect_true(grepl("uuid-a-001", origin_steps[1]))

  # Should have transfer recorded
  transfer_steps <- obj_b_transferred$steps[grepl("romicsTransferEmbeddings", obj_b_transferred$steps)]
  expect_length(transfer_steps, 2)
})

test_that("romicsTransferEmbeddings preserves embedding data", {
  # Create base object with known PCA values
  obj_a <- make_mock_romics_with_embeddings(
    "uuid-a-001",
    c("romics_object", "date|Jan_08_2025|romicsPCA", "fun|romicsPCA(ncp=5)"),
    embedding_types = c("pca")
  )

  # Store original embedding values
  original_pca <- obj_a$embeddings[grepl("pca_component_", rownames(obj_a$embeddings)), ]

  obj_b <- obj_a
  obj_b$uuid <- "uuid-a-001"
  obj_b$embeddings <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(obj_b$embeddings) <- paste0("sample_", 1:5)

  # Transfer embeddings
  obj_b_transferred <- romicsTransferEmbeddings(obj_a, obj_b, type = "pca")

  # Check that embeddings were transferred correctly
  transferred_pca <- obj_b_transferred$embeddings[grepl("pca_component_", rownames(obj_b_transferred$embeddings)), ]
  expect_equal(nrow(transferred_pca), 3)
  expect_equal(as.matrix(transferred_pca), as.matrix(original_pca))
})

test_that("romicsTransferEmbeddings handles multiple embedding types", {
  # Create objects with multiple embeddings
  obj_a <- make_mock_romics_with_embeddings(
    "uuid-a-001",
    c("romics_object", "date|Jan_08_2025|romicsPCA", "fun|romicsPCA(ncp=5)"),
    embedding_types = c("pca", "umap")
  )

  obj_b1 <- obj_a
  obj_b1$uuid <- "uuid-a-001"
  obj_b1$steps <- c(obj_a$steps, "date|Jan_08_2025|log2transform", "fun|log2transform()")

  obj_b2 <- obj_a
  obj_b2$uuid <- "uuid-a-001"
  obj_b2$embeddings <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(obj_b2$embeddings) <- paste0("sample_", 1:5)

  # Transfer both types
  obj_b2_transferred <- romicsTransferEmbeddings(obj_b1, obj_b2, type = c("pca", "umap"))

  # Check both were transferred
  pca_rows <- sum(grepl("pca_component_", rownames(obj_b2_transferred$embeddings)))
  umap_rows <- sum(grepl("umap_component_", rownames(obj_b2_transferred$embeddings)))

  expect_equal(pca_rows, 3)
  expect_equal(umap_rows, 2)
})
