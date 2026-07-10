#' checkRelationRomicsObjects()
#' @description Determine if two romics_objects are derived from one another by comparing their UUIDs and processing steps. Identifies the branching point and all subsequent divergent steps between two objects.
#' @param romics_object1 First romics_object to compare
#' @param romics_object2 Second romics_object to compare
#' @param verbose Logical (TRUE or FALSE). If TRUE, prints detailed information about the relationship between the two objects. Default: TRUE
#' @details This function implements the Findability principle by tracking analytical traceability through UUIDs and processing history.
#' It first checks if both objects share the same UUID (indicating they come from the same origin).
#' If UUIDs match, it compares the steps layer to find where the two objects diverge and what subsequent steps differ between them.
#' This enables researchers to understand how derivative datasets branch from an original object and what specific operations led to those differences.
#' @return A list containing:
#' \itemize{
#'   \item \code{related}: Logical. TRUE if objects share the same UUID, FALSE otherwise
#'   \item \code{uuid1}: UUID of first object
#'   \item \code{uuid2}: UUID of second object
#'   \item \code{common_steps}: Character vector of steps that are identical in both objects
#'   \item \code{branching_point}: Index of the last common step (or 0 if no common steps)
#'   \item \code{divergent_steps_obj1}: Character vector of steps unique to object 1
#'   \item \code{divergent_steps_obj2}: Character vector of steps unique to object 2
#'   \item \code{summary}: Character string with human-readable summary
#' }
#' @author Geremy Clair
#' @export
checkRelationRomicsObjects <- function(romics_object1, romics_object2, verbose = TRUE) {

  # Validate inputs
  if (!is.romicsObject(romics_object1)) {
    stop("romics_object1 is not a valid romics object")
  }
  if (!is.romicsObject(romics_object2)) {
    stop("romics_object2 is not a valid romics object")
  }

  # Extract UUIDs
  uuid1 <- romics_object1$uuid
  uuid2 <- romics_object2$uuid

  # Check if objects are related (same origin)
  related <- identical(uuid1, uuid2)

  # Initialize result
  result <- list(
    related = related,
    uuid1 = uuid1,
    uuid2 = uuid2,
    common_steps = character(),
    branching_point = 0,
    divergent_steps_obj1 = character(),
    divergent_steps_obj2 = character(),
    summary = ""
  )

  if (!related) {
    result$summary <- paste(
      "Objects are NOT related - they have different origins.",
      "\n  Object 1 UUID:", uuid1,
      "\n  Object 2 UUID:", uuid2
    )

    if (verbose) {
      cat(result$summary, "\n")
    }
    return(result)
  }

  # Objects are related - compare steps
  steps1 <- romics_object1$steps
  steps2 <- romics_object2$steps

  # Find common steps
  min_length <- min(length(steps1), length(steps2))
  common_count <- 0

  for (i in 1:min_length) {
    if (identical(steps1[i], steps2[i])) {
      common_count <- i
    } else {
      break
    }
  }

  result$common_steps <- steps1[1:common_count]
  result$branching_point <- common_count

  # Get divergent steps
  if (common_count < length(steps1)) {
    result$divergent_steps_obj1 <- steps1[(common_count + 1):length(steps1)]
  }

  if (common_count < length(steps2)) {
    result$divergent_steps_obj2 <- steps2[(common_count + 1):length(steps2)]
  }

  # Create summary
  if (common_count == 0) {
    result$summary <- paste(
      "Objects are related (same origin UUID) but diverged immediately.",
      "\n  No common processing steps found.",
      "\n  Object 1 has", length(steps1), "steps",
      "\n  Object 2 has", length(steps2), "steps"
    )
  } else if (common_count == length(steps1) && common_count == length(steps2)) {
    result$summary <- paste(
      "Objects are IDENTICAL - same origin and identical processing history.",
      "\n  Total common steps:", common_count
    )
  } else if (common_count == min(length(steps1), length(steps2))) {
    # One is a subset of the other
    shorter_obj <- if (length(steps1) < length(steps2)) "Object 1" else "Object 2"
    longer_obj <- if (length(steps1) > length(steps2)) "Object 1" else "Object 2"
    additional_steps <- abs(length(steps1) - length(steps2))

    result$summary <- paste(
      "Objects are related with linear evolution.",
      "\n  Common steps:", common_count,
      "\n", shorter_obj, "is a predecessor of", longer_obj,
      "\n", longer_obj, "has", additional_steps, "additional processing step(s) after the branching point"
    )
  } else {
    # True branching - both have divergent steps
    result$summary <- paste(
      "Objects are related but have DIVERGED after step", common_count, ":",
      "\n  Common processing history:", common_count, "steps",
      "\n  Branching point: After step -", steps1[common_count],
      "\n  Object 1 then performed:", length(result$divergent_steps_obj1), "additional step(s)",
      "\n  Object 2 then performed:", length(result$divergent_steps_obj2), "additional step(s)"
    )
  }

  if (verbose) {
    cat("\n════════════════════════════════════════════\n")
    cat("Relationship Analysis: Comparing Two Romics Objects\n")
    cat("════════════════════════════════════════════\n\n")

    cat("Shared UUID (Origin):", uuid1, "\n")
    cat("Objects are related: YES\n\n")

    cat(result$summary, "\n\n")

    if (length(result$common_steps) > 0) {
      cat("Common Processing Steps:\n")
      for (i in 1:length(result$common_steps)) {
        cat(sprintf("  %d. %s\n", i, result$common_steps[i]))
      }
      cat("\n")
    }

    if (length(result$divergent_steps_obj1) > 0) {
      cat("Divergent Steps in Object 1:\n")
      for (i in 1:length(result$divergent_steps_obj1)) {
        cat(sprintf("  %d. %s\n", common_count + i, result$divergent_steps_obj1[i]))
      }
      cat("\n")
    }

    if (length(result$divergent_steps_obj2) > 0) {
      cat("Divergent Steps in Object 2:\n")
      for (i in 1:length(result$divergent_steps_obj2)) {
        cat(sprintf("  %d. %s\n", common_count + i, result$divergent_steps_obj2[i]))
      }
      cat("\n")
    }

    cat("════════════════════════════════════════════\n")
  }

  return(invisible(result))
}
