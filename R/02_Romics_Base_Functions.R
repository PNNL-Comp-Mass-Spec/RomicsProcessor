#' createRomicsObject()
#' @description Create a romics_object by combining a data and a metadata data frames
#' @param data A data frame corresponding to the data.
#' @param metadata A data frame corresponding to the metadata, the columns must be the same as the data ones.
#' @param IDs A data frame corresponding to the list of alternative IDs for the romics_object, this enables to use human readable ids for certain plots. The first column contains the same type of IDs as the imported data.
#' @param main_factor The rowname of a metadata factor, by default the first row will be used as main factor.
#' @param custom_colors A character vector containing the colors you want to use for your figures generated from an romics_object.
#' @param omics_type A character vector of length 1 indicating the type of omics data used.
#' @param omics_information A table or vector containing omics information such as pre-processing details employed.
#' @return This function generate an romics_object containing the following layers : data, metadata, missingdata, original_data, original_metadata, main_factor, colors, steps, custom_colors, omics_type, omics_information, uuid
#' @examples \dontrun{
#' ROP_romics_object <-
#' createRomicsObject(ROP_data, ROP_metadata, main_factor="growth_media")
#' }
#' @author Geremy Clair
#' @export
createRomicsObject <- function(data, metadata, IDs, main_factor="none", custom_colors=ROP_colors, omics_type="unknown", omics_information="unknown"){
  arguments <- as.list(match.call())
  funName <- arguments[[1]]

  if(missing(data)){stop("<data> is missing")}
  if(missing(metadata)) {stop("<metadata> is missing")}
  if(sum(is.na(metadata))>0){stop("<metadata> should not contain NAs.")}
  if(missing(custom_colors)) {custom_colors<-ROP_colors}
  if(missing(omics_type)) {omics_type<-"unknown"}
  if(missing(omics_information)){omics_information<-"unknown"}
  if(missing(main_factor)){main_factor<-"none"}
  if(typeof(data)!="list"){data<-data.frame(data)}
  if(typeof(metadata)!="list"){metadata<-data.frame(metadata)}

  argumentsNames<-names(arguments)[2:length(arguments)]
  colnames(data)[1]<-"ID"
  colnames(metadata)[1]<-"ID"
  if(!sum(colnames(data) %in% colnames(metadata))==ncol(metadata)){stop("data and metadata columns names are different or don't have the same number of columns")}
  metadata<- metadata[, match(colnames(data), colnames(metadata))]
  data_IDs<-data[,1]
  rownames(metadata)<-metadata[,1]
  data<-data[,-1]
  metadata<-metadata[,-1]
  data <- as.data.frame(lapply(data, as.double))
  rownames(data)<-data_IDs
  original_data<-data
  original_metadata<-metadata  # Add original_metadata here
  metadata<-metadata[,order(match(colnames(metadata),colnames(data)))]
  missingdata<- data
  missingdata<-data.frame(is.na(data))
  if(!missing(IDs)){
    if(!is.data.frame(IDs) & sum(rownames(data) %in% IDs[,1])!=nrow(data)){stop("'IDs'has to be a data.frame with all the ids of the data included in the first column")}
    IDs<-merge(data.frame(original=rownames(data)),IDs,by.x=1,by.y=1)
    rownames(IDs)<-IDs[,1]
  }else{IDs<-"unknown"}
  if(main_factor=="none"){
    warning("your main_factor was missing the first row of your metadata was used as factor")
    main_factor<-metadata[1,]}else{
      if(sum(rownames(metadata)==main_factor)==1){ main_factor<-metadata[rownames(metadata)==main_factor,]}else {
        stop("Your main_factor is not present in your metadata")
      }
    }
  factor<-rownames(main_factor)
  main_factor<- as.character(as.factor(t(main_factor)[,1]))
  names(main_factor)<- colnames(metadata)
  main_factor_lvl<- data.frame(unique(main_factor))
  tp_data<-data.frame(matrix(nrow=nrow(data),ncol=0))
  tp_metadata<-data.frame(matrix(nrow=nrow(metadata),ncol=0))
  tp_missingdata<-data.frame(matrix(nrow=nrow(missingdata),ncol=0))
  tp_main_factor<-character()
  for(i in 1:nrow(main_factor_lvl)){
    tp_data<-cbind(tp_data,data[,main_factor==as.character(main_factor_lvl[i,])])
    tp_metadata<- cbind(tp_metadata,metadata[,main_factor==as.character(main_factor_lvl[i,])])
    tp_missingdata<-cbind(tp_missingdata,missingdata[,main_factor==as.character(main_factor_lvl[i,])])
    tp_main_factor<- c(tp_main_factor,main_factor[main_factor==as.character(main_factor_lvl[i,])])
  }
  colnames(tp_data)<-colnames(tp_metadata)<-colnames(tp_missingdata)<-names(tp_main_factor)
  rownames(tp_data)<-rownames(data)
  rownames(tp_metadata)<-rownames(metadata)
  data<-tp_data
  metadata<-tp_metadata
  missingdata<-tp_missingdata
  main_factor<-tp_main_factor
  remove(tp_data,tp_metadata,tp_missingdata)
  if(length(custom_colors)<nrow(main_factor_lvl)){
    warning("your color vector is shorter than the number of factors selected, some colors will be picked automatically")
    custom_colors<- c(custom_colors,ROP_colors)
  }
  main_factor_lvl$colors_romics<- custom_colors[1:nrow(main_factor_lvl)]
  colors_romics<-data.frame(main_factor_lvl[match(main_factor,main_factor_lvl[,1]),2])
  colors_romics<-t(colors_romics)
  rownames(colors_romics)<-"colors_romics"
  colnames(colors_romics)<- colnames(metadata)
  metadata<-rbind(metadata,colors_romics)
  colors_romics<- t(colors_romics)
  colors <- rep(as.character(colors_romics), each = nrow(data))
  steps<- c("romics_object",paste0("date|",gsub(" ","_",format(Sys.time(),"%b_%d_%Y_%X")),"|createRomicsObject"))
  fun<-paste0(funName,"(")
  for(i in 1:length(argumentsNames)){
    fun<-paste0(fun,argumentsNames[i],"=")
    if(typeof(arguments[[i+1]])=="character"){
      fun<-paste0(fun,"'",unlist(arguments[[i+1]]),"'")
    }else{
      fun<-paste0(fun,unlist(arguments[[i+1]]))
    }
    if(i<length(argumentsNames)){fun<-paste0(fun,",")}
  }
  fun<-paste0("fun|",fun,")")
  steps<- c(steps,fun)
  dependencies<-romicsCreateDependencies()

  # Create the list WITHOUT uuid initially (since romicsAttributeUUID will add it)
  l<-list(data, metadata, missingdata, original_data, original_metadata, IDs, factor, colors, steps, dependencies, custom_colors, omics_type, omics_information)
  names(l)<-c("data","metadata","missingdata","original_data", "original_metadata", "IDs", "main_factor","colors","steps","dependencies","custom_colors","omics_type","omics_information")
  class(l)<-"romics_object"
  romics_object<-l

  # Generate UUID using your function
  romics_object <- romicsAttributeUUID(romics_object)
  return(romics_object)
}

#' resetRomicsObject()
#' @description Reset a romics_object to its original state using createRomicsObject and preserving steps/UUID
#' @param romics_object A romics_object created using createRomicsObject()
#' @param main_factor Character string indicating the metadata row to use as main factor. If not provided, the first row will be used.
#' @details Restores the romics_object to its original state by recreating it with createRomicsObject
#' and then restoring the original steps (up to createRomicsObject) and UUID
#' @return A reset romics_object in its original state
#' @author Geremy Clair
#' @export
resetRomicsObject <- function(romics_object, main_factor = NULL) {
  # Input validation
  if(!is.romicsObject(romics_object) || missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }

  if(is.null(romics_object$original_data)) {
    stop("romics_object does not contain an 'original_data' layer")
  }

  if(is.null(romics_object$original_metadata)) {
    stop("romics_object does not contain an 'original_metadata' layer")
  }

  message("Resetting romics_object to original state...")

  # Preserve what needs to be restored
  original_uuid <- romics_object$uuid
  original_steps <- romics_object$steps
  original_IDs <- romics_object$IDs
  original_omics_information <- romics_object$omics_information
  original_custom_colors <- romics_object$custom_colors
  original_omics_type <- romics_object$omics_type

  # Get original data and metadata
  original_data <- romics_object$original_data
  original_metadata <- romics_object$original_metadata

  # Prepare data for createRomicsObject (add ID column)
  data_with_id <- cbind(ID = rownames(original_data), original_data)
  metadata_with_id <- cbind(ID = rownames(original_metadata), original_metadata)

  # Determine main_factor
  if(is.null(main_factor)) {
    main_factor_to_use <- rownames(original_metadata)[1]
    message(paste("No main_factor specified, using first row:", main_factor_to_use))
  } else {
    if(!main_factor %in% rownames(original_metadata)) {
      stop(paste("main_factor '", main_factor, "' not found in metadata. Available factors:",
                 paste(rownames(original_metadata), collapse = ", ")))
    }
    main_factor_to_use <- main_factor
    message(paste("Using main_factor:", main_factor_to_use))
  }

  # Filter IDs table if it exists and clean tracking columns
  IDs_to_use <- "unknown"
  if(!is.null(original_IDs)) {
    # Remove tracking columns
    cols_to_remove <- c("Original_IDs", "Current_IDs")
    cols_to_keep <- setdiff(colnames(original_IDs), cols_to_remove)
    IDs_filtered <- original_IDs[, cols_to_keep, drop = FALSE]

    # Filter to only include IDs present in original_data
    original_ids <- rownames(original_data)
    ids_in_original <- intersect(rownames(IDs_filtered), original_ids)

    if(length(ids_in_original) > 0) {
      IDs_to_use <- IDs_filtered[ids_in_original, , drop = FALSE]
      # Add ID column for createRomicsObject
      IDs_to_use <- cbind(ID = rownames(IDs_to_use), IDs_to_use)
      message(paste("Filtered IDs table to", nrow(IDs_to_use), "entries"))
    } else {
      warning("No IDs in IDs table match original_data, using default")
      IDs_to_use <- "unknown"
    }
  }

  # Recreate the object using createRomicsObject
  if(is.data.frame(IDs_to_use)) {
    reset_object <- createRomicsObject(
      data = data_with_id,
      metadata = metadata_with_id,
      IDs = IDs_to_use,
      main_factor = main_factor_to_use,
      custom_colors = original_custom_colors,
      omics_type = original_omics_type,
      omics_information = original_omics_information
    )
  } else {
    reset_object <- createRomicsObject(
      data = data_with_id,
      metadata = metadata_with_id,
      main_factor = main_factor_to_use,
      custom_colors = original_custom_colors,
      omics_type = original_omics_type,
      omics_information = original_omics_information
    )
  }

  # Restore original UUID (overwrite the new one)
  reset_object$uuid <- original_uuid

  # Clean up steps - keep only up to createRomicsObject from original
  createRomicsObject_index <- which(grepl("createRomicsObject", original_steps))[1]

  if(!is.na(createRomicsObject_index)) {
    original_steps_cleaned <- original_steps[1:(createRomicsObject_index + 1)]
  } else {
    warning("Could not find createRomicsObject in original steps")
    original_steps_cleaned <- original_steps[1:2]  # Keep first two steps as fallback
  }

  # Replace the new steps with cleaned original steps
  reset_object$steps <- original_steps_cleaned

  # Add the reset operation to steps
  reset_step <- paste0("date|", gsub(" ", "_", format(Sys.time(), "%b_%d_%Y_%X")), "|resetRomicsObject")
  reset_object$steps <- c(reset_object$steps, reset_step)

  message("Successfully reset romics_object to original state")
  message(paste("UUID preserved:", reset_object$uuid))
  message(paste("Steps cleaned to", length(reset_object$steps), "entries"))

  return(reset_object)
}

#' romicsAttributeUUID()
#' @description Add or update a UUID to a romics_object
#' @param romics_object A romics_object to which a UUID will be added or updated
#' @param force_new Logical, if TRUE will generate a new UUID even if one already exists (default: FALSE)
#' @return The romics_object with a UUID added to the "uuid" layer
#' @examples \dontrun{
#' # Add UUID to existing romics_object
#' romics_object <- romicsAttributeUUID(romics_object)
#'
#' # Force generation of new UUID
#' romics_object <- romicsAttributeUUID(romics_object, force_new = TRUE)
#' }
#' @author Geremy Clair
#' @export
romicsAttributeUUID <- function(romics_object, force_new = FALSE) {
  # Check if UUID already exists and force_new is FALSE
  if ("uuid" %in% names(romics_object) && !force_new) {
    message("UUID already exists. Use force_new = TRUE to generate a new UUID.")
    return(romics_object)
  }

  # Generate new UUID
  new_uuid <- uuid::UUIDgenerate()
  romics_object$uuid <- new_uuid

  message(paste0("The following Universal Unique Identifier (UUID) was attributed to the romics_object: ",romics_object$uuid))
  # Add UUID to the object

  return(romics_object)
}

#' is.romicsObject()
#' @description Enables to check if the romics_object is in the appropriate format. Stops execution if UUID is missing.
#' @param romics_object A romics_object created using createRomicsObject()
#' @return This function will return TRUE if the object is a valid romics_object, or stop execution with an error message
#' @author Geremy Clair
#' @export
is.romicsObject <- function(romics_object){
  if(!inherits(romics_object, "romics_object")){
    stop("Your romics_object was not created using the function createRomicsObject()")
  }
  if(romics_object$steps[1] != "romics_object"){
    stop("Your romics_object is not in the appropriate format")
  }
  if(!"uuid" %in% names(romics_object)){
    stop("Your romics_object does not contain a UUID, which is required for FAIR compliance and reproducibility.\n",
         "Please add a UUID using: romics_object <- romicsAttributeUUID(romics_object)\n",
         "Or recreate the object using the updated createRomicsObject() function.")
  }

  return(TRUE)
}

#' romicsChangeFactor()
#' @description Change the main factor of the romics_object
#' @param romics_object A object created using the function createRomicsObject().
#' @param main_factor Either 'none' OR a factor from the romics_object (corresponding to a row from the original metadata file), the list of factors from an romics object can be obtained using the function romicsFactorNames().
#' @details changes the main_factor of an romics_object and updates the colors to this new factor. Also reorders the embeddings layer if present.
#' @return This function returns a modified romics object, please see the create_romics_object() documentation.
#' @author Geremy Clair
#' @export
romicsChangeFactor <- function(romics_object, main_factor = "none") {
  arguments <- as.list(match.call())
  if(!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }
  if(missing(main_factor) | main_factor == "none"){
    warning("The first row of your metadata was used as factor")
    warning(paste0("main_factor='", romicsFactorNames(romics_object)[1], "'"))
    main_factor <- romicsFactorNames(romics_object)[1]
  }
  if(!main_factor %in% romicsFactorNames(romics_object)){
    warning("Your main_factor has to be present in your metadata, the following list of factor are in the romics_object:")
    print(romicsFactorNames(romics_object))
    return(romics_object)
  }

  data <- romics_object$data
  metadata <- romics_object$metadata
  missingdata <- romics_object$missingdata
  embeddings <- romics_object$embeddings  # Extract embeddings if they exist
  f_main_factor <- main_factor

  # Update romics_object$main_factor
  romics_object$main_factor <- main_factor

  # Extract the main factor
  main_factor_values <- romicsExtractFactor(romics_object, factor = main_factor)

  # Establish the levels of the main factor
  main_factor_lvl <- levels(main_factor_values)
  n <- length(main_factor_lvl)

  f_data <- list()
  f_metadata <- list()
  f_missingdata <- list()
  f_embeddings <- list()  # Add embeddings list

  for(i in seq_along(main_factor_lvl)){
    # Obtain logical index
    index <- main_factor_values == as.character(main_factor_lvl[i])
    # Extract columns directly
    tp_data <- data[, index, drop = FALSE]
    tp_metadata <- metadata[, index, drop = FALSE]
    tp_missingdata <- missingdata[, index, drop = FALSE]

    # ADDED: Reorder embeddings if they exist
    if(!is.null(embeddings)) {
      tp_embeddings <- embeddings[, index, drop = FALSE]
      f_embeddings[[i]] <- tp_embeddings
    }

    # Append to preallocated list
    f_data[[i]] <- tp_data
    f_metadata[[i]] <- tp_metadata
    f_missingdata[[i]] <- tp_missingdata
  }

  # Combine the list elements into final data frames
  f_data <- do.call(cbind, f_data)
  f_metadata <- do.call(cbind, f_metadata)
  f_missingdata <- do.call(cbind, f_missingdata)

  # ADDED: Combine embeddings if they exist
  if(!is.null(embeddings)) {
    f_embeddings <- do.call(cbind, f_embeddings)
  }

  data <- f_data
  metadata <- f_metadata
  missingdata <- f_missingdata
  if(!is.null(embeddings)) {
    embeddings <- f_embeddings
  }

  # Remove the previous color line in the metadata
  metadata <- metadata[rownames(metadata) != "colors_romics", ]

  # IMPROVED COLOR MANAGEMENT - FIXED VERSION
  custom_colors <- romics_object$custom_colors

  # Ensure we have enough colors for all factor levels
  if(length(custom_colors) < length(main_factor_lvl)){
    warning("Your color vector is shorter than the number of factors selected, some colors will be picked automatically")
    # Extend custom_colors with ROP_colors if needed
    additional_colors_needed <- length(main_factor_lvl) - length(custom_colors)
    if(exists("ROP_colors")) {
      custom_colors <- c(custom_colors, ROP_colors[1:additional_colors_needed])
    } else {
      # Fallback if ROP_colors doesn't exist
      additional_colors <- rainbow(additional_colors_needed)
      custom_colors <- c(custom_colors, additional_colors)
    }
  }

  # FIXED: Assign one unique color per factor level
  level_colors <- custom_colors[1:length(main_factor_lvl)]
  names(level_colors) <- main_factor_lvl

  # FIXED: Get factor values in the NEW ORDER (after reordering by factor)
  reordered_factor_values <- as.character(t(metadata[rownames(metadata) == f_main_factor, ]))

  # FIXED: Create sample colors by mapping each sample's factor level to its assigned color
  sample_colors <- character(length(reordered_factor_values))
  for(i in seq_along(reordered_factor_values)) {
    factor_level <- reordered_factor_values[i]
    sample_colors[i] <- level_colors[factor_level]
  }
  names(sample_colors) <- colnames(metadata)

  # VALIDATION: Check for color consistency
  unique_factor_levels <- unique(reordered_factor_values)
  for(factor_level in unique_factor_levels) {
    level_samples <- reordered_factor_values == factor_level
    level_sample_colors <- sample_colors[level_samples]
    unique_colors_for_level <- unique(level_sample_colors)

    if(length(unique_colors_for_level) > 1) {
      stop(paste("Error: Factor level", factor_level, "has multiple colors assigned:",
                 paste(unique_colors_for_level, collapse = ", ")))
    }
  }

  # Create colors_romics row for metadata
  colors_romics <- data.frame(matrix(ncol = ncol(metadata), nrow = 1))
  colors_romics[1, ] <- sample_colors
  colnames(colors_romics) <- colnames(metadata)

  # Add colors_romics to metadata
  metadata <- rbind(metadata, colors_romics)
  rownames(metadata)[nrow(metadata)] <- "colors_romics"

  # Create a vector containing the colors for the whole data points and store it in colors layer
  colors <- rep(sample_colors, each = nrow(data))

  # Replace the different transformed layers in the romics_object
  romics_object$data <- data
  romics_object$metadata <- metadata
  romics_object$missingdata <- missingdata
  romics_object$colors <- colors
  romics_object$main_factor <- f_main_factor

  # ADDED: Update embeddings layer if it exists
  if(!is.null(embeddings)) {
    romics_object$embeddings <- embeddings
    message("Embeddings layer reordered to match new factor arrangement")
  }

  # Update the steps
  romics_object <- romicsUpdateSteps(romics_object, arguments)

  # Print summary of color assignment with validation
  message("Color assignment summary:")
  for(i in seq_along(main_factor_lvl)) {
    level_name <- main_factor_lvl[i]
    level_color <- level_colors[level_name]
    level_count <- sum(reordered_factor_values == level_name)
    message(paste("  ", level_name, ":", level_color, "(", level_count, "samples)"))
  }

  # FINAL VALIDATION: Verify no duplicate colors
  if(length(unique(level_colors)) != length(level_colors)) {
    warning("Warning: Some factor levels have been assigned the same color!")
    dup_colors <- level_colors[duplicated(level_colors) | duplicated(level_colors, fromLast = TRUE)]
    message("Duplicate colors found: ", paste(unique(dup_colors), collapse = ", "))
  } else {
    message("✓ All factor levels have unique colors")
  }

  return(romics_object)
}

#' romicsChangeIDs()
#' @description Change the main IDs of the romics_object using a specified column from the IDs layer
#' @param romics_object A romics_object created using createRomicsObject()
#' @param newIDs Character string. Must match exactly one of the column names of the 'IDs' layer
#' @details Changes the row names (IDs) of the data, missingdata, and statistics layers using the specified column
#' from the IDs layer. Always stores the conversion history with "original" and "converted" columns.
#' Updates the IDs table with "Original_IDs" and "Current_IDs" columns to enable reverse mapping.
#' Missing mappings are handled with "Undefined_ID_matched_to_" prefix.
#' @return This function returns a modified romics object with updated IDs and conversion history stored.
#' @author Geremy Clair
#' @export
romicsChangeIDs <- function(romics_object, newIDs) {
  arguments <- as.list(match.call())

  # Input validation
  if(!is.romicsObject(romics_object) || missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }
  if(is.null(romics_object$IDs)) {
    stop("romics_object does not contain an 'IDs' layer")
  }
  if(missing(newIDs) || length(newIDs) != 1 || !is.character(newIDs)) {
    stop("newIDs must be a character vector of length 1")
  }
  if(!newIDs %in% colnames(romics_object$IDs)) {
    available_cols <- paste(colnames(romics_object$IDs), collapse = ", ")
    stop(paste("newIDs must be a column name from romics_object$IDs. Available columns:", available_cols))
  }

  # Get current IDs (what we're converting FROM)
  current_ids <- rownames(romics_object$data)

  # Get the IDs table
  IDs_table <- romics_object$IDs

  # Initialize Original_IDs column if this is the first conversion
  if(!"Original_IDs" %in% colnames(IDs_table)) {
    romics_object$IDs$Original_IDs <- rownames(IDs_table)
    IDs_table <- romics_object$IDs  # Update local copy
    message("Created Original_IDs column with initial row names")
  }

  # Find mapping between current_ids and the IDs table
  mapping_indices <- NULL
  mapping_method <- NULL

  # Method 1: Try to match current_ids with "Current_IDs" column (if it exists)
  if("Current_IDs" %in% colnames(IDs_table)) {
    if(all(current_ids %in% IDs_table$Current_IDs)) {
      mapping_indices <- match(current_ids, IDs_table$Current_IDs)
      mapping_method <- "Current_IDs column"
      message("Mapping using Current_IDs column")
    }
  }

  # Method 2: Try to match current_ids with row names of IDs table
  if(is.null(mapping_indices)) {
    if(all(current_ids %in% rownames(IDs_table))) {
      mapping_indices <- match(current_ids, rownames(IDs_table))
      mapping_method <- "row names"
      message("Mapping using IDs table row names")
    }
  }

  # Method 3: Try to match with Original_IDs column
  if(is.null(mapping_indices)) {
    if(all(current_ids %in% IDs_table$Original_IDs)) {
      mapping_indices <- match(current_ids, IDs_table$Original_IDs)
      mapping_method <- "Original_IDs column"
      message("Mapping using Original_IDs column")
    }
  }

  # Method 4: Try to match with each other column in IDs table
  if(is.null(mapping_indices)) {
    for(col_name in colnames(IDs_table)) {
      if(col_name %in% c("Original_IDs", "Current_IDs")) next  # Skip these, already tried
      if(all(current_ids %in% IDs_table[[col_name]])) {
        mapping_indices <- match(current_ids, IDs_table[[col_name]])
        mapping_method <- paste("column", col_name)
        message(paste("Mapping using column:", col_name))
        break
      }
    }
  }

  if(is.null(mapping_indices)) {
    # Debug information
    cat("Debug: Current IDs (first 5):", head(current_ids), "\n")
    cat("Available columns in IDs table:", colnames(IDs_table), "\n")
    if("Current_IDs" %in% colnames(IDs_table)) {
      cat("Current_IDs column (first 5):", head(IDs_table$Current_IDs), "\n")
    }
    if("Original_IDs" %in% colnames(IDs_table)) {
      cat("Original_IDs column (first 5):", head(IDs_table$Original_IDs), "\n")
    }
    stop("Cannot find mapping between current IDs and IDs table")
  }

  # Extract new IDs from the specified column
  converted_ids <- IDs_table[[newIDs]][mapping_indices]

  # Handle missing/empty values with the specific prefix you want
  missing_mask <- is.na(converted_ids) | converted_ids == "" | converted_ids == " "
  if(any(missing_mask)) {
    converted_ids[missing_mask] <- paste0("Undefined_ID_matched_to_", current_ids[missing_mask])
    message(paste("Created", sum(missing_mask), "undefined ID entries"))
  }

  # Handle duplicates with letter suffixes
  if(any(duplicated(converted_ids))) {
    dup_ids <- unique(converted_ids[duplicated(converted_ids)])
    total_duplicates_fixed <- 0

    for(dup_id in dup_ids) {
      dup_positions <- which(converted_ids == dup_id)
      n_dups <- length(dup_positions)

      if(n_dups > 1) {
        # Add letter suffixes to duplicates (keep first as-is)
        letters_needed <- LETTERS[1:(n_dups-1)]
        if(n_dups > 27) {
          # Extend with double letters if needed
          extra_letters <- paste0(rep(LETTERS, each = 26), LETTERS)[1:(n_dups-27)]
          letters_needed <- c(letters_needed, extra_letters)
        }

        for(i in 2:n_dups) {
          converted_ids[dup_positions[i]] <- paste0(converted_ids[dup_positions[i]], "_", letters_needed[i-1])
        }
        total_duplicates_fixed <- total_duplicates_fixed + (n_dups - 1)
      }
    }
    message(paste("Resolved", total_duplicates_fixed, "duplicate IDs using letter suffixes"))
  }

  # Create conversion history with "original" and "converted" columns
  conversion_history <- data.frame(
    original = current_ids,
    converted = converted_ids,
    conversion_date = Sys.Date(),
    column_used = newIDs,
    mapping_method = mapping_method,
    stringsAsFactors = FALSE
  )

  # Store conversion history in romics_object
  if(is.null(romics_object$ID_conversions)) {
    romics_object$ID_conversions <- list()
  }

  conversion_name <- paste0("conversion_", length(romics_object$ID_conversions) + 1, "_to_", newIDs)
  romics_object$ID_conversions[[conversion_name]] <- conversion_history

  # Update row names in all data layers
  rownames(romics_object$data) <- converted_ids
  rownames(romics_object$missingdata) <- converted_ids

  # Update statistics layer if it exists
  if("statistics" %in% names(romics_object) && !is.null(romics_object$statistics)) {
    if(nrow(romics_object$statistics) == length(converted_ids)) {
      rownames(romics_object$statistics) <- converted_ids
    }
  }

  # Update embeddings layer if it exists
  if("embeddings" %in% names(romics_object) && !is.null(romics_object$embeddings)) {
    if(ncol(romics_object$embeddings) == length(converted_ids)) {
      colnames(romics_object$embeddings) <- converted_ids
    }
  }

  # Update Current_IDs column in IDs table
  # Initialize Current_IDs if it doesn't exist
  if(!"Current_IDs" %in% colnames(romics_object$IDs)) {
    romics_object$IDs$Current_IDs <- romics_object$IDs$Original_IDs  # Initialize with original
    message("Created Current_IDs column initialized with Original_IDs")
  }

  # Update Current_IDs for the mapped rows
  romics_object$IDs$Current_IDs[mapping_indices] <- converted_ids

  # Success messages
  message(paste("Successfully converted", length(converted_ids), "IDs"))
  message(paste("Conversion history stored as:", conversion_name))
  message("Updated Current_IDs column in IDs table for reverse mapping")

  # Update steps
  romics_object <- romicsUpdateSteps(romics_object, arguments)

  return(romics_object)
}

#' romicsChangeColors()
#' @description Change the colors of an existing romics_object based on the main factor levels (performance optimized)
#' @param romics_object A romics_object created using createRomicsObject()
#' @param new_colors A character vector containing the new colors you want to use. The colors will be applied in the order of the main factor levels.
#' @details This function updates the colors in a romics_object by replacing the existing color scheme with new colors. Uses vectorized operations for optimal performance with large datasets.
#' @return This function returns the romics_object with updated colors in the metadata, colors, and custom_colors layers
#' @author Geremy Clair
#' @export
romicsChangeColors <- function(romics_object, new_colors) {

  # Capture arguments using your standard procedure
  arguments <- as.list(match.call())

  # Input validation
  if (!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }

  if (missing(new_colors)) {
    stop("new_colors is missing. Please provide a vector of colors.")
  }

  if (!is.character(new_colors)) {
    stop("new_colors must be a character vector")
  }

  # Extract the main factor using romicsExtractFactor function
  main_factor <- romicsExtractFactor(romics_object, factor = romics_object$main_factor)
  main_factor_levels <- levels(main_factor)
  n_levels <- length(main_factor_levels)

  # Check if enough colors are provided
  if (length(new_colors) < n_levels) {
    warning(paste("Only", length(new_colors), "colors provided for", n_levels,
                  "factor levels. Colors will be recycled."))
    new_colors <- rep(new_colors, length.out = n_levels)
  }

  #Create named vector for fast lookup
  level_colors <- setNames(new_colors[1:n_levels], main_factor_levels)

  #Vectorized color mapping - single operation
  colors_romics_vector <- level_colors[as.character(main_factor)]


  # Direct assignment to metadata row
  if ("colors_romics" %in% rownames(romics_object$metadata)) {
    romics_object$metadata["colors_romics", ] <- colors_romics_vector
  } else {
    # Pre-allocate and assign in one operation
    romics_object$metadata <- rbind(romics_object$metadata,
                                    colors_romics = colors_romics_vector)
  }

  romics_object$colors<-new_colors
  names(romics_object$colors)<-main_factor_levels

  # Update custom_colors
  romics_object$custom_colors <- new_colors

  # Update steps using your standard procedure
  romics_object <- romicsUpdateSteps(romics_object, arguments)

  return(romics_object)
}

#' romicsSubset()
#' @description Keeps or drops a subset of specific elements/columns from the romics_object
#' @param romics_object A romics_object created using createRomicsObject()
#' @param subset_vector A character vector of factor levels or colnames to keep in the object.
#' @param filter_mode Either 'keep' or 'drop' to indicate if you want to conserve or to drop the elements from a given factor
#' @param by Either 'colnames' or 'level' to indicate what you want to keep or drop.
#' @param factor A factor contained in the metadata of the romics_object, to obtain the list of factors please use the function romicsFactorNames()
#' @param handle_embeddings Either 'remove' (default) or 'subset'. How to handle embeddings when subsetting.
#' @details This function creates a new object based on a previous romics_object to include or drop a list of specified columns from the original object. The created object will have a new step object created that will indicate the name of the original object to be subsetted and the log/non-log status of the object.
#' @details Note that this function will remove the stat layer from your object
#' @return This function generates a subsetted romics_object
#' @author Geremy Clair
#' @export
romicsSubset <- function(romics_object, subset_vector, filter_mode = "keep", by = "colnames",
                         factor = "main", handle_embeddings = "remove") {
  arguments <- as.list(match.call())

  if(!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }
  if(missing(subset_vector)) {stop("Your subsetting vector is missing")}
  if(!is.character(subset_vector)) {stop("Your subset vector has to be a character vector")}
  if(missing(filter_mode)) {filter_mode <- "keep"}
  if(!filter_mode %in% c("keep", "drop")) {stop("<filter_mode> has to be either 'keep' or 'drop'")}
  if(missing(by)) {by <- "colnames"}
  if(!by %in% c("colnames", "level")) {stop("<by> has to be either 'colnames' or 'level'")}
  if(!handle_embeddings %in% c("remove", "subset")) {
    stop("<handle_embeddings> has to be either 'remove' or 'subset'")
  }
  if(missing(factor)) {factor <- "main"}

  if(factor == "main") {
    factor <- romics_object$main_factor
  } else {
    if(!factor %in% rownames(romics_object$metadata)) {
      stop("Your factor is not present in the metadata")
    }
  }

  # Store original sample count
  original_n_samples <- ncol(romics_object$data)

  # Extract factor
  fac <- romics_object$metadata[factor == rownames(romics_object$metadata), ]
  fac <- as.factor(t(fac))

  # Create a logical vector containing the columns to keep (TRUE)
  if(by == "colnames") {
    if(sum(subset_vector %in% colnames(romics_object$data)) != length(subset_vector)) {
      warning("Not all the elements of the subset_vector were present in the colnames of the data")
    }
    if(filter_mode == "keep") {
      sub_logical <- colnames(romics_object$data) %in% subset_vector
    } else {
      sub_logical <- !colnames(romics_object$data) %in% subset_vector
    }
  } else {
    if(sum(subset_vector %in% fac) != length(subset_vector)) {
      warning("Not all the elements of the subset_vector were levels of the factor")
    }
    if(filter_mode == "keep") {
      sub_logical <- fac %in% subset_vector
    } else {
      sub_logical <- !fac %in% subset_vector
    }
  }

  # Get column names to keep
  cols_to_keep <- colnames(romics_object$data)[sub_logical]
  n_samples_kept <- length(cols_to_keep)
  n_samples_removed <- original_n_samples - n_samples_kept
  pct_removed <- round((n_samples_removed / original_n_samples) * 100, 1)

  # Print summary message first
  message(paste0("Subsetting complete: ", n_samples_kept, " samples retained (",
                 pct_removed, "% of samples removed)."))

  # Remove the columns from the element data, metadata, and missingness of the romics_object
  romics_object$data <- romics_object$data[, sub_logical]
  romics_object$metadata <- romics_object$metadata[, sub_logical]
  romics_object$missingdata <- romics_object$missingdata[, sub_logical]

  # Handle embeddings
  if(!is.null(romics_object$embeddings) && length(romics_object$embeddings) > 0) {
    if(handle_embeddings == "remove") {
      romics_object$embeddings <- NULL
      message("Embeddings removed from the subsetted object.")
    } else if(handle_embeddings == "subset") {
      for(emb_name in names(romics_object$embeddings)) {
        emb_data <- romics_object$embeddings[[emb_name]]
        # Check if rownames match column names
        if(!is.null(rownames(emb_data)) && all(rownames(emb_data) %in% colnames(romics_object$metadata))) {
          romics_object$embeddings[[emb_name]] <- emb_data[rownames(emb_data) %in% cols_to_keep, , drop = FALSE]
        } else {
          warning(paste("Could not subset embedding:", emb_name, ". Removing it instead."))
          romics_object$embeddings[[emb_name]] <- NULL
        }
      }
      # Remove embeddings list if empty
      if(length(romics_object$embeddings) == 0) {
        romics_object$embeddings <- NULL
      }
      message("Embeddings subsetted to match the data.")
    }
  }

  # Remove the stat layer
  if("statistics" %in% names(romics_object)) {
    warning("The statistics layer was removed from the romics_object; the statistics were calculated on the non-subsetted object")
    romics_object <- romics_object[names(romics_object) != "statistics"]
  }

  class(romics_object) <- "romics_object"

  # Update the colors
  romics_object <- romicsUpdateColor(romics_object)

  # Update the steps
  romics_object <- romicsUpdateSteps(romics_object, arguments)

  return(romics_object)
}


#' romicsChangeSampleNames()
#' @description Changes the sample identifiers using values contained in one of the factors from the romics_object metadata.
#' The new sample names will be taken from the specified factor, which must contain unique values for each sample.
#' @param romics_object A romics_object created using createRomicsObject()
#' @param factor_name A character string specifying the name of a factor contained in the metadata of the romics_object.
#'        This factor must contain unique values that will become the new sample names.
#'        To obtain the list of available factors, use romicsFactorNames().
#' @param original_names Either 'keep' or 'drop'. Indicates what to do with the original sample names:
#'        'keep' will store them as a new factor in the metadata, 'drop' will discard them completely. Default: 'keep'
#' @param original_names_factor The name to give to the new factor that will store the original sample names
#'        (only used when original_names = 'keep'). Default: 'original_sample_names'
#' @details This function enables quick renaming of sample names using values from a factor contained in the metadata.
#'          The specified factor must contain only unique values to ensure each sample gets a distinct name.
#'          Sample names will be updated across all layers of the romics_object (data, metadata, missingdata, embeddings).
#'
#'          When original_names = 'keep', the current sample names are preserved as a new factor in the metadata
#'          before being replaced. When original_names = 'drop', the original names are permanently lost.
#'
#'          The function works with all romics_object layers including embeddings from PCA, UMAP, or t-SNE analyses.
#' @return A romics_object with updated sample names across all layers. If original_names = 'keep',
#'         a new factor containing the original sample names will be added to the metadata.
#' @author Geremy Clair
#' @export
romicsChangeSampleNames <- function(romics_object,
                                    factor_name,
                                    original_names = "keep",
                                    original_names_factor = "original_sample_names") {

  arguments <- as.list(match.call())

  # Input validation
  if (!is.romicsObject(romics_object) || missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }

  if (missing(factor_name)) {
    stop("factor_name is required. Use romicsFactorNames() to see available factors.")
  }

  if (!is.character(factor_name) || length(factor_name) != 1) {
    stop("factor_name must be a character string of length 1")
  }

  if (!factor_name %in% romicsFactorNames(romics_object)) {
    available_factors <- romicsFactorNames(romics_object)
    stop(paste("Factor '", factor_name, "' not found in romics_object metadata.\n",
               "Available factors: ", paste(available_factors, collapse = ", "), sep = ""))
  }

  if (!original_names %in% c("keep", "drop")) {
    warning("original_names must be either 'keep' or 'drop'. Using 'keep' as default.")
    original_names <- "keep"
  }

  if (!is.character(original_names_factor) || length(original_names_factor) != 1) {
    stop("original_names_factor must be a character string of length 1")
  }

  # Extract new names from the specified factor
  new_names <- as.character(romicsExtractFactor(romics_object, factor_name))

  # Check for duplicated values in the factor
  if (any(duplicated(new_names))) {
    duplicated_values <- new_names[duplicated(new_names)]
    cat("The selected factor contains duplicated values:\n")
    print(table(new_names))
    stop(paste("Factor '", factor_name, "' contains duplicated values: ",
               paste(unique(duplicated_values), collapse = ", "),
               "\nSample names must be unique.", sep = ""))
  }

  # Check for missing values
  if (any(is.na(new_names))) {
    stop(paste("Factor '", factor_name, "' contains missing (NA) values. All values must be present for sample renaming."))
  }

  # Verify length matches number of samples
  n_samples <- ncol(romics_object$data)
  if (length(new_names) != n_samples) {
    stop(paste("Length of factor values (", length(new_names),
               ") does not match number of samples (", n_samples, ")", sep = ""))
  }

  # Check which layers exist
  has_missingdata <- !is.null(romics_object$missingdata)
  has_embeddings <- !is.null(romics_object$embeddings)

  # Store original sample names if requested
  if (original_names == "keep") {
    # Check if the original_names_factor already exists
    if (original_names_factor %in% rownames(romics_object$metadata)) {
      warning(paste("Factor '", original_names_factor, "' already exists and will be replaced.", sep = ""))
      romics_object$metadata <- romics_object$metadata[rownames(romics_object$metadata) != original_names_factor, , drop = FALSE]
    }

    # Create new factor with original sample names
    original_colnames_df <- data.frame(matrix(colnames(romics_object$metadata), nrow = 1))
    colnames(original_colnames_df) <- colnames(romics_object$metadata)
    rownames(original_colnames_df) <- original_names_factor

    # Add to metadata (before renaming columns)
    romics_object$metadata <- rbind(original_colnames_df, romics_object$metadata)
  }

  # Apply new sample names to all layers
  message("Updating sample names across romics_object layers...")

  # Core layers (always present)
  colnames(romics_object$data) <- new_names
  colnames(romics_object$metadata) <- new_names

  # Optional layers (if they exist)
  if (has_missingdata) {
    colnames(romics_object$missingdata) <- new_names
    message("✓ Updated missingdata layer")
  }

  if (has_embeddings) {
    colnames(romics_object$embeddings) <- new_names
    message("✓ Updated embeddings layer (PCA/UMAP/t-SNE results)")
  }

  # Update colors and steps
  romics_object <- romicsUpdateColor(romics_object)
  romics_object <- romicsUpdateSteps(romics_object, arguments)

  # Report results
  message(paste("✓ Sample names successfully changed using factor '", factor_name, "'", sep = ""))
  if (original_names == "keep") {
    message(paste("✓ Original sample names stored in factor '", original_names_factor, "'", sep = ""))
  }
  message(paste("✓ Updated", n_samples, "sample names across all data layers"))

  # Final verification of name consistency
  layers_updated <- c("data", "metadata")
  if (has_missingdata) layers_updated <- c(layers_updated, "missingdata")
  if (has_embeddings) layers_updated <- c(layers_updated, "embeddings")

  message(paste("✓ Layers updated:", paste(layers_updated, collapse = ", ")))

  # Verify consistency
  data_names <- colnames(romics_object$data)
  metadata_names <- colnames(romics_object$metadata)

  if (!identical(data_names, metadata_names)) {
    warning("Sample name inconsistency detected between data and metadata layers!")
  }

  if (has_embeddings) {
    embedding_names <- colnames(romics_object$embeddings)
    if (!identical(data_names, embedding_names)) {
      warning("Sample name inconsistency detected between data and embeddings layers!")
    } else {
      message("✓ Sample names are consistent across all layers including embeddings")
    }
  }

  return(romics_object)
}

#' romicsNumericHeaderAdjust()
#' @description Adjusts sample names in a romics_object by adding a prefix to prevent issues with numeric-only sample names.
#' This function is particularly useful when sample names are purely numeric, which can cause issues in R data processing,
#' plotting, and statistical analysis. The function adds a specified prefix to all sample names across all layers of the romics_object.
#' @param romics_object A romics_object created using createRomicsObject()
#' @param prefix A character string to be added as a prefix to all sample names. Default is "sample"
#' @param separator A character string used to separate the prefix from the original sample name. Default is "_"
#' @param force_adjustment Logical. If TRUE, adds prefix to all sample names regardless of whether they are numeric.
#'        If FALSE (default), only adds prefix to purely numeric sample names.
#' @details This function modifies sample names (column names) across all relevant layers of the romics_object
#' @return Returns the romics_object with adjusted sample names. A message will indicate
#'         how many sample names were modified.
#' @author Geremy Clair
#' @export
romicsNumericHeaderAdjust <- function(romics_object, prefix = "sample", separator = "_", force_adjustment = FALSE) {

  arguments <- as.list(match.call())

  # Input validation
  if (!is.romicsObject(romics_object) || missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }

  if (!is.character(prefix) || length(prefix) != 1 || prefix == "") {
    stop("prefix must be a non-empty character string")
  }

  if (!is.character(separator) || length(separator) != 1) {
    stop("separator must be a character string of length 1")
  }

  if (!is.logical(force_adjustment) || length(force_adjustment) != 1) {
    stop("force_adjustment must be a logical value (TRUE or FALSE)")
  }

  # Get current sample names
  current_names <- colnames(romics_object$data)
  n_samples <- length(current_names)

  if (n_samples == 0) {
    warning("No samples found in romics_object")
    return(romics_object)
  }

  # Determine which samples need adjustment
  if (force_adjustment) {
    # Adjust all sample names
    samples_to_adjust <- rep(TRUE, n_samples)
    message("Force adjustment mode: Adding prefix to all sample names")
  } else {
    # Only adjust purely numeric sample names
    samples_to_adjust <- !is.na(suppressWarnings(as.numeric(current_names)))
    n_numeric <- sum(samples_to_adjust)

    if (n_numeric == 0) {
      message("No purely numeric sample names detected. No changes made.")
      return(romics_object)
    } else {
      message(paste("Found", n_numeric, "purely numeric sample names out of", n_samples, "total samples"))
    }
  }

  # Create new sample names
  new_names <- current_names
  new_names[samples_to_adjust] <- paste0(prefix, separator, current_names[samples_to_adjust])

  # Check for potential duplicates after adjustment
  if (any(duplicated(new_names))) {
    duplicated_names <- new_names[duplicated(new_names)]
    warning(paste("The following sample names will be duplicated after adjustment:",
                  paste(duplicated_names, collapse = ", "),
                  ". Consider using a different prefix."))
  }

  # Apply new names to all relevant layers
  # Update data layer
  colnames(romics_object$data) <- new_names

  # Update metadata layer
  colnames(romics_object$metadata) <- new_names

  # Update missingdata layer
  if (!is.null(romics_object$missingdata)) {
    colnames(romics_object$missingdata) <- new_names
  }

  # Update embeddings layer if it exists
  if (!is.null(romics_object$embeddings)) {
    colnames(romics_object$embeddings) <- new_names
  }

  # Update statistics layer if it exists (column names might reference samples)
  # Note: Statistics layer typically has features as rows, so we don't modify it
  # unless there are sample-specific statistics columns

  # Report changes
  n_changed <- sum(samples_to_adjust)
  if (n_changed > 0) {
    message(paste("Successfully adjusted", n_changed, "sample names"))
    message(paste("Example: '", current_names[which(samples_to_adjust)[1]], "' -> '",
                  new_names[which(samples_to_adjust)[1]], "'", sep = ""))
  }

  # Update processing steps
  romics_object <- romicsUpdateSteps(romics_object, arguments)

  return(romics_object)
}

#' romicsFactorNames()
#' @description Indicates the list of factor names from the romics_object
#' @param romics_object A romics_object created using createRomicsObject()
#' @details This function allows to quickly get a vector containing all the factor names present in a romics_object
#' @return A character vector containing the list of factor contained in an romics_object
#' @author Geremy Clair
#' @export
romicsFactorNames<-function(romics_object){
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  rownames(romics_object$metadata)
}

#' romicsSampleNames()
#' @description Indicates the list of sample names from the romics_object (e.g.,colnames )
#' @param romics_object A romics_object created using createRomicsObject()
#' @details This function allows to quickly get a vector containing all the factor names present in a romics_object
#' @return A character vector containing the list of factor contained in an romics_object
#' @author Geremy Clair
#' @export
romicsSampleNames<-function(romics_object){
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  colnames(romics_object$metadata)
}

#' romicsFeatureNames()
#' @description Indicates the list of feature names from the romics_object (e.g.,rownames )
#' @param romics_object A romics_object created using createRomicsObject()
#' @details This function allows to quickly get a vector containing all the factor names present in a romics_object
#' @return A character vector containing the list of features contained in an romics_object
#' @author Geremy Clair
#' @export
romicsFeatureNames<-function(romics_object){
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  rownames(romics_object$data)
}

#' romicsIDsNames()
#' @description Indicates the list of ID names from the romics_object
#' @param romics_object A romics_object created using createRomicsObject()
#' @details This function allows to quickly get a vector containing all the ID names present in a romics_object
#' @return A character vector containing the list of IDs contained in an romics_object
#' @author Geremy Clair
#' @export
romicsIDsNames<-function(romics_object){
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  colnames(romics_object$IDs)
}

#' romicsExtractFactor()
#' @description Extract a factor from the romics_object
#' @param romics_object A romics_object created using createRomicsObject()
#' @param factor has to be either 'main' to indicate that the main factor has to be used (default) or a factor contained in the romics_object, the list of factors can be obtained using the function romicsFactorNames()
#' @details This function allows to quickly extract the content of a factor present in the romics_object.
#' @return a factor from an romics_object
#' @author Geremy Clair
#' @export
romicsExtractFactor<-function(romics_object, factor = "main"){
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(factor)){
    factor<-"main"
    print("The factor was missing, the 'main' factor of the romics_object was extracted")
  }
  if(factor=="main"){factor<-romics_object$main_factor}
  if(!factor %in% romicsFactorNames(romics_object)){stop("The <factor> is not in the list of factors for this romics_object. The list of available factors can be obtained with the function romicsFactorNames()")}
  fact<-as.factor(t(romics_object$metadata[romicsFactorNames(romics_object)==factor,]))
  names(fact)<-colnames(romics_object$metadata)
  return(fact)
  }

#' romicsExtractFeature()
#' @description Extract the data for a single feature from the romics_object
#' @param romics_object A romics_object created using createRomicsObject()
#' @param feature has to be either a feature contained in the romics_object, the list of features can be obtained using the function romicsFeatureNames()
#' @details This function allows to quickly extract the data for a romics_object feature.
#' @return a feature from an romics_object
#' @author Geremy Clair
#' @export
romicsExtractFeature<-function(romics_object, feature = "feature"){
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(feature)){
    factor<-romicsFeatureNames(romics_object)[1]
    print("The factor was missing, the first of the romics_object was extracted")
  }
  if(!feature %in% romicsFeatureNames(romics_object)){stop("The <feature> is not in the list of features for this romics_object. The list of available features can be obtained with the function romicsFeatureNames()")}
  f<-as.character(t(romics_object$data[romicsFeatureNames(romics_object)==feature,]))
  names(f)<-colnames(romics_object$data)
  return(f)
}

#' romicsUpdateColor()
#' @description Updates the colors layer contained in the romics_object
#' @param romics_object A romics_object created using createRomicsObject()
#' @return This function returns a romics_object with updated colors.
#' @author Geremy Clair
#' @export
romicsUpdateColor<- function(romics_object) {
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}

  #get the colors from the metadata file
  colors_romics<- as.character(t(romics_object$metadata[grepl("colors_romics",rownames(romics_object$metadata)),]))

  #fill <- character(length = 0)
  #for (i in 1:length(colors_romics))
  #{
  #  fill<- c(fill,rep(as.character(colors_romics[i]),nrow(romics_object$data)))
  #}

  romics_object$colors<-rep(as.character(colors_romics), each = nrow(romics_object$data))

  return(romics_object)

}

#' romicsUpdateSteps()
#' @description Updates the steps of the romics_object, require to have recorded the argument in earlier steps of the function
#' @param romics_object A romics_object created using createRomicsObject()
#' @param arguments the arguments of a function are required to read the user input of a function, this user input will be used to generate the steps, the arguments are obtained by running the following code <arguments<-as.list(match.call())> in the first line of a function
#' @details The goal of Romics processor is to provide a trackable and reproducible pipeline for processing omics data. Subsequently it is necessary when a function is created to implement a way to record the user input that will be recorded in the steps layer of the Romics_object.
#' @details This function will enable to simplify the work of developers who want to contribute to Romics by simplifying this process. Only two lines of codes are then necessary to update the steps.
#' @details The first line of code has to be placed in the first line after the function declaration : <arguments<-as.list(match.call())>
#' @details The second line of code has to be <romics_object<-stepUpdater(romics_object,arguments)> placed at the end of the function code (ideally right before returning the processed romics_object or graphic generated by the function)
#' @return This function add the description of the processing to the step layer of an Romics object
#' @author Geremy Clair
#' @export
romicsUpdateSteps<-function(romics_object, arguments){
  if(missing(arguments)){
    steps<- c(paste0("date|",gsub(" ","_",format(Sys.time(),"%b_%d_%Y_%X")),"|step_not_recorded"), "note|The arguments of the function were not recorded using arguments<-as.list(match.call())")
  }
  funName<-arguments[[1]]
  argumentsNames<-names(arguments)[2:length(arguments)]

  if(romics_object$steps[1]!="romics_object"){stop(paste0("The function ",funName," was run on an object that was not an romics_object"))}

  steps<- c(paste0("date|",gsub(" ","_",format(Sys.time(),"%b_%d_%Y_%X")),"|",funName))
  fun<-paste0(funName,"(")
  for(i in 1:length(argumentsNames)){
    fun<-paste0(fun,argumentsNames[i],"=")
    if(typeof(arguments[[i+1]])=="language"){
      arguments[[i+1]]<-gsub("\"","'",deparse(arguments[[i+1]]))
      fun<-paste0(fun,unlist(arguments[[i+1]]))
      }else{
    if(typeof(arguments[[i+1]])=="character"){
      fun<-paste0(fun,"'",unlist(arguments[[i+1]]),"'")
    }else{
      fun<-paste0(fun,unlist(arguments[[i+1]]))
      }}
    if(i<length(argumentsNames)){fun<-paste0(fun,",")}
  }
  fun<-paste0("fun|",fun,")")

  steps<- c(steps,fun)
  romics_object$steps<-c(romics_object$steps,steps)
  return(romics_object)
}



#' romicsLogCheck()
#' @description Identifies if the romics_object is log_transformed or not
#' @param romics_object A romics_object created using createRomicsObject()
#' @return This function will return TRUE or FALSE indicating if the object was or not log transformed using the function log2transform
#' @author Geremy Clair
#' @export
romicsLogCheck<-function(romics_object){
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}

  log<-0
  r<-FALSE
  if((sum(romicsSteps(romics_object,show_dates = F,show_details = F)=="log2transform")-
      sum(romicsSteps(romics_object,show_dates = F,show_details = F)=="unlog2data")) %% 2 == 1 ){log<-log+2}

  if((sum(romicsSteps(romics_object,show_dates = F,show_details = F)=="log10transform")-
      sum(romicsSteps(romics_object,show_dates = F,show_details = F)=="unlog10data")) %% 2 == 1 ){log<-log+10}

  if(log==2|log==10){r<-TRUE}
  if(log>10){stop("The data was logged more than once")}

  return(r)
}

#' romicsCreateDependencies()
#' @description Creates the original dependencies of the romics_object (only when it is created to add dependencies use the romicsAddDependency() function)
#' @return This function will return a data.frame with two columns the required packages and the version at the time the code was run
#' @author Geremy Clair
#' @export
romicsCreateDependencies<-function(){
  # Get package dependencies from DESCRIPTION file
  desc <- packageDescription("RomicsProcessor")
  imports_str <- desc$Imports

  # Parse the Imports string to extract package names
  if(!is.na(imports_str)) {
    packages <- trimws(strsplit(imports_str, ",")[[1]])
    # Remove version specifications (e.g., "package (>= 1.0)" -> "package")
    packages <- gsub("\\s*\\(.*\\)$", "", packages)
  } else {
    packages <- character(0)
  }

  # Create data frame with required packages
  Required <- data.frame(Required = packages, Version_used = NA, stringsAsFactors = FALSE)

  # Get version for each package
  for(i in 1:nrow(Required)){
    tryCatch({
      Required$Version_used[i] <- as.character(packageVersion(Required$Required[i]))
    }, error = function(e) {
      Required$Version_used[i] <<- "unavailable"
    })
  }

  Required[,1] <- as.character(Required[,1])

  # Add R version
  Required <- rbind(Required, c(Required = "r", Version_used = paste0(R.Version()$major, ".", R.Version()$minor)))
  # Add RomicsProcessor itself at the top
  Required <- rbind(c(Required = "RomicsProcessor", Version_used = as.character(packageVersion("RomicsProcessor"))), Required)

  return(Required)
}

#' romicsAddDependency()
#' @description Adds a package to the list of dependencies of the romics_object. Enables developers to add automatically a dependency when their function has been applied by the user on their romics_object.
#' @param romics_object A romics_object created using createRomicsObject()
#' @param new_dependency The name of one or more R package
#' @return This function will return an romics_object updated with the new dependency.
#' @author Geremy Clair
#' @export
romicsAddDependency<-function(romics_object,new_dependency=c("package_1", "package_2")){
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(new_dependency)){stop("<new_dependencies> was empty. No dependencies were added.")}

  installed_packages<-as.data.frame(installed.packages())
  installed_packages<-installed_packages[,colnames(installed_packages) %in% c("Package", "Version")]
  colnames(installed_packages)<-c("Required", "Version_used")

  for(i in 1:length(new_dependency)){
  if(!new_dependency[i] %in% installed_packages$Required){
    stop(paste0("The package ",new_dependency[i]," is not installed it cannot be added to the dependencies of the romics_object."))}
  if(new_dependency[i] %in% romics_object$dependencies$Required){
    warning(paste0("The package ",new_dependency[i]," was in the list of dependencies of the romics_object, the version_used was updated"))
    romics_object$dependencies<-romics_object$dependencies[!romics_object$dependencies$Required %in% new_dependency[i],]
    }
    romics_object$dependencies<-rbind(romics_object$dependencies,installed_packages[installed_packages$Required==new_dependency[i],])
    }

 return(romics_object)
}

#' romicsCalculatedStats()
#' @description Indicates the stat columns calculated for the romics_object
#' @param romics_object A romics_object created using createRomicsObject()
#' @return This function will return character vector containing the stat columns previously calculated for a romics_object, if no stats were previously calculated an error message will be displayed
#' @author Geremy Clair
#' @export
romicsCalculatedStats<-function(romics_object){
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(is.null(romics_object$statistics)){
   warning("No statistics were calculated for this romics_object")
    return(FALSE)
     }else{
      return(colnames(romics_object$statistics))}
  }

#' romicsSteps()
#' @description Displays the content steps layer of the romics_object
#' @param romics_object A romics_object created using createRomicsObject()
#' @param show_dates Boolean indicating if the dates have to be displayed
#' @param show_details Boolean indicating if the details have to be displayed
#' @return This function will return the steps of an romics_object
#' @export
romicsSteps<-function(romics_object, show_dates=TRUE, show_details=TRUE){
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(show_dates)){show_dates=TRUE}
  if(!show_dates %in% c(TRUE,FALSE)){stop("show_dates has to be either TRUE or FALSE")}
  if(missing(show_details)){show_details=TRUE}
  if(!show_details %in% c(TRUE,FALSE)){stop("show_details has to be either TRUE or FALSE")}

  steps<-romics_object$steps[2:length(romics_object$steps)]
  if(show_dates==FALSE){steps[grepl("date\\|",steps)]<- gsub(".*\\|","", steps[grepl("date\\|",steps)])}
  if(show_details==FALSE){steps<-steps[!grepl("fun\\|",steps)]}
  return (steps)
}

#' createRomicsPipeline()
#' @description Extracts a pipeline from the romics_object. the pipeline can then be saved in a classic R object.
#' @param romics_object A romics_object created using createRomicsObject()
#' @return This function will return character vector containing the stat columns previously calculated for a romics_object, if no stats were previously calculated an error message will be displayed
#' @author Geremy Clair
#' @export
createRomicsPipeline<- function(romics_object){
  arguments<-as.list(match.call())
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}

  #create the text of the pipeline
  pl<- romics_object$steps
  #remove the romics creation steps
  pl<- pl[4:length(pl)]
  #keep only the elements with the recorded details
  pl<- pl[substr(pl, 1, 4)=="fun|"]
  #remove fun|
  pl<-sub("fun\\|","",pl)

  for(i in 1: length(pl)){
    if(grepl("\\,",pl[i])){
     before_parenthesis<- gsub("\\(.*","(",pl[i])
     after_first_comma <- sub("^.*?\\,",",",pl[i])
     pl[i]<-paste0(before_parenthesis,"|@|",after_first_comma)
      }else{pl[i]<-gsub("romics_object\\=[^>]+\\)" , "|@|)",pl[i])}
  }
  pl<-paste0("|@|<-",pl)

  pl<-gsub("\\|\\@\\|", "romics_object",pl)

  return(pl)

  }

#' applyRomicsPipeline()
#' @description Applies a pipeline created with the romicsCreatePipeline() function, pipelines can be edited prior to run them.
#' @param romics_object A romics_object created using createRomicsObject()
#' @param romics_pipeline A pipeline created with the romicsCreatePipeline() function.
#' @return This function will return an romics_object that has been processed through a pipeline
#' @export
applyRomicsPipeline<- function(romics_object, romics_pipeline){
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(typeof(romics_pipeline)!= "character" |sum(substr(romics_pipeline, 1, 15)=="romics_object<-")!=length(romics_pipeline)){stop("Your pipeline will not be applied to the romics_object. Please check the text of your pipeline")}

  for(i in 1:length(romics_pipeline)){
    eval(parse(text=romics_pipeline[i]))
  }

  return(romics_object)
}

#' romicsExportData()
#' @description Creates an exportable data frame from the romics_object. The generated data.frame contains the processed data, the statistics and the missing status of the data on demand.
#' @param romics_object A romics_object created using createRomicsObject()
#' @param statistics boolean, has to be TRUE or FALSE. Indicates if the statistics should be exported along with the data (FALSE by default).
#' @param missing_data boolean, has to be TRUE or FALSE. Indicates if the missing status of the data should be exported along with the data (FALSE by default).
#' @param IDs boolean, has to be TRUE or FALSE. Indicates if the IDs of the data should be exported along with the data (FALSE by default).
#' @return This function will return an data frame containing the results of the processing, the statistics and the missingness status of the data as specified by the user.
#' @author Geremy Clair
#' @export
romicsExportData<-function(romics_object, data=TRUE, statistics = FALSE, missing_data = FALSE, IDs=FALSE){
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(statistics)){statistics = FALSE}
  if(missing(missing_data)){missing_data = FALSE}
  if(!is.logical(statistics)){stop("<statistics> should be either TRUE or FALSE")}
  if(!is.logical(missing_data)){stop("<missing_data> should be either TRUE or FALSE")}
  if(!is.logical(IDs)){stop("<IDs> should be either TRUE or FALSE")}

  df<-romics_object$data

  if(data==FALSE){
  df<-df[,-1:-ncol(df)]
  }

  if(statistics==TRUE){
    if(is.null(romics_object$statistics)){stop("The selected romics object does not contain a 'statistics' layer")}else{df<-cbind(df,romics_object$statistics)}
    }

  if(missing_data==TRUE){
    md <- romics_object$missingdata
    colnames(md)<-paste0("missing_data_",colnames(md))
    df<-cbind(df,md)
    }

  if(IDs==TRUE){
    if("IDs" %in% names(romics_object)){
      ids <- romics_object$IDs
      ids <- cbind(original_ids=rownames(ids),ids)
      df<-cbind(original_ids=rownames(df),df)
      df<- merge(df,ids,by="original_ids")
      rownames(df)<-df[,colnames(df)=="original_ids"]
      df<-df[,colnames(df)!="original_ids"]
    }
  }

  return(df)
}
