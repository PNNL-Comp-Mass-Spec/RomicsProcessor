#' log2transform()
#' @description log2-tranforms the romics_object data layer.
#' @param romics_object has to be an romics_object created using romicsCreateObject()
#' @details will log2 transform the romics object
#' @return This function returns the transformed romics_object with updated data layer
#' @author Geremy Clair
#' @export
#'
log2transform<-function(romics_object){
  arguments<-as.list(match.call())
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(romicsLogCheck(romics_object)){stop("Your romics_object was previously log transformed, see your romics_object$steps layer for more details")}
  romics_object$data<- data.frame(apply(romics_object$data,2,log2))
  romics_object<-romicsUpdateSteps(romics_object,arguments)
  return(romics_object)
  }

#' log10transform()
#' @description log10-tranforms the romics_object data layer.
#' @param romics_object has to be an romics_object created using romicsCreateObject()
#' @details will log10 transform the romics object
#' @return This function returns the transformed romics_object with updated data layer
#' @author Geremy Clair
#' @export
log10transform<-function(romics_object){
  arguments<-as.list(match.call())
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(romicsLogCheck(romics_object)){stop("Your romics_object was previously log transformed, see the romics_object$steps layer for more details")}
  romics_object$data<- data.frame(apply(romics_object$data,2,log10))
  romics_object<-romicsUpdateSteps(romics_object,arguments)
  return(romics_object)
  }

#' unlog2data()
#' @description Reverses the log2 tranformation of the romics_object data layer.
#' @param romics_object has to be an romics_object created using romicsCreateObject()
#' @details will reverse the log2 transformation of the romics_object
#' @return This function returns the transformed romics_object with updated data layer
#' @author Geremy Clair
#' @export
unlog2data<-function(romics_object){
  arguments<-as.list(match.call())
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(!sum(grepl("log2transform\\(",romics_object$steps))>0){
    stop("Your romics_object was not previously log2 transformed using the function log2transform(). See the romics_object$steps layer for more details.")}
  if(sum(grepl("log2transform\\(",romics_object$steps))<=sum(grepl("unlog2data\\(",romics_object$steps))){
   stop("Your romics_object was already unlogged using the function unlog2data(). See the romics_object$steps layer for more details.")}
  romics_object$data<- 2^romics_object$data
  romics_object<-romicsUpdateSteps(romics_object,arguments)
  return(romics_object)
}

#' unlog10data()
#' @description Reverses the log10 tranformation of the romics_object data layer.
#' @param romics_object has to be an romics_object created using romicsCreateObject()
#' @details will reverse the log10 transformation of the romics_object
#' @return This function returns the transformed romics_object with updated data layer
#' @author Geremy Clair
#' @export
unlog10data<-function(romics_object){
  arguments<-as.list(match.call())
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(!sum(grepl("log10transform\\(",romics_object$steps))>0){
    stop("Your romics_object was not previously log10 transformed using the function log10transform(). See the romics_object$steps layer for more details.")}
  if(sum(grepl("log10transform\\(",romics_object$steps))<=sum(grepl("unlog10data\\(",romics_object$steps))){
    stop("Your romics_object was already unlogged using the function unlog10data(). See the romics_object$steps layer for more details.")}
  romics_object$data<- 10^romics_object$data
  romics_object<-romicsUpdateSteps(romics_object,arguments)
  return(romics_object)
  }

#' medianNormSample()
#' @description Normalizes the samples by their median. The median of the medians of all the samples is used as the alignment point.
#' @param romics_object has to be an romics_object created using romicsCreateObject()
#' @details Median normalize within each sample the median of each sample median will be used as alignment point. If you waht to center the median at 0 please use the function medianCenterSample().
#' @return This function returns the transformed romics_object with updated data layer.
#' @author Geremy Clair
#' @export
medianNormSample<-function(romics_object){
  arguments<-as.list(match.call())
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  col_median<- apply(romics_object$data,2,function(x){median(x[!is.na(x)])})
  median_median<- median(col_median)
  f<-function(row){
    if(romicsLogCheck(romics_object)==TRUE){
      row<-row-col_median+median_median}else{
      row<-row/col_median*median_median}
    return(row)}
  romics_object$data<- as.data.frame(t(apply(romics_object$data,1,f)))
  romics_object<-romicsUpdateSteps(romics_object,arguments)
  return(romics_object)
  }

#' medianCenterSample()
#' @description Normalizes the samples by their median. Zero is used as the median alignment center.
#' @param romics_object has to be an romics_object created using romicsCreateObject()
#' @details Median normalize within each sample the median will be zero centered
#' @return This function returns the transformed romics_object with updated data layer.
#' @author Geremy Clair
#' @export
medianCenterSample<-function(romics_object){
  arguments<-as.list(match.call())
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  col_median<- apply(romics_object$data,2,function(x){median(x[!is.na(x)])})
  f<-function(row){
    if(romicsLogCheck(romics_object)==TRUE){
      row<-row-col_median}else{
        row<-row/col_median}
    return(row)}
  romics_object$data<- as.data.frame(t(apply(romics_object$data,1,f)))
  romics_object<-romicsUpdateSteps(romics_object,arguments)
  return(romics_object)
}

#' medianCenterFeature()
#' @description Normalizes the features by their median. Zero is used as the median alignment center.
#' @param romics_object has to be an romics_object created using romicsCreateObject()
#' @details Median normalize within each feature so that the median will be zero centered
#' @return This function returns the transformed romics_object with updated data layer.
#' @author Geremy Clair
#' @export
medianCenterFeature <- function(romics_object) {
  arguments <- as.list(match.call())

  if (!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }

  row_median <- apply(romics_object$data, 1, function(x) { median(x[!is.na(x)]) })

  f <- function(column) {
    if (romicsLogCheck(romics_object) == TRUE) {
      column <- column - row_median
    } else {
      column <- column / row_median
    }
    return(column)
  }

  romics_object$data <- as.data.frame(apply(romics_object$data, 2, f))
  romics_object <- romicsUpdateSteps(romics_object, arguments)

  return(romics_object)
}

#' medianNormFactor()
#' @description Normalizes the samples by their median within a given factor. The median of the median within this factor will be used as factor-specific median center.
#' @param romics_object has to be an romics_object created using romicsCreateObject()
#' @param factor Character. Factor name from romics_object to use for grouping. Use "main" for the main factor. Default: "main"
#' @details median normalize the samples within a given factor the median of the median of this factor will be used as new factor median.
#' @return This function returns the transformed romics_object with updated data layer.
#' @author Geremy Clair
#' @export
medianNormFactor<-function(romics_object, factor= "main"){
  if(missing(romics_object)){stop("romics_object is missing")}
  if(!inherits(romics_object, "romics_object")){stop("your romics_object was not created using the function romicsCreateObject")}
  if(romics_object$steps[1]!="romics_object"){stop("romics_object is not in the appropriate format")}
  if(missing(factor)){factor<-"main"}
  # Resolve factor name
  if(factor == "main"){
    factor_name <- romics_object$main_factor
  } else {
    factor_name <- factor
  }
  #import the data
  data<-romics_object$data
  #calculate the median for each column
  col_median<- apply(data,2,function(x){median(x[!is.na(x)])})
  #unique factors
  factor_values<- as.character(t(romics_object$metadata[rownames(romics_object$metadata)==factor_name,]))
  unique_factor<-unique(factor_values)
  #calculate the median of the median within each factor
  median_factor<-as.numeric()
  for (i in 1:length(unique_factor)){
    median_factor[i]<-median(col_median[factor_values==unique_factor[i]])
  }
  #calculate the adjusting of the median for each sample (based on a given factor)
  adjust_median <- median_factor[match(factor_values,unique_factor)]
  #do the normalization using the - option if the data was log transformed or with the / if the data was not transformed
  if(sum(romics_object$steps=="log10transform()")+sum(romics_object$steps=="log2transform()")==0){
    print("Your data was median transformed")
    for(i in 1:ncol(data)){
      data[,i]<- data[,i]/col_median[i]*adjust_median[i]
    }
  } else {
    print("Your log transformed data was median transformed")
    for(i in 1:ncol(data)){
      data[,i]<- data[,i]-col_median[i]+adjust_median[i]
    }
  }
  #place this in your romics_object
  romics_object$data<-data
  #append the steps
  romics_object$steps<- c(romics_object$steps,paste0(gsub(" ","_",format(Sys.time(),"%b_%d_%Y_%X")),": the data was was median normalized for each column"))
  romics_object$steps<- c(romics_object$steps,"medianNormFactor()")
  #return object
  return(romics_object)
  }

#' quantileNormSample()
#' @description Quantile-normalizes the samples.
#' @param romics_object has to be an romics_object created using romicsCreateObject()
#' @details Quantile-normalizes within each sample.
#' @return This function returns the transformed romics_object with updated data layer.
#' @author Geremy Clair
#' @export
quantileNormSample<-function(romics_object){
  arguments<-as.list(match.call())
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  d <- romics_object$data
  original_rownames <- rownames(d)

  r <- apply(d,2,rank,ties.method="min")

  s<-data.frame(apply(d, 2, sort, na.last=TRUE))

  m <- rowMeans(s,na.rm=T)

  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }

  t <- apply(r, 2, index_to_mean, my_mean=m)
  t[is.na(d)]<-NA

  # Convert to data frame and restore original rownames
  normalized_data <- data.frame(t)
  rownames(normalized_data) <- original_rownames

  romics_object$data<-normalized_data

  romics_object<-romicsUpdateSteps(romics_object,arguments)

  return(romics_object)
}

#' upperQuartileNormSample()
#' @description This function Normalizes the data using the UpperQuartile value (i.e., third quartile).
#' @param romics_object has to be an romics_object created using romicsCreateObject() that has not been previously log-transformed using log10transform()
#' @details UpperQuartile-normalizes the dataset.
#' @return This function returns the transformed romics_object with updated data layer.
#' @author Geremy Clair
#' @export
upperQuartileNormSample<-function(romics_object){
  arguments <- as.list(match.call())
  if (!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}

  d <- romics_object$data
  upper_quartile <- apply(d, 2, function(x) quantile(x, 0.75, na.rm = TRUE))
  normalized_data <- sweep(d, 2, upper_quartile, FUN = "/", check.margin = FALSE)
  romics_object$data <- data.frame(normalized_data)
  romics_object <- romicsUpdateSteps(romics_object, arguments)

  return(romics_object)
}

#' romicsMeanByLevel()
#' @description The function will condensate the values based on a romics factor, it will reduce the number of sample to one by level using the mean for each feature.
#' @param romics_object has to be an romics_object created using romicsCreateObject()
#' @param averaging_factor as to be a factor of an romics_object factors can be identified using the function romicsFactorNames()
#' @param main_factor as to be either "main" or a factor of the romics_object that can be identified using the function romicsFactorNames()
#' @details The function will average the data inside the factor, NA will be not considered. a new metadata file will be created, factors with discrepancies while the averaging is done will be removed.
#' @return This function returns the transformed romics_object with updated data layer
#' @author Geremy Clair, Igor Estevao
#' @export
#'
romicsMeanByLevel<-function(romics_object, averaging_factor, main_factor){
  arguments <- as.list(match.call())
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("The <romics_object> is missing or is not in the appropriate format")}
  if(!averaging_factor %in% romicsFactorNames(romics_object) | missing(averaging_factor)){stop("The <averaging_factor> is missing or is not a factor of your romics_object, use the function romicsFactorNames() to identify the factors")}
  if(missing(main_factor)){main_factor<-"main"}
  if(main_factor=="main"){main_factor<-romics_object$main_factor}
  if(!main_factor %in% romicsFactorNames(romics_object)){stop("The <main_factor> is not a factor of your romics_object, use the function romicsFactorNames() to identify the factors")}

  data <-data.frame(romics_object$data)
  missingdata<-data.frame(romics_object$missingdata)
  metadata<-romics_object$metadata
  f<-romicsExtractFactor(romics_object,factor = averaging_factor)

  resulting_data<-data.frame(matrix(nrow=nrow(data),ncol=0))
  resulting_missingdata<-data.frame(matrix(nrow=nrow(missingdata),ncol=0))
  resulting_metadata<-data.frame(matrix(nrow=nrow(metadata),ncol=0))
  rownames(resulting_data)<- rownames(resulting_missingdata)<-rownames(data)
  rownames(resulting_metadata)<- rownames(metadata)
  for(i in 1:length(levels(f))){
    l<-as.matrix(data[,f==levels(f)[i]])
    t<-data.frame(t=as.numeric(l[,1]))
    lm<-data.frame(missingdata[,f==levels(f)[i]])
    tm<-data.frame(t=as.character(lm[,1]))
    lmeta<-as.matrix(metadata[,f==levels(f)[i]])
    tmeta<-data.frame(t=as.character(lmeta[,1]))
    for(j in 1:nrow(t)){
      if(sum(is.na(l[j,]))==ncol(l)){
        t[j,1]<-NA
      }else{
        t[j,1]<-mean(l[j,],na.rm=T)
      }
      if(sum(lm[j,])==ncol(lm)){
        tm[j,1]<-TRUE
      }else{
        tm[j,1]<-FALSE
      }
    }

    for (j in 1:nrow(lmeta)){
      if(sum(lmeta[j,]==lmeta[j,1])==ncol(lmeta)){
        tmeta[j,1]<-lmeta[j,1]
      }else{
        tmeta[j,1]<-paste(lmeta[j,],collapse = "@")
      }
    }

    colnames(t)<- colnames(tm)<-colnames(tmeta)<-paste0(levels(f)[i],"_mean")
    rownames(t)<-rownames(tm)<-rownames(l)
    rownames(tmeta)<-rownames(metadata)
    resulting_data<-cbind(resulting_data,t)
    resulting_missingdata<-cbind(resulting_missingdata,tm)
    resulting_metadata<-cbind(resulting_metadata,tmeta)
  }
  romics_object$data<-resulting_data
  romics_object$missingdata<-resulting_missingdata
  romics_object$metadata<-resulting_metadata
  romics_object<- romicsChangeFactor(romics_object,main_factor = main_factor)

  if("statistics" %in% names(romics_object)){
    warning("The statistics layer was removed from the romics_object the statistics were calculated on the non-subsetted object")
    romics_object<-romics_object["statistics" != names(romics_object)]}

    romics_object <- romicsUpdateSteps(romics_object, arguments)
    class(romics_object)<-"romics_object"
    return(romics_object)
  }

#' createBulkMetadata()
#' @description Create enhanced pseudobulk metadata with optimized processing
#' @param romics_object Original romics object
#' @param bulk_factor Factor values used for bulking (as factor or character vector)
#' @param bulk_levels Unique levels of the bulking factor
#' @param sample_names Names for the new pseudobulk samples
#' @param keep_factors Character vector of factors to preserve in metadata. If NULL, preserves all factors.
#' @param bulk_by Name of the bulking factor
#' @return Enhanced metadata data frame with optimized processing
#' @author Geremy Clair
#' @export
createBulkMetadata <- function(romics_object, bulk_factor, bulk_levels, sample_names, keep_factors = NULL, bulk_by) {

  # Input validation
  if (!is.romicsObject(romics_object)) {
    stop("romics_object must be a valid romics object")
  }

  if (length(bulk_levels) != length(sample_names)) {
    stop("bulk_levels and sample_names must have the same length")
  }

  # Get all available factors
  base_factors <- romicsFactorNames(romics_object)

  # Validate bulk_by factor
  if (!bulk_by %in% base_factors) {
    stop(paste("bulk_by factor '", bulk_by, "' not found in romics object", sep = ""))
  }

  # Validate keep_factors if specified
  if (!is.null(keep_factors)) {
    invalid_factors <- keep_factors[!keep_factors %in% base_factors]
    if (length(invalid_factors) > 0) {
      warning(paste("The following factors in keep_factors were not found:", paste(invalid_factors, collapse = ", ")))
      keep_factors <- keep_factors[keep_factors %in% base_factors]
    }
  }

  # Define factors to process
  factors_to_process <- if (is.null(keep_factors)) base_factors else keep_factors

  # Add special metadata factors
  new_factors <- c("original_factor", "n_observations")
  all_factors <- c(factors_to_process, new_factors)

  # Pre-allocate metadata matrix for efficiency
  bulk_metadata <- matrix(NA_character_,
                          nrow = length(all_factors),
                          ncol = length(sample_names),
                          dimnames = list(all_factors, sample_names))

  # Convert bulk_factor to character for consistent processing
  bulk_factor_char <- as.character(bulk_factor)

  # Fill special metadata first (most efficient operations)
  bulk_metadata["original_factor", ] <- as.character(bulk_levels)

  # Calculate group sizes efficiently using table
  group_counts <- table(bulk_factor_char)
  bulk_metadata["n_observations", ] <- as.character(group_counts[bulk_levels])

  # Pre-calculate group indices for efficiency
  group_indices_list <- vector("list", length(bulk_levels))
  names(group_indices_list) <- bulk_levels

  for (i in seq_along(bulk_levels)) {
    group_indices_list[[i]] <- which(bulk_factor_char == bulk_levels[i])
  }

  # Process each factor efficiently
  for (factor_name in factors_to_process) {

    if (factor_name == bulk_by) {
      # For the bulking factor, just use the levels directly
      bulk_metadata[factor_name, ] <- as.character(bulk_levels)

    } else {
      # Extract factor values once
      factor_values <- romicsExtractFactor(romics_object, factor = factor_name)
      factor_values_char <- as.character(factor_values)

      # Process each group
      for (i in seq_along(bulk_levels)) {
        group_indices <- group_indices_list[[i]]

        if (length(group_indices) > 0) {
          group_values <- factor_values_char[group_indices]

          # Remove NA values for cleaner processing
          group_values_clean <- group_values[!is.na(group_values) & group_values != ""]

          if (length(group_values_clean) > 0) {
            # Find most common value efficiently
            value_counts <- table(group_values_clean, useNA = "no")

            if (length(value_counts) > 0) {
              most_common <- names(value_counts)[which.max(value_counts)]
              bulk_metadata[factor_name, i] <- most_common

              # Warn if there's high heterogeneity in the group
              if (length(value_counts) > 1) {
                max_prop <- max(value_counts) / sum(value_counts)
                if (max_prop < 0.7) {  # Less than 70% consensus
                  warning(paste("High heterogeneity in factor '", factor_name,
                                "' for group '", bulk_levels[i],
                                "' (", round(max_prop * 100, 1), "% consensus)", sep = ""))
                }
              }
            }
          }
        }
      }
    }
  }

  # Convert to data frame with proper structure
  bulk_metadata_df <- as.data.frame(bulk_metadata, stringsAsFactors = FALSE)

  # Convert numeric-looking columns back to numeric where appropriate
  for (factor_name in rownames(bulk_metadata_df)) {
    values <- bulk_metadata_df[factor_name, ]
    # Check if all non-NA values can be converted to numeric
    non_na_values <- values[!is.na(values)]
    if (length(non_na_values) > 0) {
      numeric_test <- suppressWarnings(as.numeric(non_na_values))
      if (!any(is.na(numeric_test))) {
        bulk_metadata_df[factor_name, ] <- as.character(suppressWarnings(as.numeric(values)))
      }
    }
  }

  # Add summary information as attributes
  attr(bulk_metadata_df, "bulk_summary") <- list(
    original_samples = length(bulk_factor),
    bulk_samples = length(sample_names),
    bulk_by = bulk_by,
    factors_processed = factors_to_process,
    group_sizes = as.list(group_counts[bulk_levels])
  )

  message(paste("Metadata created:", length(bulk_factor), "->", length(sample_names),
                "samples across", length(factors_to_process), "factors"))

  return(bulk_metadata_df)
}

#' romicsPseudoBulk
#' @description Creates pseudobulk samples by aggregating observations within groups
#' @param romics_object A romics_object created using romicsCreateObject()
#' @param bulk_by Factor name to aggregate by
#' @param method Aggregation method: "mean", "median", or "sum". Default: "mean"
#' @param min_observations Minimum observations required per group. Default: 3
#' @param keep_factors Character vector of factor names to preserve in metadata. Default: NULL
#' @param suffix Suffix to add to sample names. Default: "_bulk"
#' @return A new romics_object with pseudobulk data
#' @author Geremy Clair
#' @export
romicsPseudoBulk <- function(romics_object,
                             bulk_by,
                             method = "mean",
                             min_observations = 3,
                             keep_factors = NULL,
                             suffix = "_bulk") {

  arguments <- as.list(match.call())

  # Basic validation
  if (!is.romicsObject(romics_object) || missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }

  if (missing(bulk_by)) {
    stop("bulk_by is required")
  }

  if (!method %in% c("mean", "median", "sum")) {
    stop("method must be 'mean', 'median', or 'sum'")
  }

  # Handle main factor
  if (bulk_by == "main") {
    bulk_by <- romics_object$main_factor
    message(paste("Using main factor:", bulk_by))
  }

  if (!bulk_by %in% romicsFactorNames(romics_object)) {
    stop(paste("Factor", bulk_by, "not found in metadata"))
  }

  # Extract data as matrices
  data_matrix <- as.matrix(romics_object$data)
  missing_matrix <- as.matrix(romics_object$missingdata)

  # Extract bulking factor
  bulk_factor <- romicsExtractFactor(romics_object, factor = bulk_by)
  bulk_levels <- levels(bulk_factor)
  n_groups <- length(bulk_levels)

  message(paste("Creating pseudobulk:", ncol(data_matrix), "->", n_groups, "samples using", method))

  start_time <- Sys.time()

  # Initialize output matrices
  sample_names <- paste0(bulk_levels, suffix)
  bulk_data <- matrix(NA,
                      nrow = nrow(data_matrix),
                      ncol = n_groups,
                      dimnames = list(rownames(data_matrix), sample_names))

  bulk_missing <- matrix(FALSE,
                         nrow = nrow(data_matrix),
                         ncol = n_groups,
                         dimnames = list(rownames(data_matrix), sample_names))

  # Process each group
  for (i in seq_along(bulk_levels)) {
    level <- bulk_levels[i]
    group_idx <- which(bulk_factor == level)

    if (length(group_idx) < min_observations) {
      warning(paste("Group", level, "has only", length(group_idx), "observations"))
    }

    # Extract group data
    group_data <- data_matrix[, group_idx, drop = FALSE]
    group_missing <- missing_matrix[, group_idx, drop = FALSE]

    # Set missing values to NA
    group_data[group_missing] <- NA

    # Fast aggregation
    if (method == "mean") {
      bulk_data[, i] <- rowMeans(group_data, na.rm = TRUE)
    } else if (method == "median") {
      bulk_data[, i] <- apply(group_data, 1, median, na.rm = TRUE)
    } else if (method == "sum") {
      bulk_data[, i] <- rowSums(group_data, na.rm = TRUE)
    }

    # Check for insufficient data
    valid_counts <- rowSums(!is.na(group_data))
    insufficient <- valid_counts < min_observations
    bulk_data[insufficient, i] <- NA
    bulk_missing[insufficient, i] <- TRUE
  }

  end_time <- Sys.time()
  message(paste("Completed in", round(as.numeric(end_time - start_time, units = "secs"), 1), "seconds"))

  # Get metadata
  original_metadata <- romics_object$metadata

  # Initialize new metadata with sample names
  bulk_metadata <- data.frame(matrix(nrow = 0, ncol = n_groups))
  colnames(bulk_metadata) <- sample_names

  # Add the bulking factor
  bulk_factor_row <- data.frame(matrix(bulk_levels, nrow = 1))
  colnames(bulk_factor_row) <- sample_names
  rownames(bulk_factor_row) <- bulk_by
  bulk_metadata <- rbind(bulk_metadata, bulk_factor_row)

  # Add kept factors if specified
  if (!is.null(keep_factors)) {
    for (factor_name in keep_factors) {
      if (factor_name %in% rownames(original_metadata)) {
        # For each bulk level, get the most common value of this factor
        factor_row <- character(n_groups)

        for (i in seq_along(bulk_levels)) {
          level <- bulk_levels[i]
          group_idx <- which(bulk_factor == level)

          if (length(group_idx) > 0) {
            # Get values for this group
            group_values <- as.character(original_metadata[factor_name, group_idx])
            # Remove NAs
            group_values <- group_values[!is.na(group_values)]

            if (length(group_values) > 0) {
              # Get most common value
              factor_table <- table(group_values)
              most_common <- names(factor_table)[which.max(factor_table)]
              factor_row[i] <- most_common
            } else {
              # No valid values - use NA or "Unknown"
              factor_row[i] <- "Unknown"
              warning(paste("No valid values for factor", factor_name, "in group", level))
            }
          } else {
            factor_row[i] <- "Unknown"
            warning(paste("No observations in group", level))
          }
        }

        factor_df <- data.frame(matrix(factor_row, nrow = 1))
        colnames(factor_df) <- sample_names
        rownames(factor_df) <- factor_name
        bulk_metadata <- rbind(bulk_metadata, factor_df)
      } else {
        warning(paste("Factor", factor_name, "not found in original metadata"))
      }
    }
  }

  bulk_romics <- list()

  # CORE DATA LAYERS (required)
  bulk_romics$data <- as.data.frame(bulk_data)
  bulk_romics$metadata <- bulk_metadata
  bulk_romics$missingdata <- bulk_missing

  # MAIN FACTOR
  if (!is.null(keep_factors) && length(keep_factors) > 0) {
    # Use first kept factor as main factor
    bulk_romics$main_factor <- keep_factors[1]
  } else {
    # Use bulking factor as main factor
    bulk_romics$main_factor <- bulk_by
  }

  # IDs LAYER
  if (!is.null(romics_object$IDs)) {
    bulk_romics$IDs <- romics_object$IDs
  } else {
    bulk_romics$IDs <- data.frame(
      names = rownames(bulk_data),
      stringsAsFactors = FALSE
    )
  }

  # COLORS (will be generated)
  bulk_romics$colors <- NULL

  # STEPS (ensure compliant format)
  bulk_romics$steps <- list("romics_object")

  if (!is.null(romics_object$steps) && length(romics_object$steps) > 0) {
    previous_steps <- romics_object$steps
    if (previous_steps[[1]] != "romics_object") {
      bulk_romics$steps <- c(bulk_romics$steps, previous_steps)
    } else {
      bulk_romics$steps <- previous_steps
    }
  }

  # PRESERVE ORIGINAL DATA
  if (!is.null(romics_object$original_data)) {
    bulk_romics$original_data <- romics_object$original_data
  } else {
    bulk_romics$original_data <- romics_object$data
  }

  if (!is.null(romics_object$original_metadata)) {
    bulk_romics$original_metadata <- romics_object$original_metadata
  } else {
    bulk_romics$original_metadata <- romics_object$metadata
  }

  # PRESERVE OTHER LAYERS
  bulk_romics$dependencies <- romics_object$dependencies
  bulk_romics$custom_colors <- romics_object$custom_colors
  bulk_romics$omics_type <- romics_object$omics_type
  bulk_romics$omics_information <- romics_object$omics_information

  # UUID
  bulk_romics$uuid <- NULL

  # REMOVE INVALID LAYERS
  if (!is.null(romics_object$embeddings)) {
    message("Removing embeddings layer - no longer valid after pseudobulk aggregation")
  }
  if (!is.null(romics_object$statistics)) {
    message("Removing statistics layer - no longer valid after pseudobulk aggregation")
  }

  # SET CLASS
  class(bulk_romics) <- "romics_object"

  # ENSURE UUID EXISTS
  bulk_romics <- romicsAttributeUUID(bulk_romics)

  # UPDATE STEPS
  bulk_romics <- romicsUpdateSteps(bulk_romics, arguments)

  # GENERATE COLORS
  tryCatch({
    bulk_romics <- romicsChangeFactor(bulk_romics, main_factor = bulk_romics$main_factor)
    message("Colors generated for pseudobulk samples")
  }, error = function(e) {
    message("Note: Could not generate colors. Use romicsChangeFactor() if needed.")
  })

  # Summary message
  message(paste("\nPseudobulk completed:", ncol(data_matrix), "->", ncol(bulk_data), "samples"))
  message(paste("Method:", method))
  message(paste("Min observations:", min_observations))
  message(paste("Main factor:", bulk_romics$main_factor))
  if (!is.null(keep_factors)) {
    message(paste("Kept factors:", paste(keep_factors, collapse = ", ")))
  }

  return(bulk_romics)
}
