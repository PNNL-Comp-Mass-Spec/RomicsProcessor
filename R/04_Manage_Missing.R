#' romicsZeroToMissing()
#' @description Replaces zeros in the romics_object by NA values
#' @param romics_object has to be an romics_object created using romicsCreateObject()
#' @details This function will convert 0 values to NA in the data and missingdata layers
#' @return This function returns the transformed romics_object with updated data and missingdata layers
#' @author Geremy Clair
#' @export
romicsZeroToMissing<-function(romics_object){
  arguments<-as.list(match.call())
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  coln<-colnames(romics_object$data)
  rown<-rownames(romics_object$data)
  romics_object$data <- data.frame(
    lapply(romics_object$data, function(x) {
      x <- as.numeric(as.character(x))
      ifelse(x == 0, NA, x)}))
  colnames(romics_object$data)<-coln
  rownames(romics_object$data)<-rown
  romics_object$missingdata<-data.frame(is.na(romics_object$data))
  romics_object<-romicsUpdateSteps(romics_object,arguments)
  return(romics_object)
  }

#' romicsFilterMissing()
#' @description Filters out the variables of the romics_object below the choosen percentage of completeness. The percentage of completeness can either be global or by factor of a given factor (in this later case if the percentage of completeness is achieved for at least one level of the factor it the variable will be kept).
#' @param romics_object A romics_object created using romicsCreateObject()
#' @param percentage_completeness Numerical vector indicating the minimum percentage of data to be considered
#' @param factor has to be either "main", "none" or a factor of an romics_object created using romicsCreateObject() the list of factors can be obtained by running the function romicsFactorNames() on the romics_object
#' @param all_groups if this parameter is TRUE the completeness requirement is for each and every group, if not, the completeness requirement is for at least one group.
#' @details  This function will use the completeness of the features in the overall samples (when none is used as factor), or of a given level of a specific defined factor (in this case the factor has to be set to either "main" or to the given factor of filtering). By default factor is the main factor of the object, the percentage_completeness is set at 50%
#' @return  The function will return a filtered romics_object with the rows of the data and missing data object removed when appropriate.
#' @author Geremy Clair, Nicholas Day
#' @export
romicsFilterMissing<-function(romics_object, percentage_completeness=50, factor = "main", all_groups=FALSE){
  arguments<-as.list(match.call())
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(factor)){
    warning("your factor was missing the main factor of the romics_object was used")
    factor<-"main"}
  if(missing(percentage_completeness)){
    warning("The acceptable percentage of completeness was not set by the user it was set at 50% by default")
    percentage_completeness<-50}
  if(percentage_completeness<0 & percentage_completeness>100){stop("the completeness has to be comprised between 0 and 100 %")}
  if(sum(is.na(romics_object$data))==0){warning("There is no missing values in this dataset data to be removed")}
  if(missing(all_groups)){all_groups=FALSE}
  #extract factor
  if(factor=="none"){
    selected_factor<-rep("overal_sample_number",ncol(romics_object$metadata))}
  if(factor %in% romicsFactorNames(romics_object) ){
    selected_factor<-romics_object$metadata[romicsFactorNames(romics_object)==factor,]}
  if(factor=="main"){
    selected_factor<-romics_object$main_factor
    selected_factor<-romics_object$metadata[romicsFactorNames(romics_object)==selected_factor,]
  }else{
    if(factor!="none"){
    stop("The selected <factor> was not present in the list of factor of this romics_object use the function romicsFactorNames() to identify the usable factors.")}}
  #transform in character
  selected_factor<-as.character(t(selected_factor))
  #create a table counting each factor level
  table<-table(selected_factor)
  #verify if none of the levels had 0 or 1 member only and warn if it is the case
  if(sum(table %in% 0:1)>0){
    warning("One or more of the factor levels had only 1 member.")
    warning(table)
    warning("You should consider subsetting the object using the function romicsSubset() to remove the levels with only one columns.")
    stop()
    }
  #establish the levels of this factor
  level_factor<-levels(as.factor(selected_factor))
  #calculate the number of conditions
  length_factor <- length(level_factor)
  #create a count table for each factor
  replicates_factor <- as.double(table)
  names(replicates_factor) <- level_factor
  #calculate the quantity of cells to be full in each given condition
  if(percentage_completeness==0){
    max_empty<-replicates_factor*0
  }else{max_empty <- ceiling((replicates_factor)*(1-percentage_completeness/100))}
  min_full<-replicates_factor-max_empty
  #calculate if the missingness maximum pass for each level of the factor
  list_usable <- data.frame(matrix(nrow=nrow(romics_object$data),ncol=0))
  rownames(list_usable)<-rownames(romics_object$data)
  for (i in 1:length_factor) {
    usable.df <- data.frame()
    usable.df <- romics_object$missingdata[,selected_factor %in% level_factor[i]]
    vec<-rowSums(usable.df!=TRUE)>=min_full[[i]]
    list_usable[,i] <- vec
  }
  colnames(list_usable)<-level_factor
  if (all_groups == "TRUE") {
    usable_groups <- (length_factor - 1) #if TRUE, then apply filter to all groups
  } else {
    if (all_groups == "FALSE") {
      usable_groups <- 0 #if FALSE, filter applies to at least one group at a minimum
    }
}
  usable <- rowSums(list_usable) > usable_groups
  #remove the rows based on this usable vector
  romics_object$data<-romics_object$data[usable,]
  #update the missingness
  romics_object$missingdata<-romics_object$missingdata[usable,]
  #print the number info
  print(paste(sum(usable==FALSE),"rows were removed for the data"))
  print(paste0("Based on the minimum completeness set at ",percentage_completeness,"%"))
  print("at least the following number of sample(s) containing data was required:")
  print(min_full)
  #message with the number of features removed
  print(paste0(nrow(romics_object$data),"/", nrow(romics_object$original_data)," features remained after filtering", " (",round(nrow(romics_object$data)/nrow(romics_object$original_data)*100,2),"%)."))
  romics_object<- romicsUpdateColor(romics_object)
  romics_object<-romicsUpdateSteps(romics_object,arguments)
  #return object
  return(romics_object)
}

#' romicsFilterMissingSamples()
#' @description Filters out the samples of the romics_object$data layer based on their level of missingness below a chosen percentage of completeness for the different features
#' @param romics_object A romics_object created using romicsCreateObject().
#' @param percentage_completeness Numerical value indicating the minimum percentage of completeness required for a sample to be kept.
#' @details This function removes columns (samples) from the romics_object where the percentage of missing data exceeds the chosen threshold. By default, percentage_completeness is set at 50%.
#' @return The function will return a filtered romics_object with samples (columns) removed based on their level of missingness.
#' @author Geremy Clair, Nicholas Day, Modified by AI Assistant
#' @export
romicsFilterMissingSamples <- function(romics_object, percentage_completeness = 50) {
  # Validate the romics_object and parameters
  arguments <- as.list(match.call())
  if (!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }
  if (missing(percentage_completeness)) {
    warning("The acceptable percentage of completeness was not set by the user; it was set at 50% by default")
    percentage_completeness <- 50
  }
  if (percentage_completeness < 0 | percentage_completeness > 100) {
    stop("The completeness has to be between 0 and 100%")
  }
  if (sum(is.na(romics_object$data)) == 0) {
    warning("There are no missing values in this dataset; no samples will be removed.")
    return(romics_object)
  }

  # Calculate the minimum number of non-NA values required per column
  min_full <- ceiling((percentage_completeness / 100) * nrow(romics_object$data))

  # Identify columns to keep based on missing data
  completeness <- colSums(!is.na(romics_object$data)) >= min_full

  # Filter samples (columns) from the data and missingness matrices
  romics_object$data <- romics_object$data[, completeness==TRUE]
  romics_object$missingdata <- romics_object$missingdata[, completeness==TRUE]
  romics_object$metadata <- romics_object$metadata[, completeness==TRUE]

  # Print summary information
  removed_count <- sum(!completeness)
  retained_count <- sum(completeness)
  total_samples <- ncol(romics_object$data) + removed_count # Total original samples

  print(paste(removed_count, "samples were removed from the dataset"))
  print(paste(retained_count, "/", total_samples, "samples remain after filtering",
              "(", round(retained_count / total_samples * 100, 2), "%)."))

  # Update romics_object metadata
  romics_object <- romicsUpdateSteps(romics_object, arguments)

  # Return the filtered object
  return(romics_object)
}

#' romicsPlotMissing()
#' @description Plots the missingness of each sample contained in the romics_object. The colors used for the plotting will correspond to the main_factor of the romics_object.
#' @param romics_object has to be an romics_object created using romicsCreateObject()
#' @param main_factor has to be either "main", "none" or a factor of an romics_object created using romicsCreateObject() the list of factors can be obtained by running the function romicsFactorNames()
#' @details This function does not alter the romics_object, it plots the the missingness of each sample in a barplot.
#' @return This function will return a ggplot2 geom_bar plot. it can then be further visually adjusted using the ggplot2 commands
#' @author Geremy Clair
#' @export
romicsPlotMissing<-function(romics_object,custom_colors= "colorlist"){
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(custom_colors)){custom_colors<-romics_object$custom_colors}
  #load data
  data_NA<-romics_object$missingdata
  #calculate the missing per column
  data_NA<- colSums(data_NA)
  #return a message if there is no missing data
  if(sum(data_NA)==0){stop("Your data has no missing values")}
  #create a color object
  color<-romics_object$metadata[grep("colors_romics",rownames(romics_object$metadata)),]
  color<-as.character(t(color))
  #calculate the percentage of missingness
  percent<- data_NA/nrow(romics_object$missingdata)*100
  #place in a df
  percent_missing<-data.frame(Samples=colnames(romics_object$missingdata),Percent_missing=percent,count=data_NA,color=color)
  percent_missing$Samples<-factor(percent_missing$Samples, levels=unique(percent_missing$Samples))
  #calculate the overall percentage of missingness
  overall_missing<-round(sum(percent_missing$count)/(nrow(romics_object$missingdata)*ncol(romics_object$missingdata))*100,2)
  #create a maximum for the Y_scale
  Y_scale<-percent_missing$Percent_missing
  Y_scale<-round(max(Y_scale),1)+2
  Y_scale<- as.numeric(c(0,Y_scale))
  breaks<-function(x) unique(floor(pretty(seq(0, (max(as.numeric(percent_missing$Percent_missing)) + 1) * 1.1))))
  #plot the result
  p<- ggplot2::ggplot(data=percent_missing,ggplot2::aes(x=Samples,y=Percent_missing))+
    ggplot2::geom_bar(stat="identity",fill=color,alpha=.8)+
    ggplot2::scale_y_continuous(name="Percent_missing", limits=c(0,round(max(percent_missing$Percent_missing),0)+5))+
    ggplot2::ggtitle(paste("Missingness = ", overall_missing,"% of the values"))+
    theme_ROP()
  return(p)
  }

#' romicsPlotPresent()
#' @description Plots the percentage of proteins with valid (non-missing) values for each sample in the romics_object. The colors used for the plotting will correspond to the main_factor of the romics_object.
#' @param romics_object has to be an romics_object created using romicsCreateObject()
#' @param custom_colors has to be either "colorlist" or a character vector containing the colors to be used for the plotting
#' @details This function does not alter the romics_object, it plots the percentage of valid data for each sample in a barplot. This is complementary to romicsPlotMissing() and shows data availability instead of data absence.
#' @return This function will return a ggplot2 geom_bar plot. it can then be further visually adjusted using the ggplot2 commands
#' @author Geremy Clair
#' @export
romicsPlotPresent<-function(romics_object, custom_colors="colorlist"){
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(custom_colors)){custom_colors<-romics_object$custom_colors}

  # Load data
  data_NA <- romics_object$missingdata

  # Calculate the valid (non-missing) per column
  data_valid <- nrow(data_NA) - colSums(data_NA)

  # Return a message if all data is missing (no valid data)
  if(sum(data_valid)==0){stop("Your data has no valid values")}

  # Create a color object
  color <- romics_object$metadata[grep("colors_romics", rownames(romics_object$metadata)), ]
  color <- as.character(t(color))

  # Calculate the percentage of valid data
  percent <- data_valid / nrow(data_NA) * 100

  # Place in a data frame
  percent_present <- data.frame(
    Samples = colnames(data_NA),
    Percent_present = percent,
    count = data_valid,
    color = color
  )
  percent_present$Samples <- factor(percent_present$Samples, levels = unique(percent_present$Samples))

  # Calculate the overall percentage of valid data
  overall_present <- round(sum(percent_present$count) / (nrow(data_NA) * ncol(data_NA)) * 100, 2)

  # Plot the result
  p <- ggplot2::ggplot(data = percent_present, ggplot2::aes(x = Samples, y = Percent_present)) +
    ggplot2::geom_bar(stat = "identity", fill = color, alpha = 0.8) +
    ggplot2::scale_y_continuous(name = "Percent_present", limits = c(0, round(max(percent_present$Percent_present), 0) + 5)) +
    ggplot2::ggtitle(paste("Data Completeness = ", overall_present, "% of the values")) +
    theme_ROP()

  return(p)
}

#' romicsPlotMissingFeatures()
#' @description Plots the missingness of each feature contained in the romics_object.
#' @param romics_object has to be an romics_object created using romicsCreateObject()
#' @param features has to be either empty or a character vector containing a list of features. if a character list of vector is used only the requested features will be shown.
#' @param reorder has to be either empty (original order will be used), 'ascending', or 'descending' to indicate in what order the features should be displayed on the plot.
#' @details This function does not alter the romics_object, it plots the the missingness of each feature for all the samples in a barplot.
#' @return This function will return a ggplot2 geom_bar plot. it can then be further visually adjusted using the ggplot2 commands.
#' @author Geremy Clair
#' @export
romicsPlotMissingFeatures<-function(romics_object=romics_object,features=c(feature_a,feature_b),reorder=c(NULL,"descending","ascending")){
  if(!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")}
  if(missing(features)){features<-NULL}
  if(missing(reorder)){reorder<-NULL}
  m<-data.frame(Features=rownames(romics_object$missingdata),Percentage_missing=rowSums(romics_object$missingdata)/ncol(romics_object$missingdata)*100)

  if(!is.null(features)){
    if(sum(features %in% m)==length(features)){
      stop("Not all 'features' are in the feature list of the romics_object")}
    m<-m[m[,1] %in% features,]}
  if(!is.null(reorder)){
    if(reorder=="ascending"){m<-m[order(m$Percentage_missing),]}
    if(reorder=="descending"){m<-m[order(m$Percentage_missing,decreasing = T),]}
  }else{warning("<reorder> was either empty or neither 'ascending', or 'descending', the data was not reordered.")}
  m$Features<-factor(m$Features,levels=m$Features)

  ggplot(m, aes(x=Features,y=Percentage_missing))+ geom_bar(stat="identity")+theme_ROP()
}

#' romicsMissingFeatures()
#' @description returns a table containing the missingness of each feature contained in the romics_object.
#' @param romics_object has to be an romics_object created using romicsCreateObject()
#' @param reorder has to be either empty (original order will be used), 'ascending', or 'descending' to indicate in what order the features should be displayed on the plot.
#' @details This function does not alter the romics_object, it generates a table containing the percentage of missingness for each feature of the romics_object
#' @return This function will return a ggplot2 geom_bar plot. it can then be further visually adjusted using the ggplot2 commands.
#' @author Geremy Clair
#' @export
romicsMissingFeatures<-function(romics_object=romics_object,reorder=c(NULL,"descending","ascending")){
  if(!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")}
  if(missing(reorder)){reorder<-NULL}
  m<-data.frame(Features=rownames(romics_object$missingdata),Percentage_missing=rowSums(romics_object$missingdata)/ncol(romics_object$missingdata)*100)
  if(!is.null(reorder)){
    if(reorder=="ascending"){m<-m[order(m$Percentage_missing),]}
    if(reorder=="descending"){m<-m[order(m$Percentage_missing,decreasing = T),]}
  }
  rownames(m)<-1:nrow(m)
 return(m)
}

#' romicsFreqMissingFeatures()
#' @description Plots the frequency of the feature missingness for all features contained in an romics_object.
#' @param romics_object has to be an romics_object created using romicsCreateObject()
#' @param binwidth has to be a numeric value indicating the width of the desired frequency bins, 1 will be used by default.
#' @details This function does not alter the romics_object, it generates a table containing the percentage of missingness for each feature of the romics_object.
#' @return This function will return a ggplot2 geom_bar plot. it can then be further visually adjusted using the ggplot2 commands.
#' @author Geremy Clair
#' @export
romicsFreqMissingFeatures<-function(romics_object=romics_object,binwidth=1){
  if(!is.romicsObject(romics_object)){stop("The object was not created with the function 'romicsCreateObject()'.")}
  if(missing(binwidth)){binwidth=1}
  m<-data.frame(Features=rownames(romics_object$missingdata),Percentage_missing=rowSums(romics_object$missingdata)/ncol(romics_object$missingdata)*100)
  ggplot(m,aes(x=Percentage_missing))+geom_bar(stat="bin",binwidth=binwidth)+theme_ROP()+ylab(label = "Frequency")
}

#' romicsDeleteRowMissing()
#' @description Removes the rows that have any missing values from the romics_object
#' @param romics_object has to be an romics_object created using romicsCreateObject()
#' @details  This function will remove rows that have any missing values in the data and missingdata layers
#' @return  The function will return a filtered romics_object with the rows of the data and missing data object removed when appropriate.
#' @author Geremy Clair
#' @export
romicsDeleteRowMissing<-function(romics_object){
  arguments<-as.list(match.call())
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  #load data and missingness
  data<- romics_object$data
  missing<- romics_object$missingdata
  #remove those from the data and missing
  romics_object$data<-data[rowSums(missing)==0,]
  romics_object$missingdata<-missing[rowSums(missing)==0,]
  #update the colors
  romics_object<- romicsUpdateColor(romics_object)
  #indicate how many rows were removed
  print(paste(sum(rowSums(missing)!=0),"row(s) removed for the data because it contained missing values"))
  romics_object<-romicsUpdateSteps(romics_object,arguments)
  return(romics_object)
  }

#' imputeMissingEval()
#' @description Plots the distribution of the data to be imputed using the imputeMissing() function. Enables  to optimize the parameters prior to apply the imputeMissing() function.
#' @param romics_object has to be an romics_object created using romicsCreateObject()
#' @param nb_stdev number of standard deviations away from the median to place the missingness. Default: 1.8
#' @param width_stdev width of the random values used for imputation expressed in standard deviation from the data. Default: 0.5
#' @param bin binwidth for the histogram. Default: 1
#' @param scale_x numeric vector of length 2 for x-axis limits. Default: c(-10,10)
#' @param sample_size maximum number of data points to use for visualization (for performance with large datasets). Default: 100000
#' @param set.seed logical indicating whether to set a seed for reproducible results. Default: TRUE
#' @param seed numeric value for the random seed when set.seed=TRUE. Default: 42
#' @details  This function does not alter the romics_object, it plots the distribution of the whole data and of the imputed values using the method described in the Perseus paper by Tyranova et al. 2016 (Nature Method). By default the data is imputed with values in a normal distribution 1.8 standard deviation away from the median  and a width of distribution of 0.5.
#' @return This function will return a ggplot2 geom_histogram plot. it can then be further visually adjusted using the ggplot2 commands
#' @author Geremy Clair
#' @export
imputeMissingEval <- function(romics_object, nb_stdev = 1.8, width_stdev = 0.5, bin = 1,
                              scale_x = "auto", sample_size = 100000, set.seed = TRUE, seed = 42) {

  if (!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }

  # Set seed for reproducibility if requested
  if (set.seed) {
    set.seed(seed)
  }

  # Check if there are missing values (FIXED)
  if (sum(romics_object$missingdata == TRUE) == 0) {
    warning("No values were missing according to your missingdata layer")
  }

  # Extract data efficiently
  data <- romics_object$data
  missing_mask <- romics_object$missingdata

  # Get non-missing values more efficiently
  non_missing_values <- data[!missing_mask & !is.na(data)]

  # Sample data if too large for performance
  if (length(non_missing_values) > sample_size) {
    message(paste("Sampling", sample_size, "values from", length(non_missing_values), "non-missing values for visualization"))
    non_missing_values <- sample(non_missing_values, sample_size)
  }

  # Calculate statistics efficiently
  data_median <- median(non_missing_values, na.rm = TRUE)
  data_sd <- sd(non_missing_values, na.rm = TRUE)

  # Generate imputed values
  n_missing <- sum(missing_mask == TRUE)
  imputation <- rnorm(n_missing,
                      mean = data_median - nb_stdev * data_sd,
                      sd = width_stdev * data_sd)

  # Sample imputed values if too many
  if (length(imputation) > sample_size) {
    message(paste("Sampling", sample_size, "values from", length(imputation), "imputed values for visualization"))
    imputation <- sample(imputation, sample_size)
  }

  # Create combined dataset efficiently
  combined_data <- data.frame(
    combined = c(non_missing_values, imputation),
    data_type = factor(c(rep("Data", length(non_missing_values)),
                         rep("Imputed values", length(imputation))),
                       levels = c("Data", "Imputed values")),
    stringsAsFactors = FALSE
  )

  # Create plot
  p <- ggplot2::ggplot(combined_data, ggplot2::aes(x = combined, fill = data_type)) +
    ggplot2::geom_histogram(position = "identity", alpha = 0.8, binwidth = bin) +
    ggplot2::xlab("data distribution") +
    theme_ROP() +
    ggplot2::scale_fill_manual(values = c("Data" = "gray", "Imputed values" = "goldenrod1")) +
    ggplot2::theme(legend.position = "right")

  # Handle x-axis limits intelligently
  if (is.character(scale_x) && scale_x == "auto") {
    # Automatically determine limits based on data
    x_min <- min(combined_data$combined, na.rm = TRUE)
    x_max <- max(combined_data$combined, na.rm = TRUE)

    # Add some padding (5% on each side)
    x_range <- x_max - x_min
    x_padding <- x_range * 0.05

    p <- p + ggplot2::scale_x_continuous(limits = c(x_min - x_padding, x_max + x_padding))

  } else if (is.numeric(scale_x) && length(scale_x) == 2) {
    # Use custom limits
    p <- p + ggplot2::scale_x_continuous(limits = scale_x)

  } else if (!is.null(scale_x) && scale_x != "auto") {
    warning("scale_x should be 'auto' or a numeric vector of length 2. Using automatic limits.")
    # Fall back to auto
    x_min <- min(combined_data$combined, na.rm = TRUE)
    x_max <- max(combined_data$combined, na.rm = TRUE)
    x_range <- x_max - x_min
    x_padding <- x_range * 0.05
    p <- p + ggplot2::scale_x_continuous(limits = c(x_min - x_padding, x_max + x_padding))
  }
  # If scale_x is NULL, ggplot will use default automatic scaling

  return(p)
}

#' imputeMissing()
#' @description Imputes the data using a normal distribution down-shifted from the median by a user defined number of standard deviations and a user defined width. the distribution of the imputed data can be evaluated using the function imputeMissingEval().
#' @param romics_object has to be an romics_object created using romicsCreateObject()
#' @param nb_stdev Numeric. Number of standard deviations to shift the imputation distribution from the median. Default: 1.8
#' @param width_stdev Numeric. Width of the imputation distribution in standard deviations. Default: 0.5
#' @param seed Integer. Random seed for reproducible results. Default: 42
#' @details This function will impute the data using the method described in the Perseus paper by Tyranova et al. 2016 (Nature Method). By default the data is imputed with values in a normal distribution 1.8 standard deviation away from the median and a width of distribution of 0.5.
#' @return The function will return a modified romics_object that will have imputed data, however the missingdata layer will conserve the location of the missingness, the missingness can subsequently be restored using the function romicsRestoreMissing().
#' @author Geremy Clair
#' @export
imputeMissing <- function(romics_object, nb_stdev=1.8, width_stdev=0.5, seed=42) {
  arguments <- as.list(match.call())

  # Input validation
  if(!is.romicsObject(romics_object) || missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }
  if(missing(nb_stdev)) {nb_stdev <- 1.8}
  if(missing(width_stdev)) {width_stdev <- 0.5}
  if(missing(seed)) {seed <- 42}

  # Check if there are missing values to impute
  if(sum(romics_object$missingdata) == 0) {
    warning("No values were imputed as the missingdata layer was indicating no existing missing values")
    return(romics_object)
  }

  # Extract data and preserve structure
  data <- romics_object$data
  missingdata <- romics_object$missingdata
  original_rownames <- rownames(data)
  original_colnames <- colnames(data)
  was_dataframe <- is.data.frame(data)

  # Convert to matrices and ensure proper dimensions
  if(was_dataframe) {
    data <- data.matrix(data)
  }
  if(!is.matrix(data)) {
    data <- as.matrix(data)
  }
  if(!is.matrix(missingdata)) {
    missingdata <- as.matrix(missingdata)
  }

  # Ensure proper types
  storage.mode(data) <- "double"
  missingdata <- missingdata == TRUE  # Ensure logical

  # Verify dimensions match
  if(!identical(dim(data), dim(missingdata))) {
    stop("Data and missingdata dimensions do not match")
  }

  # Restore names
  rownames(data) <- original_rownames
  colnames(data) <- original_colnames

  # Calculate statistics on non-missing values only
  if(requireNamespace("matrixStats", quietly = TRUE)) {
    col_medians <- matrixStats::colMedians(data, na.rm = TRUE)
  } else {
    col_medians <- apply(data, 2, median, na.rm = TRUE)
  }
  global_median <- median(col_medians, na.rm = TRUE)

  col_sds <- apply(data, 2, sd, na.rm = TRUE)
  global_sd <- median(col_sds, na.rm = TRUE)

  # Validate statistics
  if(is.na(global_median) || is.na(global_sd) || global_sd <= 0) {
    stop("Unable to calculate valid statistics for imputation")
  }

  # Set seed for reproducibility
  set.seed(seed)

  # Calculate imputation parameters
  imputation_mean <- global_median - nb_stdev * global_sd
  imputation_sd <- width_stdev * global_sd

  # Count missing values and perform imputation
  missing_count <- sum(missingdata)

  if(missing_count > 0) {
    message(paste("Imputing", missing_count, "missing values..."))

    # Generate random values
    random_values <- rnorm(missing_count, mean = imputation_mean, sd = imputation_sd)

    # Find positions to impute (use the missingdata matrix)
    missing_positions <- which(missingdata, arr.ind = TRUE)

    # Impute values one by one to maintain order (fast enough for most cases)
    for(i in seq_len(nrow(missing_positions))) {
      row_idx <- missing_positions[i, 1]
      col_idx <- missing_positions[i, 2]
      data[row_idx, col_idx] <- random_values[i]
    }

    # Verify imputation worked
    remaining_na <- sum(is.na(data))
    if(remaining_na > 0) {
      warning(paste("Warning:", remaining_na, "NA values remain after imputation"))
    }

    message("Imputation completed successfully")
  }

  # Convert back to original format - FAST VERSION
  if(was_dataframe) {
    # Fast conversion - since matrix is already numeric, this preserves it
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    # No need for the slow loop - the data.frame will inherit the numeric type from the matrix
  }

  # Update romics object
  romics_object$data <- data
  romics_object <- romicsUpdateSteps(romics_object, arguments)

  return(romics_object)
}

#' imputeMinDiv2()
#' @description Imputes the data layer of the romics_object using the minimum value of the table divided by 2. this function will work if the data only contains positive values.
#' @param romics_object has to be an romics_object created using romicsCreateObject()
#' @details This function will use the minimum value of the data table divided by 2 to impute the missing values of the data layer
#' @return  The function will return a modified romics_object that will have imputed data, however the missingdata layer will conserve the location of the missingness, the missingness can subsequently be restored using the function romicsRestoreMissing().
#' @author Geremy Clair
#' @export
imputeMinDiv2<- function(romics_object){
  arguments<-as.list(match.call())
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(min(romics_object$data,na.rm =T)<0){stop("To apply this function all the values of the romics_object data layer have to be either positive or NAs")}
  if(sum(romics_object$missingdata==TRUE)==0){warning("No values were imputed as the missingdata layer was indicating no existing missing values")}
  min_div2<-min(romics_object$data, na.rm=TRUE)/2
  romics_object$data[romics_object$missingdata==TRUE]<-min_div2
  romics_object<-romicsUpdateSteps(romics_object,arguments)
  return(romics_object)
  }

#' romicsRestoreMissing()
#' @description Restores the missingness of the data layer of the romics_object based on the content of the missingdata layer.
#' @param romics_object has to be an romics_object created using romicsCreateObject()
#' @details  This function will renmove any imputed value based on the content of the missingdata layer.
#' @return  The function will return a modified romics_object that will have NA instead of the imputed data.
#' @author Geremy Clair
#' @export
romicsRestoreMissing<-function(romics_object){
  arguments<-as.list(match.call())
  if(!is.romicsObject(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  romics_object$data[romics_object$missingdata==TRUE] <- NA
  romics_object<-romicsUpdateSteps(romics_object,arguments)
  return(romics_object)
}

#' romicsValidFeatures()
#' @description Calculates and returns the number of features with valid (non-missing) values for each sample in the romics_object.
#' @param romics_object has to be a romics_object created using romicsCreateObject()
#' @param return_type Character string specifying the output format. Options are: "table" (default) returns a data frame with sample names and counts; "vector" returns a named numeric vector; "plot" returns a ggplot2 bar plot.
#' @param custom_colors Character string or vector. If "colorlist" (default), uses the colors from romics_object. If a vector of colors is provided, it will be used for the plot (only applicable when return_type = "plot").
#' @details This function counts the number of non-missing values (valid features) for each sample (column) in the romics_object$data layer. The results can be returned as a table, vector, or plot depending on the return_type parameter.
#' @return Depending on return_type:
#' @author Geremy Clair
#' @export
romicsValidFeatures <- function(romics_object, return_type = "table", custom_colors = "colorlist") {
  # Input validation
  if(!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }

  # Set defaults
  if(missing(return_type)) {return_type <- "table"}
  if(missing(custom_colors)) {custom_colors <- "colorlist"}

  # Validate return_type
  if(!return_type %in% c("table", "vector", "plot")) {
    stop("return_type must be one of: 'table', 'vector', or 'plot'")
  }

  # Calculate valid features per sample (column)
  valid_counts <- colSums(!romics_object$missingdata)

  # Calculate total features
  total_features <- nrow(romics_object$data)

  # Calculate percentages
  valid_percent <- (valid_counts / total_features) * 100

  # Return based on type
  if(return_type == "vector") {
    names(valid_counts) <- colnames(romics_object$data)
    return(valid_counts)
  }

  if(return_type == "table") {
    result_table <- data.frame(
      Sample = colnames(romics_object$data),
      Valid_features = valid_counts,
      Percent_valid = round(valid_percent, 2),
      stringsAsFactors = FALSE
    )
    rownames(result_table) <- NULL
    return(result_table)
  }

  if(return_type == "plot") {
    # Get colors
    if(length(custom_colors) == 1 && custom_colors == "colorlist") {
      colors <- as.character(t(romics_object$metadata[grep("colors_romics", rownames(romics_object$metadata)), ]))
    } else {
      colors <- custom_colors
    }

    # Create data frame for plotting
    plot_data <- data.frame(
      Sample = factor(colnames(romics_object$data), levels = colnames(romics_object$data)),
      Valid_features = valid_counts,
      Percent_valid = valid_percent,
      Color = colors,
      stringsAsFactors = FALSE
    )

    # Calculate overall percentage of valid values
    overall_valid <- round(sum(valid_counts) / (total_features * ncol(romics_object$data)) * 100, 2)

    # Create the plot
    p <- ggplot(data = plot_data, aes(x = Sample, y = Valid_features)) +
      geom_bar(stat = "identity", fill = colors, alpha = 0.8) +
      scale_y_continuous(name = "Number of valid features",
                         limits = c(0, total_features * 1.05)) +
      ggtitle(paste("Valid values =", overall_valid, "% of the data")) +
      theme_ROP() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    return(p)
  }
}

#' romicsValidFeaturesPercent()
#' @description Calculates and displays the percentage of valid (non-missing) samples for each feature in a romics_object.
#' @param romics_object A romics_object created with createRomicsObject()
#' @param return_type Character string specifying the format for the output: "table" (default), "vector", or "plot"
#'   - "table": Returns a data frame with Feature, Valid_samples, and Percent_valid columns
#'   - "vector": Returns a named numeric vector with feature names and their completeness percentages
#'   - "plot": Returns a ggplot2 bar plot visualization of feature completeness
#' @param reorder Character string or NULL specifying whether to sort by completeness percentage: NULL (original order, default), "ascending", or "descending"
#' @param features Character vector of specific feature names to include in the analysis. If NULL, all features are analyzed.
#' @details This function calculates the percentage of samples with valid (non-missing) values for each feature using the missingdata layer. Results can be returned as a summary table, a named vector, or a visualization. Useful for identifying features with varying levels of data completeness.
#' @return Depends on return_type: data frame (table), named numeric vector (vector), or ggplot2 object (plot)
#' @author Geremy Clair
#' @export
romicsValidFeaturesPercent <- function(romics_object, return_type = "table", reorder = NULL, features = NULL) {
  # Input validation
  if(!is.romicsObject(romics_object) | missing(romics_object)) {
    stop("romics_object is missing or is not in the appropriate format")
  }

  # Set defaults
  if(missing(return_type)) {return_type <- "table"}
  if(missing(reorder)) {reorder <- NULL}
  if(missing(features)) {features <- NULL}

  # Validate return_type
  if(!return_type %in% c("table", "vector", "plot")) {
    stop("return_type must be one of: 'table', 'vector', or 'plot'")
  }

  # Validate reorder
  if(!is.null(reorder) && !reorder %in% c("ascending", "descending")) {
    warning("reorder was neither NULL, 'ascending', nor 'descending'. Using original order.")
    reorder <- NULL
  }

  # Calculate valid samples per feature (row)
  valid_counts <- rowSums(!romics_object$missingdata)
  total_samples <- ncol(romics_object$data)
  valid_percent <- (valid_counts / total_samples) * 100

  # Create data frame
  result_df <- data.frame(
    Feature = rownames(romics_object$data),
    Valid_samples = valid_counts,
    Percent_valid = round(valid_percent, 2),
    stringsAsFactors = FALSE
  )

  # Filter for specific features if requested
  if(!is.null(features)) {
    if(!all(features %in% result_df$Feature)) {
      missing_features <- features[!features %in% result_df$Feature]
      stop(paste("The following features are not in the romics_object:",
                 paste(missing_features, collapse = ", ")))
    }
    result_df <- result_df[result_df$Feature %in% features, ]
  }

  # Reorder if requested
  if(!is.null(reorder)) {
    if(reorder == "ascending") {
      result_df <- result_df[order(result_df$Percent_valid), ]
    } else if(reorder == "descending") {
      result_df <- result_df[order(result_df$Percent_valid, decreasing = TRUE), ]
    }
  }

  # Return based on type
  if(return_type == "vector") {
    result_vector <- result_df$Percent_valid
    names(result_vector) <- result_df$Feature
    return(result_vector)
  }

  if(return_type == "table") {
    rownames(result_df) <- NULL
    return(result_df)
  }

  if(return_type == "plot") {
    # Make feature a factor with levels in the current order
    result_df$Feature <- factor(result_df$Feature, levels = result_df$Feature)

    # Create the plot
    p <- ggplot(result_df, aes(x = Feature, y = Percent_valid)) +
      geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
      scale_y_continuous(name = "Percent of valid samples (%)",
                         limits = c(0, 105)) +
      xlab("Features") +
      ggtitle("Data completeness per feature") +
      theme_ROP() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

    return(p)
  }
}
