#' romicsZeroToMissing()
#' @description Replaces zeros in the romics_object by NA values
#' @param romics_object has to be an romics_object created using romicsCreateObject()
#' @details This function will convert 0 values to NA in the data and missingdata layers
#' @return This function returns the transformed romics_object with updated data and missingdata layers
#' @author Geremy Clair
#' @export
romicsZeroToMissing<-function(romics_object){
  arguments<-as.list(match.call())
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}

  romics_object$data[romics_object$data==0]<-NA
  romics_object$missingdata<-data.frame(is.na(romics_object$data))

  romics_object<-romicsUpdateSteps(romics_object,arguments)

  return(romics_object)
}

#' romicsFilterMissing()
#' @description Filters out the variables of the romics_object below the choosen percentage of completeness. The percentage of completeness can either be global or by factor of a given factor (in this later case if the percentage of completeness is achieved for at least one level of the factor it the variable will be kept).
#' @param romics_object A romics_object created using romicsCreateObject()
#' @param percentage_completeness Numerical vector indicating the minimum percentage of data to be considered
#' @param main_factor has to be either "main", "none" or a factor of an romics_object created using romicsCreateObject() the list of factors can be obtained by running the function romicsFactorNames() on the romics_object
#' @param all_groups if this parameter is TRUE the completeness requirement is for each and every group, if not, the completeness requirement is for at least one group.
#' @details  This function will use the completeness of the protein in the overall samples (when none is used as factor), or of a given level of a specific defined factor (in this case the factor has to be set to either "main" or to the given factor of filtering). By default main_factor is the main factor of the object, the percentage_completeness is set at 50%
#' @return  The function will return a filtered romics_object with the rows of the data and missing data object removed when appropriate.
#' @author Geremy Clair, Nicholas Day
#' @export
romicsFilterMissing<-function(romics_object, percentage_completeness=50, main_factor = "main",all_groups=FALSE){
  arguments<-as.list(match.call())
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(main_factor)){
    warning("your main_factor was missing the main_factor of the romics_object was used")
    main_factor<-"main"}
  if(missing(percentage_completeness)){
    warning("The acceptable percentage of completeness was not set by the user it was set at 50% by default")
    percentage_completeness<-50}
  if(percentage_completeness<0 & percentage_completeness>100){stop("the completeness has to be comprised between 0 and 100 %")}
  if(sum(is.na(romics_object$data))==0){warning("There is no missing values in this dataset data to be removed")}
  if(missing(all_groups)){all_groups=FALSE}
  #extract main factor
  if(main_factor=="none"){
    selected_factor<-rep("overal_sample_number",ncol(romics_object$metadata))}
  if(main_factor %in% romicsFactorNames(romics_object) ){
    selected_factor<-romics_object$metadata[romicsFactorNames(romics_object)==main_factor,]}
  if(main_factor=="main"){
    selected_factor<-romics_object$main_factor
    selected_factor<-romics_object$metadata[romicsFactorNames(romics_object)==selected_factor,]
  }else{
    if(main_factor!="none"){
    stop("The selected <main_factor> was not present in the list of factor of this romics_object use the function romicsFactorNames() to identify the usable factors.")}}

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
  }else{max_empty <- floor((replicates_factor)*(1-percentage_completeness/100))}
  min_full<-replicates_factor-max_empty

  #calculate if the missingness maximum pass for each level of the factor
  list_usable <- data.frame(matrix(nrow=nrow(romics_object$data),ncol=0))
  rownames(list_usable)<-rownames(romics_object$data)
  for (i in 1:length_factor) {
    usable.df <- data.frame()
    usable.df <- romics_object$data[,grepl(level_factor[i],selected_factor)]
    vec<-rowSums(!is.na(usable.df))>=min_full[[i]]
    list_usable[,i] <- vec
  }
  colnames(list_usable)<-level_factor

  if (all_groups == "TRUE") {
    usable_groups <- (length_factor - 1) #if TRUE, then apply filter to all groups
  }
  else {
    if (all_groups == "FALSE") {
      usable_groups <- 0 #if FALSE, filter applies to at least one group at a minimum
    }
  }
  usable <- rowSums(list_usable) > usable_groups

  #remove the rows based on this usable vector
  romics_object$data<-romics_object$data[usable,]

  #update the missingness
  romics_object$missingdata<-data.frame(is.na(romics_object$data))

  #print the number info
  print(paste(sum(usable==FALSE),"rows were removed for the data"))
  print(paste0("Based on the minimum completeness set at ",percentage_completeness,"%"))
  print("at least the following number of sample(s) containing data was required:")
  print(min_full)

  romics_object<- romicsUpdateColor(romics_object)
  romics_object<-romicsUpdateSteps(romics_object,arguments)

  #return object
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
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
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
  p<- ggplot(data=percent_missing,aes(x=Samples,y=Percent_missing))+
    geom_bar(stat="identity",fill=color,alpha=.8)+
    scale_y_continuous(name="Percent_missing", limits=c(0,round(max(percent_missing$Percent_missing),0)+5))+
    ggtitle(paste("Missingness = ", overall_missing,"% of the values"))+
    theme_ROP()

  return(p)
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
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}

  #load data and missingness
  data<- romics_object$data
  missing<- romics_object$missingdata

  #establish the rows to remove
  uncomplete_row<-rowSums(missing)>0

  #remove those from the data and missing
  romics_object$data<-data[!uncomplete_row,]
  romics_object$missingdata<-missing[!uncomplete_row,]

  #update the colors
  romics_object<- romicsUpdateColor(romics_object)

  #indicate how many rows were removed
  print(paste(sum(uncomplete_row),"row(s) removed for the data because it contained missing values"))

  romics_object<-romicsUpdateSteps(romics_object,arguments)

  return(romics_object)
  }

#' imputeMissingEval()
#' @description Plots the distribution of the data to be imputed using the imputeMissing() function. Enables  to optimize the parameters prior to apply the imputeMissing() function.
#' @param romics_object has to be an romics_object created using romicsCreateObject()
#' @details  This function does not alter the romics_object, it plots the distribution of the whole data and of the imputed values using the method described in the Perseus paper by Tyranova et al. 2016 (Nature Method). By default the data is imputed with values in a normal distribution 1.8 standard deviation away from the median  and a width of distribution of 0.5.
#' @return This function will return a ggplot2 geom_bar plot. it can then be further visually adjusted using the ggplot2 commands
#' @author Geremy Clair
#' @export
imputeMissingEval<-function(romics_object,nb_stdev=1.8,width_stdev=0.5,bin=1,scale_x=c(-10,10)){
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(nb_stdev)){nb_stdev<-1.8}
  if(missing(width_stdev)){width_stdev<-0.5}
  if(missing(bin)){bin=1}
  if(sum(romics_object$missingdata==TRUE)==0){warning("No values were missing in according to your missingdata layer")}

  #extract data
  data<- romics_object$data
  #melt all the columns
  melted_data<-melt(data[romics_object$missingdata==F])

  #remove missing
  melted_data<-melted_data[!is.na(melted_data)]

  #generate a randomly distributed
  imputation<- rnorm(sum(romics_object$missingdata==T),(median(melted_data,na.rm=TRUE)-nb_stdev*sd(melted_data)), width_stdev*sd(melted_data))

  #combined
  data<-c(melted_data,imputation)

  #data type attribution
  data_type<- c(rep("Data",length(melted_data)) ,rep("Imputed values",length(imputation)))

  #colors
  color<- c(rep("gray",length(melted_data)) ,rep("goldenrod1",length(imputation)))

  #combine in one object
  melted_data<-data.frame(cbind(data,data_type,color))
  colnames(melted_data)<-c("combined","data_type","color")


  p<-ggplot(melted_data, aes(x = data , fill=data_type)) +
    geom_histogram(position="identity", alpha=0.8,binwidth = bin)+
    xlab("data distribution")+
    theme_ROP()+
    scale_fill_manual(values=c("gray","goldenrod1"))+
    theme(legend.position="right")

  if (is.numeric(scale_x)&&length(scale_x)==2){plot<-plot+scale_x_continuous(limits = scale_x)}

  return(p)

}

#' imputeMissing()
#' @description Imputes the data using a normal distribution down-shifted from the median by a user defined number of standard deviations and a user defined width. the distribution of the imputed data can be evaluated using the function imputeMissingEval().
#' @param romics_object has to be an romics_object created using romicsCreateObject()
#' @details  This function will impute the data using the method described in the Perseus paper by Tyranova et al. 2016 (Nature Method). By default the data is imputed with values in a normal distribution 1.8 standard deviation away from the median  and a width of distribution of 0.5.
#' @return  The function will return a modified romics_object that will have imputed data, however the missingdata layer will conserve the location of the missingness, the missingness can subsequently be restored using the function romicsRestoreMissing().
#' @author Geremy Clair
#' @export
imputeMissing<-function(romics_object,nb_stdev=1.8,width_stdev=0.5, seed=42){
  arguments<-as.list(match.call())
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(nb_stdev)){nb_stdev<-1.8}
  if(missing(width_stdev)){width_stdev<-0.5}
  if(missing(seed)){seed<-42}
  if(sum(romics_object$missingdata==TRUE)==0){warning("No values were imputed as the missingdata layer was indicating no existing missing values")}

  #extract data
  data<- romics_object$data
  missingdata<-romics_object$missingdata


  #calculate the median
  median <-numeric()
  for (i in 1:ncol(data)){
    median[i]<-median(data[,i],na.rm=TRUE)
  }
  median<-median(median)

  #calculate the sd
  sd<-numeric()
  for (i in 1:ncol(data)){
    sd[i]<-sd(data[,i],na.rm=TRUE)
  }
  sd<-median(sd)

  set.seed(seed)
  for (i in 1:ncol(data)){
    for(j in 1:nrow(data)){
    if(missingdata[j,i]==T){data[j,i]<-rnorm(1,mean= median-nb_stdev*sd, sd=width_stdev*sd)}
  }}

  #update the data
  romics_object$data<-data

  #update the steps
  romics_object<-romicsUpdateSteps(romics_object,arguments)

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
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(min(romics_object$data,na.rm =T)<0){stop("To apply this function all the values of the romics_object data layer have to be either positive or NAs")}
  if(sum(romics_object$missingdata==TRUE)==0){warning("No values were imputed as the missingdata layer was indicating no existing missing values")}

  #calculate the min / 2 of the table (disregarding missing data)
  min_div2<-min(romics_object$data, na.rm=TRUE)/2

  #replace missing with this value
  romics_object$data[romics_object$missingdata==TRUE]<-min_div2

  #update the steps
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
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}

  #replace the originally missing data by NAs
  romics_object$data[romics_object$missingdata==TRUE] <- NA

  #update the steps
  romics_object<-romicsUpdateSteps(romics_object,arguments)

  return(romics_object)
}
