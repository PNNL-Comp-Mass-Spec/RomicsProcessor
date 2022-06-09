#' romicsCreateObject()
#' @description Create a romics_object by combining a data and a metadata data frames
#' @param data A data frame corresponding to the data.
#' @param metadata A data frame corresponding to the metadata, the columns must be the same as the data ones.
#' @param main_factor The rowname of a metadata factor, by default the first row will be used as main factor.
#' @param custom_colors A character vector containing the colors you want to use for your figures generated from an romics_object.
#' @param omics_type A character vector of length 1 indicating the type of omics data used.
#' @param quantif_type A character vector of length 1 indicating the type of quantification performed.
#' @return This function generate an romics_object containing the following layers : data, metadata, missingdata, original_data, main_factor, colors, steps, custom_colors, omics_type, quantif_type
#' @examples ROP_romics_object <- romicsCreateObject(ROP_data,ROP_metadata,main_factor="condition")
#' @author Geremy Clair
#' @export
romicsCreateObject<-function(data, metadata, main_factor="none", custom_colors=ROP_colors, omics_type="unknown", quantif_type="unknown"){
  arguments<-as.list(match.call())
  funName<-arguments[[1]]
  if(missing(data)){stop("<data> is missing")}
  if(missing(metadata)) {stop("<metadata> is missing")}
  if(sum(is.na(metadata))>0){stop("<metadata> should not contain NAs.")}
  if(missing(custom_colors)) {custom_colors<-ROP_colors}
  if(missing(omics_type)) {omics_type<-"unknown"}
  if(missing(quantif_type)){quantif_type<-"unknown"}
  if(missing(main_factor)){main_factor<-"none"}
  if(typeof(data)!="list"){data<-data.frame(data)}
  if(typeof(metadata)!="list"){data<-data.frame(metadata)}

  argumentsNames<-names(arguments)[2:length(arguments)]

#check if the IDs of the two objects are the same
  colnames(data)[1]<-"ID"
  colnames(metadata)[1]<-"ID"
  if(!sum(colnames(data) %in% colnames(metadata))==ncol(metadata)){stop("data and metadata columns names are different or don't have the same number of columns")}

#order the columns of metadata to match data
  metadata<- metadata[, match(colnames(data), colnames(metadata))]

# from first column make rownames for data and metadata
  rownames(data)<-data[,1]
  rownames(metadata)<-metadata[,1]

# remove first column
  data<-data[,2:ncol(data)]
  metadata<-metadata[,2:ncol(metadata)]

#ensure that the data is numeric
  for(i in 1:ncol(data)){
    data[,i]<-as.double(data[,i])
  }

# original_data : save the state of the data at its loading this will allow for re-running the analysis form the original data
# nothing else than the ID column setting was done on this data
  original_data<-data

#order metadata columns based on data column order
  metadata<-metadata[,order(match(colnames(metadata),colnames(data)))]

#define missing table based on NAs in the data object
  missingdata<- data
  missingdata<-data.frame(is.na(data))

# select the main_factor to use based in the user input this factor will be used as default if not indicated by user in a variety of functions
  if(main_factor=="none"){
    warning("your main_factor was missing the first row of your metadata was used as factor")
    main_factor<-metadata[1,]}else{
      if(sum(rownames(metadata)==main_factor)==1){ main_factor<-metadata[rownames(metadata)==main_factor,]}else {
        stop("Your main_factor is not present in your metadata")
      }
      }

  factor<-rownames(main_factor)

  #data, metadata,and missingdata based on selected factor, if none is selected, use the first row as main factor

  # find level order in main_factor
  main_factor<- as.character(as.factor(t(main_factor)[,1]))

  names(main_factor)<- colnames(metadata)
  main_factor_lvl<- data.frame(unique(main_factor))

  tp_data<-data.frame(matrix(nrow=nrow(data),ncol=0))
  tp_metadata<-data.frame(matrix(nrow=nrow(metadata),ncol=0))
  tp_missingdata<-data.frame(matrix(nrow=nrow(missingdata),ncol=0))
  tp_main_factor<-character()

  #reorder data, metadata and missingdata based on the factors which will make visualization more pretty
  for(i in 1:nrow(main_factor_lvl))
  {
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

  #add a color_romics line in metadata
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

 #create a vector containing the colors for the whole data points and store it in colors layer
  colors_romics<- t(colors_romics)
  fill <- character(length = 0)
  for (i in 1:length(colors_romics))
    {
    fill<- c(fill,rep(as.character(colors_romics[i]),nrow(data)))
    }
  colors<-fill
  remove(fill,colors_romics, main_factor_lvl,i)

  #create a steps character vector with the date time of the creation of the object
  steps<- c("romics_object",paste0("date|",gsub(" ","_",format(Sys.time(),"%b_%d_%Y_%X")),"|romicsCreateObject"))
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

  #create the final list
  l<-list(data, metadata, missingdata, original_data, factor, colors, steps, dependencies, custom_colors, omics_type, quantif_type)
  names(l)<-c("data","metadata","missingdata","original_data", "main_factor","colors","steps","dependencies","custom_colors","omics_type","quantif_type")
  class(l)<-"romics_object"
  return(l)
}

#' romicsChangeFactor()
#' @description Change the main factor of the romics_object
#' @param romics_object A object created using the function romicsCreateObject().
#' @param main_factor Either 'none' OR a factor from the romics_object (corresponding to a row from the original metadata file), the list of factors from an romics object can be obtained using the function romicsFactorNames().
#' @details changes the main_factor of an romics_object and updates the colors to this new factor.
#' @return This function returns a modified romics object, please see the create_romics_object() documentation.
#' @author Geremy Clair
#' @export
romicsChangeFactor<- function(romics_object , main_factor = "none" ) {
  arguments<-as.list(match.call())
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}

  if(missing(main_factor) | main_factor=="none"){
    warning("The first row of your metadata was used as factor")
    warning(paste0("main_factor='",romicsFactorNames(romics_object)[1]),"'")
    main_factor<-romicsFactorNames(romics_object)[1]
  }

  if(!main_factor %in% romicsFactorNames(romics_object)){
    warning("Your main_factor has to be present in your metadata, the following list of factor are in the romics_object:")
    print(romicsFactorNames(romics_object))
  }

  data<-romics_object$data
  metadata<-romics_object$metadata
  missingdata<-romics_object$missingdata

  #update romics_object$main_factor
  romics_object$main_factor<-main_factor

  #extract the main factor
  main_factor<-  romicsExtractFactor(romics_object,factor = main_factor)
  #establish the lvls of the main factor
  main_factor_lvl<-levels(main_factor)

  tp_data<-data.frame(matrix(nrow=nrow(data),ncol=0))
  tp_metadata<-data.frame(matrix(nrow=nrow(metadata),ncol=0))
  tp_missingdata<-data.frame(matrix(nrow=nrow(missingdata),ncol=0))
  tp_main_factor<-character()

  #reorder data, metadata and missingdata based on the factors which will make visualization more pretty
  for(i in 1:length(main_factor_lvl))
  {
    tp_data<-cbind(tp_data,data[,main_factor==as.character(main_factor_lvl[i])])
    tp_metadata<- cbind(tp_metadata,metadata[,main_factor==as.character(main_factor_lvl[i])])
    tp_missingdata<-cbind(tp_missingdata,missingdata[,main_factor==as.character(main_factor_lvl[i])])
    tp_main_factor<- c(tp_main_factor,main_factor[main_factor==as.character(main_factor_lvl[i])])
  }

  data<-tp_data
  metadata<-tp_metadata
  missingdata<-tp_missingdata
  main_factor<-tp_main_factor
  remove(tp_data,tp_metadata,tp_missingdata)

  #remove the previous color line in the metadata
  metadata<-metadata[rownames(metadata)!="colors_romics",]

  #establish a custom_color vector of the same lenght of the number of levels
  custom_colors<-romics_object$custom_colors

  #add a color_romics line in metadata
  if(length(custom_colors)<length(main_factor_lvl)){
    warning("your color vector is shorter than the number of factors selected, some colors will be picked automatically")
    custom_colors<- c(custom_colors,ROP_colors)
  }

  colors_romics<-data.frame(matrix(ncol=ncol(metadata),nrow=0))
  colors_romics<-rbind(colors_romics,custom_colors[as.numeric(main_factor)])
  colnames(colors_romics)<-colnames(metadata)
  metadata<-rbind(metadata, colors_romics)
  rownames(metadata)[nrow(metadata)]<-"colors_romics"

  #create a vector containing the colors for the whole data points and store it in colors layer
  colors_romics<- t(colors_romics)
  fill <- character(length = 0)
  for (i in 1:length(colors_romics))
  {
    fill<- c(fill,rep(as.character(colors_romics[i]),nrow(data)))
  }
  colors<-fill
  remove(fill,colors_romics, main_factor_lvl,i)

  #replace the different transformed layer in the romics_object
  romics_object$data<-data
  romics_object$metadata<-metadata
  romics_object$missingdata<-missingdata
  romics_object$colors<-colors

  #Update the steps
  romics_object<-romicsUpdateSteps(romics_object,arguments)

  #return the romics_object
  return(romics_object)
}

#' romicsSubset()
#' @description Keeps or drop a subset of specific elements/columns from the romics_object
#' @param romics_object A romics_object created using romicsCreateObject()
#' @param subset_vector A character vector of factor levels or colnames to keep in the object.
#' @param type Either 'keep' or 'drop' to indicate if you want to conserve or to drop the elements from a given factor
#' @param by Either 'colnames' or 'level' to indicate what you want to keep or drop.
#' @param factor A factor contained in the metadata of the romics_object, to obtain the list of factors please use the function romicsFactorNames()
#' @details This function create a new object based on a previous romics_object to include or drops a list of specified columns from the original object. The created object will have a new step object created that will indicate the name of the original object to be subsetted and the log/non-log status of the object.
#' @details Note that this function will remove the stat layer from your object
#' @return This function generate a subseted romics_object
#' @author Geremy Clair
#' @export
romicsSubset<-function(romics_object, subset_vector,type= "keep", by= "colnames", factor="main"){
  arguments<-as.list(match.call())
 # if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(subset_vector)){stop("Your subseting vector is missing")}
  if(!is.character(subset_vector)){stop("Your subset vector has to be a character vector")}
  if(missing(type)){type="keep"}
  if(!type %in% c("keep", "drop")){stop("<type> has to be either 'keep' or 'drop'")}
  if(missing(by)){by="colnames"}
  if(!by %in% c("colnames", "level")){stop("<by> has to be either 'colnames' or 'level'")}
  if(missing(factor)){factor="main"}

  if(factor=="main"){
    factor<-romics_object$main_factor
  }else{
    if(!factor %in% rownames(romics_object$metadata)){stop("Your factor is not present in the metadata")}}

  #extract factor
  fac<-romics_object$metadata[factor==rownames(romics_object$metadata),]
  fac<-as.factor(t(fac))

  #create a logical vector containing the columns to keep (TRUE)
  if(by=="colnames"){
    if(sum(subset_vector %in% colnames(romics_object$data))==length(subset_vector)){warnings("not all the elements of the subset_vector were present in the colnames of the data")}
    if(type=="keep"){
      sub_logical <-  colnames(romics_object$data) %in% subset_vector
    }else{
      sub_logical <-  !colnames(romics_object$data) %in% subset_vector
    }
  }else{
    if(sum(subset_vector %in% fac)==length(subset_vector)){warnings("not all the elements of the subset_vector were levels of the factor ")}
    if(type=="keep"){
      sub_logical <- fac %in% subset_vector
    }else{
      sub_logical <- !fac %in% subset_vector
    }
  }

  #Remove the columns from the element data, metadata, and missingness of the romics_object
    romics_object$data<-romics_object$data[,sub_logical]
    romics_object$metadata<-romics_object$metadata[,sub_logical]
    romics_object$missingdata<-romics_object$missingdata[,sub_logical]

  #remove the stat layer
    if("statistics" %in% names(romics_object)){
      warning("The statistics layer was removed from the romics_object the statistics were calculated on the non-subsetted object")
      romics_object<-romics_object["statistics" != names(romics_object)]}

  #update the colors
    romics_object<-romicsUpdateColor(romics_object)
  #update the steps
    romics_object<-romicsUpdateSteps(romics_object,arguments)

    return(romics_object)
}

#' romicsSampleNameFromFactor()
#' @description Changes the samples identifiers using values contained in one of the factors of the metadata of the romics_object
#' @param romics_object A romics_object created using romicsCreateObject()
#' @param factor A factor contained in the metadata of the romics_object, to obtain the list of factors please use the function romicsFactorNames().
#' @param original_name Either 'keep' or 'drop', this indicate what to do with the original name either store it in a row of the metadata, OR drop it completely.
#' @param colname_factor The name of a factor from the romics object of which the values are unique for each sample. the factors names can be obtain using the function romicsFactorNames()
#' @details enaables the quick renaming of the sample names from a factor contained in metadata, the factor has to contain only unique values.
#' @return a romics_object with its columns of the layers data, metadata, and missingdata renamed.
#' @author Geremy Clair
#' @export
romicsSampleNameFromFactor<-function(romics_object, factor, original_colnames="keep", factor_name="original_colnames"){
  arguments<-as.list(match.call())
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(!factor %in% romicsFactorNames(romics_object)){print("Your Factor should be a row from your metadata")}
  if(sum(duplicated(romics_object$metadata[rownames(romics_object$metadata)==factor,]))>0){
    print("Your selected <factor> was containing duplicated values")
    print(romics_object$metadata[rownames(romics_object$metadata)==factor,])
    stop()
  }
  if(!original_colnames %in% c("keep", "drop")){
    print("You did not selected properly if you wanted to keep or drop the original colnames, by defaults those will be kept")
    original_colnames="keep"
  }
  if(original_colnames %in% c("keep", "drop") & (!is.character(factor_name)) & length(factor_name)!=1){stop("Your <factor_name> as to be a character object of lenght 1.")}
  if(missing(factor_name)){factor_name="original_colnames"}

  new_names<-as.character(t(romics_object$metadata[rownames(romics_object$metadata)==factor,]))

 if(original_colnames=="keep"){
 original_colnames<- data.frame(t(colnames(romics_object$metadata)))
 colnames(original_colnames)<-colnames(romics_object$metadata)
 rownames(original_colnames)<-factor_name
 romics_object$metadata<- rbind(original_colnames,romics_object$metadata)
 }

 colnames(romics_object$data)<- colnames(romics_object$metadata)<-colnames(romics_object$missingdata)<-new_names

 #update the steps
 romics_object<-romicsUpdateColor(romics_object)
 romics_object<-romicsUpdateSteps(romics_object,arguments)

 return(romics_object)
}

#' romicsFactorNames()
#' @description Indicates the list of factor names from the romics_object
#' @param romics_object A romics_object created using romicsCreateObject()
#' @details This function allows to quickly get a vector containing all the factor names present in a romics_object
#' @return A character vector containing the list of factor contained in an romics_object
#' @author Geremy Clair
#' @export
romicsFactorNames<-function(romics_object){
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  rownames(romics_object$metadata)
}

#' romicsExtractFactor()
#' @description Extract a factor from the romics_object
#' @param romics_object A romics_object created using romicsCreateObject()
#' @param factor A factor contained in the romics_object, the list of factors can be obtained using the function romicsFactorNames()
#' @details This function allows to quickly extract the content of a factor present in the romics_object.
#' @return a factor contained in an romics_object the order is the same as the columns of the romics_object$data.
#' @author Geremy Clair
#' @export
romicsExtractFactor<-function(romics_object, factor = "factor"){
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(factor) | !factor %in% romicsFactorNames(romics_object)){stop("The <factor> is missing or is not in the list of factors for this romics_object. The list of available factors can be obtained with the function romicsFactorNames()")}
  fact<-as.factor(t(romics_object$metadata[romicsFactorNames(romics_object)==factor,]))
  names(fact)<-colnames(romics_object$metadata)
  return(fact)
  }

#' romicsUpdateColor()
#' @description Updates the colors layer contained in the romics_object
#' @param romics_object A romics_object created using romicsCreateObject()
#' @return This function returns a romics_object with updated colors.
#' @author Geremy Clair
#' @export
romicsUpdateColor<- function(romics_object) {
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}

  #get the colors from the metadata file
  colors_romics<- as.character(t(romics_object$metadata[grepl("colors_romics",rownames(romics_object$metadata)),]))

  fill <- character(length = 0)
  for (i in 1:length(colors_romics))
  {
    fill<- c(fill,rep(as.character(colors_romics[i]),nrow(romics_object$data)))
  }

  romics_object$colors<-fill

  return(romics_object)

}

#' stepUpdater()
#' @description Updates the steps of the romics_object, require to have recorded the argument in earlier steps of the function
#' @param romics_object A romics_object created using romicsCreateObject()
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

#' is.romics_object()
#' @description Enables to check if the romics_object is in the appropriate format
#' @param romics_object A romics_object created using romicsCreateObject()
#' @return This function will return TRUE or FALSE indicating if the object is or not an romics_object
#' @author Geremy Clair
#' @export
is.romics_object<-function(romics_object){
  if(class(romics_object)!="romics_object"){
    warning("Your romics_object was not created using the function romicsCreateObject")
    return(FALSE)
    }
  if(romics_object$steps[1]!="romics_object"){
    warning("Your romics_object is not in the appropriate format")
    return(FALSE)
  }else{
      return(TRUE)
    }
}

#' romicsLogCheck()
#' @description Identifies if the romics_object is log_transformed or not
#' @param romics_object A romics_object created using romicsCreateObject()
#' @return This function will return TRUE or FALSE indicating if the object was or not log transformed using the function log2transform
#' @author Geremy Clair
#' @export
romicsLogCheck<-function(romics_object){
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(sum(grepl("log2transform\\(",romics_object$steps))>0 | sum(grepl("log10transform\\(",romics_object$steps))>0){return(TRUE)}else{return(FALSE)}
}

#' romicsCreateDependencies()
#' @description Creates the original dependencies of the romics_object (only when it is created to add dependencies use the romicsAddDependency() function)
#' @return This function will return a data.frame with two columns the required packages and the version at the time the code was run
#' @author Geremy Clair
#' @export
romicsCreateDependencies<-function(){
  Required<-data.frame(Required=as.character(getDependencies("RomicsProcessor")), Version_used=NA)
  for(i in 1:nrow(Required)){
    Required$Version_used[i]<- as.character(packageVersion(Required$Required[i]))
  }
  Required[,1]<-as.character(Required[,1])

  Required<-rbind(Required,c( Required= "r", Version_used=paste0(R.Version()$major,".",R.Version()$minor)))
  Required<-rbind(c( Required= "RomicsProcessor", Version_used=as.character(packageVersion("RomicsProcessor"))),Required )
  }

#' romicsAddDependency()
#' @description Adds a package to the list of dependencies of the romics_object. Enables developpers to add automatically a dependency when their function has been applied by the user on their romics_object.
#' @param romics_object A romics_object created using romicsCreateObject()
#' @param new_dependency The name of one or more R package
#' @return This function will return an romics_object updated with the new dependency.
#' @author Geremy Clair
#' @export
romicsAddDependency<-function(romics_object,new_dependency=c("package_1", "package_2")){
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
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
#' @param romics_object A romics_object created using romicsCreateObject()
#' @return This function will return character vector containing the stat columns previously calculated for a romics_object, if no stats were previously calculated an error message will be displayed
#' @author Geremy Clair
#' @export
romicsCalculatedStats<-function(romics_object){
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(is.null(romics_object$statistics)){
   warning("No statistics were calculated for this romics_object")
    return(FALSE)
     }else{
      return(colnames(romics_object$statistics))}
  }

#' romicsSteps()
#' @description Displays the content steps layer of the romics_object
#' @param romics_object A romics_object created using romicsCreateObject()
#' @param show_dates Boolean indicating if the dates have to be displayed
#' @param show_details Boolean indicating if the details have to be displayed
#' @return This function will return the steps of an romics_object
#' @export
romicsSteps<-function(romics_object, show_dates=TRUE, show_details=TRUE){
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(show_dates)){show_dates=TRUE}
  if(!show_dates %in% c(TRUE,FALSE)){stop("show_dates has to be either TRUE or FALSE")}
  if(missing(show_details)){show_details=TRUE}
  if(!show_details %in% c(TRUE,FALSE)){stop("show_details has to be either TRUE or FALSE")}

  steps<-romics_object$steps[2:length(romics_object$steps)]
  if(show_dates==FALSE){steps[grepl("date\\|",steps)]<- gsub(".*\\|","", steps[grepl("date\\|",steps)])}
  if(show_details==FALSE){steps<-steps[!grepl("fun\\|",steps)]}
  return (steps)
}

#' romicsCreatePipeline()
#' @description Extracts a pipeline from the romics_object. the pipeline can then be saved in a classic R object.
#' @param romics_object A romics_object created using romicsCreateObject()
#' @return This function will return character vector containing the stat columns previously calculated for a romics_object, if no stats were previously calculated an error message will be displayed
#' @author Geremy Clair
#' @export
romicsCreatePipeline<- function(romics_object){
  arguments<-as.list(match.call())
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}

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

#' romicsApplyPipeline()
#' @description Applies a pipeline created with the romicsCreatePipeline() function, pipelines can be edited prior to run them.
#' @param romics_object A romics_object created using romicsCreateObject()
#' @param romics_pipeline A pipeline created with the romicsCreatePipeline() function.
#' @return This function will return an romics_object that has been processed through a pipeline
#' @export
romicsApplyPipeline<- function(romics_object, romics_pipeline){
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(typeof(romics_pipeline)!= "character" |sum(substr(romics_pipeline, 1, 15)=="romics_object<-")!=length(romics_pipeline)){stop("Your pipeline will not be applied to the romics_object. Please check the text of your pipeline")}

  for(i in 1:length(romics_pipeline)){
    eval(parse(text=romics_pipeline[i]))
  }

  return(romics_object)
}

#' romicsOutputData()
#' @description Creates an exportable data frame from the romics_object. The generated data.frame contains the processed data, the statistics and the missing status of the data on demand.
#' @param romics_object A romics_object created using romicsCreateObject()
#' @param statistics boolean, has to be TRUE or FALSE. Indicates if the statistics should be exported along with the data (FALSE by default).
#' @param missing_data boolean, has to be TRUE or FALSE. Indicates if the missing status of the data should be exported along with the data (FALSE by default).
#' @return This function will return an data frame containing the results of the processing, the statistics and the missingness status of the data as specified by the user.
#' @author Geremy Clair
#' @export
romicsExportData<-function(romics_object, statistics = FALSE, missing_data = FALSE){
  if(!is.romics_object(romics_object) | missing(romics_object)) {stop("romics_object is missing or is not in the appropriate format")}
  if(missing(statistics)){statistics = FALSE}
  if(missing(missing_data)){missing_data = FALSE}
  if(!is.logical(statistics)){stop("<statistics> should be either TRUE or FALSE")}
  if(!is.logical(missing_data)){stop("<statistics> should be either TRUE or FALSE")}

  df<-romics_object$data

  if(statistics==TRUE){
    if(is.null(romics_object$statistics)){stop("The selected romics object does not contain a 'statistics' layer")}else{df<-cbind(df,romics_object$statistics)}
    }

  if(missing_data==TRUE){
    md <- romics_object$missingdata
    colnames(md)<-paste0("missing_data_",colnames(md))
    df<-cbind(df,md)
    }

  return(df)
}
