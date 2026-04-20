#' extractMaxQuant()
#' @description Extracts the quantification information from a MaxQuant ProteinGroup.txt file
#' @param file This has to be the file path and file name of the maxQuant proteinGroup.txt file from which the information has to be extracted
#' @param quantification_type has to be one of the following options : 'LFQ','Intensity','iBAQ','MS.MS', "TMT","TMT.corrected" Indicate what type of quantification needs to be extracted from the ProteinGroup table, can be either
#' @param remove_contaminants has to be TRUE or FALSE, indicates if the contaminant have to be removed
#' @param remove_site_only has to be TRUE or FALSE, indicates if proteins identified by site only have to be removed
#' @param remove_reverse_matches has to be TRUE or FALSE, indicates if the False Positive entries have to be removed
#' @param min_peptides has to be numerical, indicates the minimum number of peptides for protein imporation
#' @details This function will extracts the quantification information from one quantification time of a MaxQuant ProteinGroup.txt file.
#' @return it will return a data.frame with a first column containing the Protein IDs as first column, the other columns will be the Quantitative columns corresponding to the quantitation mode selected.
#' @author Geremy Clair
#' @export
extractMaxQuant<-function(file= "/filepath/proteinGroups.txt",quantification_type="LFQ", remove_contaminants=TRUE, remove_site_only=TRUE, remove_reverse_matches=TRUE, min_peptides=1, min_unique_peptides=1, min_razor_peptides=1){
# ensure that the file is named ProteinGroups.txt
if(substr(file,nchar(file)-16, nchar(file))!="proteinGroups.txt"){stop("The specified file is not a ProteinGroup file")}
if(!quantification_type %in% c("LFQ","Intensity","iBAQ","MS.MS", "Identification.Type","TMT", "TMT.corrected")){stop("The quantification_type is not appropriate")}
if(missing(remove_contaminants)){remove_contaminants<-TRUE}
if(missing(remove_site_only)){remove_site_only<-TRUE}
if(missing(remove_reverse_matches)){remove_reverse_matches<-TRUE}
if(missing(min_peptides)){min_peptides=1}
if(!is.numeric(min_peptides) & !is.numeric(min_peptides)){stop("'min_peptides' has to be numerical and of lenght 1 indicating the proteins with a given minimum of peptides.")}
if(missing(min_unique_peptides)){min_unique_peptides=1}
if(!is.numeric(min_unique_peptides) & !is.numeric(min_unique_peptides)){stop("'min_razor_peptides' has to be numerical and of lenght 1 indicating the proteins with a given minimum of peptides.")}
if(missing(min_razor_peptides)){min_razor_peptides=1}
if(!is.numeric(min_razor_peptides) & !is.numeric(min_razor_peptides)){stop("'min_razor_peptides' has to be numerical and of lenght 1 indicating the proteins with a given minimum of peptides.")}

# read the protein group file
  prGR <- read.delim(file)
# lower the case of the column names
  colnames(prGR)<- tolower(colnames(prGR))

# remove proteins with less than the min peptide count
if (min_peptides>1){
  pep_count<-prGR$peptide.counts..all.
  pep_count<-sub("\\;.*","",pep_count)
  pep_count<-as.numeric(pep_count)
  print(paste0(sum(pep_count<=min_peptides)," protein had less than ",min_peptides," peptides and were removed."))
  prGR<-prGR[pep_count>=min_peptides,]
}

if (min_unique_peptides>1){
  pep_count<-prGR$peptide.counts..razor.unique.
  pep_count<-sub("\\;.*","",pep_count)
  pep_count<-as.numeric(pep_count)
  print(paste0(sum(pep_count<=min_unique_peptides)," protein had less than ",min_unique_peptides," peptides and were removed."))
  prGR<-prGR[pep_count>=min_unique_peptides,]
}

if (min_razor_peptides>1){
  pep_count<-prGR$peptide.counts..razor.unique.
  pep_count<-sub("\\;.*","",pep_count)
  pep_count<-as.numeric(pep_count)
  print(paste0(sum(pep_count<=min_razor_peptides)," protein had less than ",min_razor_peptides," peptides and were removed."))
  prGR<-prGR[pep_count>=min_razor_peptides,]
}

#set the list of ids to keep
  listIDs <-c("protein.ids")

#Create a counts df
  counts <- data.frame(total.protein.count=nrow(prGR),
                  only.identified.by.site= sum(as.character(prGR$only.identified.by.site)=="+",na.rm=T),
                  potential.contaminants= sum(as.character(prGR$potential.contaminant)=="+",na.rm=T),
                  reverse= sum(as.character(prGR$reverse)=="+",na.rm=T))

#Remove the proteins identified by site potential contaminants and reverse hits
if(counts$only.identified.by.site>0&&remove_site_only==TRUE){prGR <- prGR[prGR$only.identified.by.site!="+",]}
if(counts$potential.contaminants>0&&remove_contaminants==TRUE){prGR <- prGR[as.character(prGR$potential.contaminant)!="+",]}
if(counts$reverse>0&&remove_reverse_matches==TRUE){prGR <- prGR[as.character(prGR$reverse)!="+",]}

#Protein Count after removal
  counts <- cbind(counts,after.removal=nrow(prGR))

  type.rm<-"("
  if(remove_site_only==TRUE){type.rm<-paste0(type.rm,"protein(s) only identified by site,")}
  if(remove_contaminants==TRUE){type.rm<-paste0(type.rm,"contaminant(s),")}
  if(remove_reverse_matches==TRUE){type.rm<-paste0(type.rm,"reverse hit(s),")}
  type.rm<- substr(type.rm,1,nchar(type.rm)-1)
  type.rm<-paste0(type.rm,")")

  if(remove_site_only+remove_contaminants+remove_reverse_matches==0){
    print("All proteins including contaminants, reverse hits and site only were conserved")
    }else{
    print(paste(counts$total.protein.count-counts$after.removal," Proteins were removed",type.rm))
    }

#store the names of the groups in a proteinIDs df
  proteinIDs <- prGR[,colnames(prGR) %in% listIDs]

#Keep only the columns corresponding the the appropriate quantification type
  if(quantification_type=="iBAQ"){
    protein_quantification <- prGR[,grepl("ibaq.",names(prGR))]
    colnames(protein_quantification) <- gsub("ibaq.","iBAQ.",colnames(protein_quantification))
    print("iBAQ quantification was used")
  }

  if(quantification_type=="LFQ"){
    protein_quantification <- prGR[,grepl("lfq.intensity.",names(prGR))]
    colnames(protein_quantification) <- gsub("lfq.intensity.","",colnames(protein_quantification))
    print("LFQ quantification was used")
  }

  if(quantification_type=="Intensity"){
    protein_quantification <- prGR[,grepl("intensity.",names(prGR))]
    colnames(protein_quantification) <- gsub("intensity.","",colnames(protein_quantification))
    print("Intensities were used")
  }

  if(quantification_type=="MS.MS"){
    protein_quantification <- prGR[,grepl("ms.ms.",names(prGR))]
    colnames(protein_quantification) <- gsub("ms.ms.","",colnames(protein_quantification))
    print("Spectral count was used")
  }

  if(quantification_type=="Identification.Type"){
    protein_quantification <- prGR[,grepl("identification.type.",names(prGR))]
    colnames(protein_quantification) <- gsub("identification.type.","",colnames(protein_quantification))
    print("The resulting object is containing the identification type")
  }

  if(quantification_type=="TMT.corrected"){
    protein_quantification <- prGR[,grepl("reporter.intensity.corrected.",names(prGR))]
    colnames(protein_quantification) <- gsub("reporter.intensity.corrected.","",colnames(protein_quantification))
    print("TMT reporter intensity was used.")
    }

  if(quantification_type=="TMT"){
    protein_quantification <- prGR[,grepl("reporter.intensity.",names(prGR))]
    protein_quantification <-  protein_quantification[,!grepl("reporter.intensity.corrected.|reporter.intensity.count.",names( protein_quantification))]
    colnames(protein_quantification) <- gsub("reporter.intensity.","",colnames(protein_quantification))
    print("TMT reporter intensity was used.")
  }

  return(cbind(proteinIDs,protein_quantification))

}

#' extractMaxQuantIDs()
#' @description Extracts the IDs information from a MaxQuant ProteinGroup.txt file
#' @param file This has to be the file path and file name of the maxQuant proteinGroup.txt file from which the information has to be extracted
#' @param remove_contaminants has to be TRUE or FALSE, indicates if the contaminant have to be removed
#' @param remove_site_only has to be TRUE or FALSE, indicates if proteins identified by site only have to be removed
#' @param remove_reverse_matches has to be TRUE or FALSE, indicates if the False Positive entries have to be removed
#' @details This function will extracts the IDs information from one quantification time of a MaxQuant ProteinGroup.txt file.
#' @return it will return a data.frame with a first column containing the Protein IDs the following columns will be the following : 'majority.protein.ids', 'fasta.headers', 'peptide.counts.all','peptide.counts.razor.unique','peptide.counts..unique.','fasta.headers','number.of.proteins','peptides','razor...unique.peptides', 'unique.peptides'
#' @author Geremy Clair
#' @export
extractMaxQuantIDs<-function(file= "/filepath/proteinGroups.txt", remove_contaminants=TRUE, remove_site_only=TRUE, remove_reverse_matches=TRUE, min_peptides=1,min_unique_peptides=1, min_razor_peptides=1){
  #ensure that the file is named ProteinGroups.txt
  if(substr(file,nchar(file)-16, nchar(file))!="proteinGroups.txt"){stop("The specified file is not a ProteinGroup file")}
  if(missing(remove_contaminants)){remove_contaminants<-TRUE}
  if(missing(remove_site_only)){remove_site_only<-TRUE}
  if(missing(remove_reverse_matches)){remove_reverse_matches<-TRUE}
  if(missing(min_peptides)){min_peptides=1}
  if(!is.numeric(min_peptides) & !is.numeric(min_peptides)){stop("'min_peptides' has to be numerical and of lenght 1 indicating the proteins with a given minimum of peptides.")}
  if(missing(min_unique_peptides)){min_unique_peptides=1}
  if(!is.numeric(min_unique_peptides) & !is.numeric(min_unique_peptides)){stop("'min_razor_peptides' has to be numerical and of lenght 1 indicating the proteins with a given minimum of peptides.")}
  if(missing(min_razor_peptides)){min_razor_peptides=1}
  if(!is.numeric(min_razor_peptides) & !is.numeric(min_razor_peptides)){stop("'min_razor_peptides' has to be numerical and of lenght 1 indicating the proteins with a given minimum of peptides.")}

  # read the protein group file
  prGR <- read.delim(file)

  # remove proteins with less than the min peptide count
  if (min_peptides>1){
    pep_count<-prGR$peptide.counts..all.
    pep_count<-sub("\\;.*","",pep_count)
    print(paste0(sum(pep_count<=min_peptides)," protein had less than ",min_peptides," peptides and were removed."))
    prGR<-prGR[pep_count<=min_peptides,]
  }

  if (min_unique_peptides>1){
    pep_count<-prGR$peptide.counts..razor.unique.
    pep_count<-sub("\\;.*","",pep_count)
    print(paste0(sum(pep_count<=min_unique_peptides)," protein had less than ",min_unique_peptides," peptides and were removed."))
    prGR<-prGR[pep_count<=min_unique_peptides,]
  }

  if (min_razor_peptides>1){
    pep_count<-prGR$peptide.counts..razor.unique.
    pep_count<-sub("\\;.*","",pep_count)
    print(paste0(sum(pep_count<=min_razor_peptides)," protein had less than ",min_razor_peptides," peptides and were removed."))
    prGR<-prGR[pep_count<=min_razor_peptides,]
  }

  #lower the case of the column names
  colnames(prGR)<- tolower(colnames(prGR))
  #set the list of ids to keep
  listIDs <-c("protein.ids",
              "majority.protein.ids",
              "fasta.headers",
              "peptide.counts..all.",
              "peptide.counts..razor.unique.",
              "peptide.counts..unique.",
              "fasta.headers",
              "number.of.proteins",
              "peptides",
              "razor...unique.peptides",
              "unique.peptides")

  #Create a counts df
  counts <- data.frame(total.protein.count=nrow(prGR),
                       only.identified.by.site= sum(as.character(prGR$only.identified.by.site)=="+",na.rm=T),
                       potential.contaminants= sum(as.character(prGR$potential.contaminant)=="+",na.rm=T),
                       reverse= sum(as.character(prGR$reverse)=="+",na.rm=T))

  #Remove the proteins identified by site potential contaminants and reverse hits
  if(counts$only.identified.by.site>0&&remove_site_only==TRUE){prGR <- prGR[prGR$only.identified.by.site!="+",]}
  if(counts$potential.contaminants>0&&remove_contaminants==TRUE){prGR <- prGR[as.character(prGR$potential.contaminant)!="+",]}
  if(counts$reverse>0&&remove_reverse_matches==TRUE){prGR <- prGR[as.character(prGR$reverse)!="+",]}

  #Protein Count after removal
  counts <- cbind(counts,after.removal=nrow(prGR))
  type.rm<-"("
  if(remove_site_only==TRUE){type.rm<-paste0(type.rm,"protein(s) only identified by site,")}
  if(remove_contaminants==TRUE){type.rm<-paste0(type.rm,"contaminant(s),")}
  if(remove_reverse_matches==TRUE){type.rm<-paste0(type.rm,"reverse hit(s),")}
  type.rm<- substr(type.rm,1,nchar(type.rm)-1)
  type.rm<-paste0(type.rm,")")

  if(remove_site_only+remove_contaminants+remove_reverse_matches==0){
    print("All proteins including contaminants, reverse hits and site only were conserved")
  }else{
    print(paste(counts$total.protein.count-counts$after.removal," Proteins were removed",type.rm))
  }

  #store the names of the groups in a proteinIDs df
  proteinIDs <- prGR[,colnames(prGR) %in% listIDs]

  return(proteinIDs)

}

#' extractFragPipeData()
#' @description Extracts the quantification information from a FragPipe combined_protein.csv file
#' @param file This has to be the file path and file name of the FragPipe combined_protein.csv file from which the information has to be extracted
#' @param quantification_type has to be one of the following options : 'MaxLFQ','Intensity','Total Spectral Count','Unique Spectral Count' Indicate what type of quantification needs to be extracted from the combined_protein.csv table
#' @param remove_contaminants has to be TRUE or FALSE, indicates if the contaminants have to be removed
#' @details This function will extracts the quantification information from one quantification time of a FragPipe combined_protein.csv file.
#' @return it will return a data.frame with a first column containing the Protein IDs as first column, the other columns will be the Quantitative columns corresponding to the quantitation mode selected.
#' @author Geremy Clair
#' @export
extractFragPipeData<-function(file= "/filepath/proteinGroups.txt",quantification_type="MaxLFQ", remove_contaminants=TRUE){
  # ensure that the file is named ombined_protein.csv
  #if(substr(file,nchar(file)-16, nchar(file))!="proteinGroups.txt"){stop("The specified file is not a ProteinGroup file")}
  if(!quantification_type %in% c("MaxLFQ","Intensity","Total Spectral Count","Unique Spectral Count")){stop("The quantification_type is not appropriate")}
  if(missing(remove_contaminants)){remove_contaminants<-TRUE}

  if(quantification_type=="MaxLFQ"){quantification_type="MaxLFQ Intensity"}
  # read the protein group file
  prGR <- read.delim(file)
  n<-nrow(prGR)
  colnames(prGR)<-gsub("\\."," ",colnames(prGR))
  print(paste(n,"proteins were detected."))

  if(remove_contaminants==TRUE){
    prGR<-prGR[!grepl("contam_", prGR$Protein),]
    m<-nrow(prGR)
    print(paste(n-m,"contaminants were removed."))
    print(paste(m,"protein remained."))
    }

  #Create a counts df
  count<- data.frame(cbind(prGR[,colnames(prGR)=="Protein ID"],prGR[,grep(quantification_type,colnames(prGR))]))
  colnames(count)<-gsub("\\."," ",colnames(count))
  colnames(count)<-gsub(paste0(" ",quantification_type),"",colnames(count))
  print(paste(ncol(count)-1,"samples had",quantification_type, "values."))

  return(count)
}

#' extractFragPipeIDs()
#' @description Extracts the IDs from a FragPipe combined_protein.csv file
#' @param file This has to be the file path and file name of the FragPipe combined_protein.csv file from which the information has to be extracted
#' @param remove_contaminants has to be TRUE or FALSE, indicates if the contaminants have to be removed
#' @details This function will extracts the quantification information from one quantification time of a FragPipe combined_protein.csv file.
#' @return it will return a data.frame with a first column containing the Protein IDs as first column, the other columns will be the Quantitative columns corresponding to the quantitation mode selected.
#' @author Geremy Clair
#' @export
extractFragPipeIDs<-function(file= "/filepath/proteinGroups.txt",remove_contaminants=TRUE){
  # ensure that the file is named ombined_protein.csv
  #if(substr(file,nchar(file)-16, nchar(file))!="proteinGroups.txt"){stop("The specified file is not a ProteinGroup file")}
  if(missing(remove_contaminants)){remove_contaminants<-TRUE}
  # read the protein group file
  prGR <- read.delim(file)
  n<-nrow(prGR)
  colnames(prGR)<-gsub("\\."," ",colnames(prGR))
  print(paste(n,"proteins were detected."))

  if(remove_contaminants==TRUE){
    prGR<-prGR[!grepl("contam_", prGR$Protein),]
    m<-nrow(prGR)
    print(paste(n-m,"contaminants were removed."))
    print(paste(m,"protein remained."))
  }
  #Create a counts df
  IDs<- data.frame(prGR[,colnames(prGR) %in% c("Protein",
                                               "Protein ID",
                                               "Entry Name",
                                               "Gene",
                                               "Protein Length",
                                               "Organism",
                                               "Protein Existence",
                                               "Description",
                                               "Protein Probability",
                                               "Top Peptide Probability",
                                               "Combined Total Peptides",
                                               "Indistinguishable Proteins")])

  IDs<- cbind(Protein.ID=IDs$Protein.ID, IDs[,colnames(IDs)!="Protein.ID"])
  return(IDs)
}

#' extractDIANN()
#' @description Extracts the quantification information from a DIANN report.pg_matrix.tsv file
#' @param file This has to be the file path and file name of the FragPipe report.pg_matrix.tsv file from which the information has to be extracted
#' @details This function will extracts the quantification information from one quantification time of a DIA-NN report.pg_matrix.tsv file.
#' @return it will return a data.frame with a first column containing the Protein IDs as first column, the other columns will be the Quantitative columns corresponding to the quantitation mode selected.
#' @author Geremy Clair
#' @export
extractDIANN<-function(file= "/filepath/report.pg_matrix.tsv"){
  # read the protein group file
  prGR <- read.delim(file)
  n<-nrow(prGR)
  colnames(prGR)<-gsub("\\."," ",colnames(prGR))
  print(paste(n,"proteins were detected."))

  #Create a counts df
  count<- data.frame(cbind(Protein_Group=prGR[,colnames(prGR)=="Protein Group"],prGR[,!colnames(prGR) %in% c("Protein Group","Protein Ids","Protein Names", "Genes","First Protein Description")]))
  return(count)
}

#' extractDIANNIDs()
#' @description Extracts the IDs information from a DIANN report.pg_matrix.tsv file
#' @param file This has to be the file path and file name of the FragPipe report.pg_matrix.tsv file from which the information has to be extracted
#' @details This function will extracts the quantification information from one quantification time of a DIA-NN report.pg_matrix.tsv file.
#' @return it will return a data.frame with a first column containing the Protein IDs as first column, the other columns will be the Quantitative columns corresponding to the quantitation mode selected.
#' @author Geremy Clair
#' @export
extractDIANNIDs<-function(file= "/filepath/report.pg_matrix.tsv"){
  # read the protein group file
  prGR <- read.delim(file)
  n<-nrow(prGR)
  colnames(prGR)<-gsub("\\."," ",colnames(prGR))
  print(paste(n,"proteins were detected."))

  #Create a counts df
  IDs<- data.frame(prGR[,colnames(prGR) %in% c("Protein Group","Protein Ids","Protein Names", "Genes","First Protein Description")])
  return(IDs)
}

#' extractScils()
#' @description Extracts intensity data, metadata, and feature IDs from a SCiLS '.slx' file.
#' The function efficiently processes large datasets common in mass spectrometry imaging.
#' @param file Character string. Path to the SCiLS '.slx' file.
#' @param scils_executable Character string. Path to the SCiLS executable.
#'        Default: 'C:/Program Files/SCiLS/SCiLS Lab/scilsMsServer.exe'
#' @param port Numeric. Port to use for the connection with the SCiLS environment.
#'        Default: 8082
#' @param feature_list Character string. Name of the feature list to import.
#'        Default: "All Features"
#' @param normId Character string. Normalization ID to apply during import.
#'        Default: "" (no normalization)
#' @param metadata_fields Character vector. Names of metadata fields to extract.
#'        Default: NULL (extracts all available metadata fields)
#' @details
#' This function extracts data from a SCiLS '.slx' file through the SCiLS Lab API.
#' It requires both the '.slx' file and corresponding '.sbd' file to be present in the same directory.
#' The necessary companion '.sbd' file must be in the same directory as the '.slx' file.
#' The function was optimized for handling large datasets with potentially hundreds of thousands of
#' data points. It might still be slow to operate due to some limitation of the SCIls API.
#'
#' It imports: the intensity matrix for the selected features, the Spot IDs and associated metadata, the Detailed information about each feature
#'
#' @return A list with three components:
#'   \item{data}{Data frame containing the intensity matrix with features as rows and samples as columns}
#'   \item{metadata}{Data frame containing spot IDs and other spot attributes}
#'   \item{IDs}{Data frame containing detailed information about each feature}
#'
#'
#' @author Geremy Clair, Brittney Gorman
#' @export
extractScils <- function(file, scils_executable, port, feature_list, normId, metadata_fields = NULL) {
  # Input validation with detailed error messages
  if(missing(file)) {
    stop("Parameter 'file' is required: please provide the path to your SCiLS .slx file.")
  }

  if(!is.character(file) || length(file) != 1) {
    stop("Parameter 'file' must be a character string specifying the path to your SCiLS .slx file.")
  }

  if(!file.exists(file)) {
    stop(sprintf("File not found: '%s'. Please verify the file path is correct.", file))
  }

  # Check for the basedata file
  basedata_file <- gsub("\\.slx$", ".sbd", file)
  if(!file.exists(basedata_file)) {
    stop(sprintf("Required companion file not found: '%s'.\n  This file must be in the same directory as your .slx file.", basedata_file))
  }

  # Set default parameters if missing
  if(missing(scils_executable)) {
    scils_executable <- "C:/Program Files/SCiLS/SCiLS Lab/scilsMsServer.exe"
  }

  if(!file.exists(scils_executable)) {
    stop(sprintf("SCiLS executable not found: '%s'.\n  Please provide the correct path to scilsMsServer.exe.", scils_executable))
  }

  if(missing(port)) {
    port <- 8082
  } else if(!is.numeric(port) || length(port) != 1) {
    stop("Parameter 'port' must be a single numeric value.")
  }

  if(missing(feature_list)) {
    feature_list <- "All Features"
  }

  if(missing(normId)) {
    normId <- ""
  }

  # Ensure proper cleanup on exit or error
  ScilsFileEnv <- NULL
  on.exit({
    if(!is.null(ScilsFileEnv)) {
      message(sprintf("Closing SCiLS Lab session (%s).", Sys.time()))
      try(close(ScilsFileEnv), silent = TRUE)
    }
  })

  # Load the SCiLS environment and import the data
  tryCatch({
    ScilsFileEnv <- SCiLSLabOpenLocalSession(file, port = port, executable = scils_executable,timeout = 120)
    message(sprintf("%s was opened to import the SCiLS data in R (%s).",
                    getServerVersion(ScilsFileEnv), Sys.time()))
  }, error = function(e) {
    stop(sprintf("Failed to open SCiLS session: %s", e$message))
  })

  # Load available feature lists
  tryCatch({
    Flists <- as.character(t(SCiLSLabClient::getFeatureLists(ScilsFileEnv)[, 2]))

    # Verify if the feature_list exists
    if(!feature_list %in% Flists) {
      stop(sprintf("Feature list '%s' not found. Available lists are: %s",
                   feature_list, paste(Flists, collapse = ", ")))
    }

    # Import the list of features
    selected_features <- SCiLSLabClient::getFeatures(
      ScilsFileEnv,
      name = feature_list,
      includeVisibleUserColumns = TRUE
    )
  }, error = function(e) {
    stop(sprintf("Failed to retrieve feature lists: %s", e$message))
  })

  # Import intensity matrix
  message(sprintf("Beginning data matrix importation (%s).", Sys.time()))
  tryCatch({
    data <- do.call(rbind, lapply(selected_features$id, function(x) {
      intensities <- SCiLSLabClient::getFeatureIntensities(
        ScilsFileEnv,
        x,
        regionId = 'Regions',
        normId = normId
      )$intensity
      as.numeric(intensities)  # Convert to numeric
    }))
    message(sprintf("Data matrix imported (%s).", Sys.time()))
  }, error = function(e) {
    stop(sprintf("Failed to import data matrix: %s", e$message))
  })

  # Get spot IDs
  spot_ids <- getRegionSpots(ScilsFileEnv, regionId = 'Regions')$spotId

  # Import metadata using the more efficient approach
  message(sprintf("Beginning metadata importation (%s).", Sys.time()))
  tryCatch({
    # Get region tree for metadata extraction
    regionTree <- getRegionTree(ScilsFileEnv)

    # Ensure regionTree is of the correct class
    if (!"RegionTree" %in% class(regionTree)) {
      stop("RegionTree argument is not of class 'RegionTree'")
    }

    # Flatten the region tree
    allRegions <- flattenRegionTree(regionTree)

    # Initialize spots data frame with coordinates
    spots_data <- allRegions[[1]]$spots
    base_fields <- c("spotId", "x", "y", "z")
    base_fields <- base_fields[base_fields %in% colnames(spots_data)]

    # Extract only base fields initially
    return_df <- spots_data[, base_fields, drop = FALSE]

    # Build attribute dictionary
    attributes <- list()
    for (region in allRegions) {
      # Skip regions with no attributes or only single attributes
      if (is.null(region$attributes) ||
          (is.character(region$attributes$name) && length(region$attributes$name) <= 1)) {
        next
      }

      # Process each attribute in the region
      for (att in seq_len(nrow(region$attributes))) {
        attribute_name <- region$attributes[att, ]$name

        # Skip "Date" attributes
        if (attribute_name == "Date") {
          next
        }

        attribute_level <- region$attributes[att, ]$value

        # Handle "Class" attribute specially
        if (attribute_name == "Class") {
          attribute_name <- attribute_level
          attribute_level <- attribute_name
        }

        # Skip if metadata_fields is specified and this attribute isn't in it
        if (!is.null(metadata_fields) && !(attribute_name %in% metadata_fields)) {
          next
        }

        # Initialize attribute_name and levels if not already in the list
        if (!attribute_name %in% names(attributes)) {
          attributes[[attribute_name]] <- list()
        }

        if (!attribute_level %in% names(attributes[[attribute_name]])) {
          attributes[[attribute_name]][[attribute_level]] <- region$spots$spotId
        } else {
          # Append unique values using union
          attributes[[attribute_name]][[attribute_level]] <- union(
            attributes[[attribute_name]][[attribute_level]], region$spots$spotId
          )
        }
      }
    }

    # If metadata_fields is NULL, use all attributes found
    if (is.null(metadata_fields)) {
      metadata_fields <- names(attributes)
    } else {
      # Ensure we only process fields that actually exist
      metadata_fields <- intersect(metadata_fields, names(attributes))
    }

    # Vectorized function to determine attribute levels for a spot
    get_region_levels <- function(spot_id, attr) {
      levels <- names(attr)[sapply(attr, function(level_spots) {
        spot_id %in% level_spots
      })]
      if (length(levels) == 0) {
        return("absent")
      } else {
        return(paste(levels, collapse = ","))
      }
    }

    # Process each metadata field
    for (attribute_name in metadata_fields) {
      if (!attribute_name %in% names(attributes)) {
        # Add column of "absent" if attribute doesn't exist
        return_df[[attribute_name]] <- "absent"
        next
      }

      attr <- attributes[[attribute_name]]

      # Use vectorized apply for efficiency
      return_df[[attribute_name]] <- vapply(
        return_df$spotId,
        FUN = get_region_levels,
        FUN.VALUE = character(1),
        attr = attr
      )
    }

    # Clean up any NA values
    return_df[is.na(return_df)] <- "absent"

    # Clean up duplicate positive values
    for (col in colnames(return_df)) {
      if (!col %in% c("spotId", "x", "y", "z")) {  # Skip coordinate columns
        # Remove duplicate values in comma-separated strings
        return_df[[col]] <- sapply(return_df[[col]], function(x) {
          if (is.na(x)) return(NA)
          unique_vals <- unique(unlist(strsplit(as.character(x), ",")))
          paste(unique_vals, collapse = ",")
        })
      }
    }

    # Format for consistency with the original function
    metadata <- t(return_df)
    rownames(metadata) <- c(base_fields, metadata_fields)

    message(sprintf("Metadata imported (%s).", Sys.time()))
  }, error = function(e) {
    stop(sprintf("Failed to import metadata: %s", e$message))
  })



  # Process feature IDs
  message(sprintf("Beginning of IDs importation (%s).", Sys.time()))
  tryCatch({
    # Create feature names
    features <- makeFeatureNames(selected_features, useName = FALSE, digits = 5)

    # Handle the "m/z:NA" case
    na_rows <- features == "m/z:NA"
    if (any(na_rows)) {
      features[na_rows] <- makeFeatureNames(selected_features[na_rows, ], useName = TRUE)
    }

    # Drop duplicated rows
    duplicated_mask <- duplicated(features) | (features == "m/z:NA")
    features <- features[!duplicated_mask]
    data <- data[!duplicated_mask, ]

    # Combine features with intensity matrix and set column names
    spot_id_cols <- paste0("spotId_", spot_ids)
    colnames(data) <- spot_id_cols
    data <- as.data.frame(cbind(features, data))

    # Process metadata to match data format
    colnames(metadata) <- spot_id_cols
    metadata <- cbind(data.frame(factor = rownames(metadata)), metadata)

    # Process IDs with deduplication
    IDs <- selected_features[!duplicated_mask, ]
    IDs <- cbind(features, IDs)
    IDs <- IDs[!duplicated(IDs$features), ]
    rownames(IDs) <- seq_len(nrow(IDs))
  }, error = function(e) {
    stop(sprintf("Failed to process feature IDs: %s", e$message))
  })

  # Create the result list
  result <- list(
    data = data,
    metadata = metadata,
    IDs = IDs
  )
  close(ScilsFileEnv)
  message(sprintf("SCiLS data has been successfully extracted (%s).", Sys.time()))

  return(result)
}

#' extractScilsOLD()
#' @description Extracts intensity data, metadata, and feature IDs from a SCiLS '.slx' file.
#' The function efficiently processes large datasets common in mass spectrometry imaging.
#' @param file Character string. Path to the SCiLS '.slx' file.
#' @param scils_executable Character string. Path to the SCiLS executable.
#'        Default: 'C:/Program Files/SCiLS/SCiLS Lab/scilsMsServer.exe'
#' @param port Numeric. Port to use for the connection with the SCiLS environment.
#'        Default: 8082
#' @param feature_list Character string. Name of the feature list to import.
#'        Default: "All Features"
#' @param normId Character string. Normalization ID to apply during import.
#'        Default: "" (no normalization)
#' @param metadata_fields Character vector. Names of metadata fields to extract.
#'        Default: NULL (extracts all available metadata fields)
#' @details
#' This function extracts data from a SCiLS '.slx' file through the SCiLS Lab API.
#' It requires both the '.slx' file and corresponding '.sbd' file to be present in the same directory.
#'
#' The function is optimized for handling large datasets with potentially hundreds of thousands of
#' data points. It imports:
#' \itemize{
#'   \item Intensity matrix for the selected features
#'   \item Spot IDs and associated metadata
#'   \item Detailed information about each feature
#' }
#'
#' The necessary companion '.sbd' file must be in the same directory as the '.slx' file.
#'
#' @return A list with three components:
#' \describe{
#'   \item{data}{Data frame containing the intensity matrix with features as rows and samples as columns}
#'   \item{metadata}{Data frame containing spot IDs and other spot attributes}
#'   \item{IDs}{Data frame containing detailed information about each feature}
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage with default parameters
#' scils_data <- extractScils("path/to/your/file.slx")
#'
#' # Custom parameters
#' scils_data <- extractScils(
#'   file = "path/to/your/file.slx",
#'   scils_executable = "C:/Custom/Path/scilsMsServer.exe",
#'   port = 8083,
#'   feature_list = "My Custom Feature List",
#'   normId = "TIC"
#' )
#'
#' # Access the components
#' data <- scils_data$data
#' metadata <- scils_data$metadata
#' feature_ids <- scils_data$IDs
#' }
#'
#' @author Geremy Clair, Brittney Gorman
#' @export
extractScilsOLD <- function(file, scils_executable, port, feature_list, normId, metadata_fields = NULL) {
  # Input validation with detailed error messages
  if(missing(file)) {
    stop("Parameter 'file' is required: please provide the path to your SCiLS .slx file.")
  }

  if(!is.character(file) || length(file) != 1) {
    stop("Parameter 'file' must be a character string specifying the path to your SCiLS .slx file.")
  }

  if(!file.exists(file)) {
    stop(sprintf("File not found: '%s'. Please verify the file path is correct.", file))
  }

  # Check for the basedata file
  basedata_file <- gsub("\\.slx$", ".sbd", file)
  if(!file.exists(basedata_file)) {
    stop(sprintf("Required companion file not found: '%s'.\n  This file must be in the same directory as your .slx file.", basedata_file))
  }

  # Set default parameters if missing
  if(missing(scils_executable)) {
    scils_executable <- "C:/Program Files/SCiLS/SCiLS Lab/scilsMsServer.exe"
  }

  if(!file.exists(scils_executable)) {
    stop(sprintf("SCiLS executable not found: '%s'.\n  Please provide the correct path to scilsMsServer.exe.", scils_executable))
  }

  if(missing(port)) {
    port <- 8082
  } else if(!is.numeric(port) || length(port) != 1) {
    stop("Parameter 'port' must be a single numeric value.")
  }

  if(missing(feature_list)) {
    feature_list <- "All Features"
  }

  if(missing(normId)) {
    normId <- ""
  }

  # Ensure proper cleanup on exit or error
  ScilsFileEnv <- NULL
  on.exit({
    if(!is.null(ScilsFileEnv)) {
      message(sprintf("Closing SCiLS Lab session (%s).", Sys.time()))
      try(close(ScilsFileEnv), silent = TRUE)
    }
  })

  # Load the SCiLS environment and import the data
  tryCatch({
    ScilsFileEnv <- SCiLSLabOpenLocalSession(file, port = port, executable = scils_executable)
    message(sprintf("%s was opened to import the SCiLS data in R (%s).",
                    getServerVersion(ScilsFileEnv), Sys.time()))
  }, error = function(e) {
    stop(sprintf("Failed to open SCiLS session: %s", e$message))
  })

  # Load available feature lists
  tryCatch({
    Flists <- as.character(t(SCiLSLabClient::getFeatureLists(ScilsFileEnv)[, 2]))

    # Verify if the feature_list exists
    if(!feature_list %in% Flists) {
      stop(sprintf("Feature list '%s' not found. Available lists are: %s",
                   feature_list, paste(Flists, collapse = ", ")))
    }

    # Import the list of features
    selected_features <- SCiLSLabClient::getFeatures(
      ScilsFileEnv,
      name = feature_list,
      includeVisibleUserColumns = TRUE
    )
  }, error = function(e) {
    stop(sprintf("Failed to retrieve feature lists: %s", e$message))
  })

  # Import intensity matrix
  message(sprintf("Beginning data matrix importation (%s).", Sys.time()))
  tryCatch({
    data <- do.call(rbind, lapply(selected_features$id, function(x) {
      intensities <- SCiLSLabClient::getFeatureIntensities(
        ScilsFileEnv,
        x,
        regionId = 'Regions',
        normId = normId
      )$intensity
      as.numeric(intensities)  # Convert to numeric
    }))
    message(sprintf("Data matrix imported (%s).", Sys.time()))
  }, error = function(e) {
    stop(sprintf("Failed to import data matrix: %s", e$message))
  })

  # Get spot IDs
  spot_ids <- getRegionSpots(ScilsFileEnv, regionId = 'Regions')$spotId

  # Import metadata using the more efficient approach
  message(sprintf("Beginning metadata importation (%s).", Sys.time()))
  tryCatch({
    # Get region tree for metadata extraction
    regionTree <- getRegionTree(ScilsFileEnv)

    # Ensure regionTree is of the correct class
    if (!"RegionTree" %in% class(regionTree)) {
      stop("RegionTree argument is not of class 'RegionTree'")
    }

    # Flatten the region tree
    allRegions <- flattenRegionTree(regionTree)

    # Initialize spots data frame with coordinates
    spots_data <- allRegions[[1]]$spots
    base_fields <- c("spotId", "x", "y", "z")
    base_fields <- base_fields[base_fields %in% colnames(spots_data)]

    # Extract only base fields initially
    return_df <- spots_data[, base_fields, drop = FALSE]

    # Build attribute dictionary
    attributes <- list()
    for (region in allRegions) {
      # Skip regions with no attributes or only single attributes
      if (is.null(region$attributes) ||
          (is.character(region$attributes$name) && length(region$attributes$name) <= 1)) {
        next
      }

      # Process each attribute in the region
      for (att in seq_len(nrow(region$attributes))) {
        attribute_name <- region$attributes[att, ]$name

        # Skip "Date" attributes
        if (attribute_name == "Date") {
          next
        }

        attribute_level <- region$attributes[att, ]$value

        # Handle "Class" attribute specially
        if (attribute_name == "Class") {
          attribute_name <- attribute_level
          attribute_level <- "positive"
        }

        # Skip if metadata_fields is specified and this attribute isn't in it
        if (!is.null(metadata_fields) && !(attribute_name %in% metadata_fields)) {
          next
        }

        # Initialize attribute_name and levels if not already in the list
        if (!attribute_name %in% names(attributes)) {
          attributes[[attribute_name]] <- list()
        }

        if (!attribute_level %in% names(attributes[[attribute_name]])) {
          attributes[[attribute_name]][[attribute_level]] <- region$spots$spotId
        } else {
          # Append unique values using union
          attributes[[attribute_name]][[attribute_level]] <- union(
            attributes[[attribute_name]][[attribute_level]], region$spots$spotId
          )
        }
      }
    }

    # If metadata_fields is NULL, use all attributes found
    if (is.null(metadata_fields)) {
      metadata_fields <- names(attributes)
    } else {
      # Ensure we only process fields that actually exist
      metadata_fields <- intersect(metadata_fields, names(attributes))
    }

    # Vectorized function to determine attribute levels for a spot
    get_region_levels <- function(spot_id, attr) {
      levels <- names(attr)[sapply(attr, function(level_spots) {
        spot_id %in% level_spots
      })]
      if (length(levels) == 0) {
        return("unknown")
      } else {
        return(paste(levels, collapse = ","))
      }
    }

    # Process each metadata field
    for (attribute_name in metadata_fields) {
      if (!attribute_name %in% names(attributes)) {
        # Add column of "unknown" if attribute doesn't exist
        return_df[[attribute_name]] <- "unknown"
        next
      }

      attr <- attributes[[attribute_name]]

      # Use vectorized apply for efficiency
      return_df[[attribute_name]] <- vapply(
        return_df$spotId,
        FUN = get_region_levels,
        FUN.VALUE = character(1),
        attr = attr
      )
    }

    # Clean up any NA values
    return_df[is.na(return_df)] <- "unknown"

    # Clean up duplicate positive values
    for (col in colnames(return_df)) {
      return_df[[col]] <- gsub("positive,positive(,positive)*", "positive", return_df[[col]])
    }

    # Format for consistency with the original function
    metadata <- t(return_df)
    rownames(metadata) <- c(base_fields, metadata_fields)

    message(sprintf("Metadata imported (%s).", Sys.time()))
  }, error = function(e) {
    stop(sprintf("Failed to import metadata: %s", e$message))
  })



  # Process feature IDs
  message(sprintf("Beginning of IDs importation (%s).", Sys.time()))
  tryCatch({
    # Create feature names
    features <- makeFeatureNames(selected_features, useName = FALSE, digits = 5)

    # Handle the "m/z:NA" case
    na_rows <- rownames(features) == "m/z:NA"
    if (any(na_rows)) {
      features[na_rows, ] <- makeFeatureNames(selected_features[na_rows, ], useName = TRUE)
    }

    # Drop duplicated rows
    duplicated_mask <- duplicated(features) | (features == "m/z:NA")
    features <- features[!duplicated_mask]
    data <- data[!duplicated_mask, ]

    # Combine features with intensity matrix and set column names
    spot_id_cols <- paste0("spotId_", spot_ids)
    colnames(data) <- spot_id_cols
    data <- as.data.frame(cbind(features, data))

    # Process metadata to match data format
    colnames(metadata) <- spot_id_cols
    metadata <- cbind(data.frame(factor = rownames(metadata)), metadata)

    # Process IDs with deduplication
    IDs <- selected_features[!duplicated_mask, ]
    IDs <- cbind(features, IDs)
    IDs <- IDs[!duplicated(IDs$features), ]
    rownames(IDs) <- seq_len(nrow(IDs))
  }, error = function(e) {
    stop(sprintf("Failed to process feature IDs: %s", e$message))
  })

  # Create the result list
  result <- list(
    data = data,
    metadata = metadata,
    IDs = IDs
  )
  close(ScilsFileEnv)
  message(sprintf("SCiLS data has been successfully extracted (%s).", Sys.time()))

  return(result)
}
