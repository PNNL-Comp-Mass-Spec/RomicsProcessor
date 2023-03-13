#' extractMaxQuant()
#' @description Extracts the quantification information from a MaxQuant ProteinGroup.txt file
#' @param file This has to be the file path and file name of the maxQuant proteinGroup.txt file from which the information has to be extracted
#' @param quantification_type has to be one of the following options : 'LFQ','Intensity','iBAQ','MS.MS', "TMT","TMT.corrected" Indicate what type of quantification needs to be extracted from the ProteinGroup table, can be either
#' @param cont.rm has to be TRUE or FALSE, indicates if the contaminant have to be removed
#' @param site.rm has to be TRUE or FALSE, indicates if the identification by site only have to be removed
#' @param rev.rm has to be TRUE or FALSE, indicates if the False Positive entries have to be removed
#' @param min_peptides has to be numerical, indicates the minimum number of peptides for protein imporation
#' @details This function will extracts the quantification information from one quantification time of a MaxQuant ProteinGroup.txt file.
#' @return it will return a data.frame with a first column containing the Protein IDs as first column, the other columns will be the Quantitative columns corresponding to the quantitation mode selected.
#' @author Geremy Clair
#' @export
extractMaxQuant<-function(file= "/filepath/proteinGroups.txt",quantification_type="LFQ", cont.rm=TRUE,site.rm=TRUE, rev.rm=TRUE, min_peptides=1, min_unique_peptides=1, min_razor_peptides=1){
# ensure that the file is named ProteinGroups.txt
if(substr(file,nchar(file)-16, nchar(file))!="proteinGroups.txt"){stop("The specified file is not a ProteinGroup file")}
if(!quantification_type %in% c("LFQ","Intensity","iBAQ","MS.MS", "Identification.Type","TMT", "TMT.corrected")){stop("The quantification_type is not appropriate")}
if(missing(cont.rm)){cont.rm<-TRUE}
if(missing(site.rm)){site.rm<-TRUE}
if(missing(rev.rm)){rev.rm<-TRUE}
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
if(counts$only.identified.by.site>0&&site.rm==TRUE){prGR <- prGR[prGR$only.identified.by.site!="+",]}
if(counts$potential.contaminants>0&&cont.rm==TRUE){prGR <- prGR[as.character(prGR$potential.contaminant)!="+",]}
if(counts$reverse>0&&rev.rm==TRUE){prGR <- prGR[as.character(prGR$reverse)!="+",]}

#Protein Count after removal
  counts <- cbind(counts,after.removal=nrow(prGR))

  type.rm<-"("
  if(site.rm==TRUE){type.rm<-paste0(type.rm,"protein(s) only identified by site,")}
  if(cont.rm==TRUE){type.rm<-paste0(type.rm,"contaminant(s),")}
  if(rev.rm==TRUE){type.rm<-paste0(type.rm,"reverse hit(s),")}
  type.rm<- substr(type.rm,1,nchar(type.rm)-1)
  type.rm<-paste0(type.rm,")")

  if(site.rm+cont.rm+rev.rm==0){
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
#' @param cont.rm has to be TRUE or FALSE, indicates if the contaminant have to be removed
#' @param site.rm has to be TRUE or FALSE, indicates if the identification by site only have to be removed
#' @param rev.rm has to be TRUE or FALSE, indicates if the False Positive entries have to be removed
#' @details This function will extracts the IDs information from one quantification time of a MaxQuant ProteinGroup.txt file.
#' @return it will return a data.frame with a first column containing the Protein IDs the following columns will be the following : 'majority.protein.ids', 'fasta.headers', 'peptide.counts.all','peptide.counts.razor.unique','peptide.counts..unique.','fasta.headers','number.of.proteins','peptides','razor...unique.peptides', 'unique.peptides'
#' @author Geremy Clair
#' @export
extractMaxQuantIDs<-function(file= "/filepath/proteinGroups.txt", cont.rm=TRUE,site.rm=TRUE, rev.rm=TRUE, min_peptides=1,min_unique_peptides=1, min_razor_peptides=1){
  #ensure that the file is named ProteinGroups.txt
  if(substr(file,nchar(file)-16, nchar(file))!="proteinGroups.txt"){stop("The specified file is not a ProteinGroup file")}
  if(missing(cont.rm)){cont.rm<-TRUE}
  if(missing(site.rm)){site.rm<-TRUE}
  if(missing(rev.rm)){rev.rm<-TRUE}
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
  if(counts$only.identified.by.site>0&&site.rm==TRUE){prGR <- prGR[prGR$only.identified.by.site!="+",]}
  if(counts$potential.contaminants>0&&cont.rm==TRUE){prGR <- prGR[as.character(prGR$potential.contaminant)!="+",]}
  if(counts$reverse>0&&rev.rm==TRUE){prGR <- prGR[as.character(prGR$reverse)!="+",]}

  #Protein Count after removal
  counts <- cbind(counts,after.removal=nrow(prGR))
  type.rm<-"("
  if(site.rm==TRUE){type.rm<-paste0(type.rm,"protein(s) only identified by site,")}
  if(cont.rm==TRUE){type.rm<-paste0(type.rm,"contaminant(s),")}
  if(rev.rm==TRUE){type.rm<-paste0(type.rm,"reverse hit(s),")}
  type.rm<- substr(type.rm,1,nchar(type.rm)-1)
  type.rm<-paste0(type.rm,")")

  if(site.rm+cont.rm+rev.rm==0){
    print("All proteins including contaminants, reverse hits and site only were conserved")
  }else{
    print(paste(counts$total.protein.count-counts$after.removal," Proteins were removed",type.rm))
  }

  #store the names of the groups in a proteinIDs df
  proteinIDs <- prGR[,colnames(prGR) %in% listIDs]

  return(proteinIDs)

}

