#' Example_romics_data
#'
#' This example is containing proteomics results of Bacillus cereus grown in 4 different media:
#' - Luria Bertani (lb)
#' - modified AOAC (maoac)
#' - autoclaved soil extract from INRAE avignon France (soil)
#' - autloclaved zucchinni puree (zucchinni)
#' 
#' All raw data are available using on MassIVE under the identifier MSV000085696 (ftp://massive.ucsd.edu/MSV000085696/)
#' The data was analyzed in MaxQuant and the iBAQ data was extracted using the function extractMaxQuant() on the proteinGroup.txt file
#' The first column contain protein identifiers while the other columns are the iBAQ intensities for the different samples
#' The row represent the different proteins identified.
#'
#' @format A data frame with 2459 rows and 13 columns
#'
"Example_romics_data"

#' Example_romics_metadata
#'
#' This example is the metadata in data.frame format for the Example_romics_data
#' This dataset contains proteomics results of Bacillus cereus grown in 4 different media:
#' - Luria Bertani (lb)
#' - modified AOAC (maoac)
#' - autoclaved soil extract from INRAE avignon France (soil)
#' - autloclaved zucchinni puree (zucchinni)
#' 
#' Metadata in romics have to be imputed in the following format:
#' The first column should contain the list of factors known/used for the analysis
#' Each column correspond to a given sample.
#' Note that the colnames of the metadata have to exactly match the ones of the data except for the first one.
#'
#' @format A data frame with 2 rows and 13 columns
#'
"Example_romics_metadata"

#' Example_unprocessed_romics_object
#'
#' This dataset contains proteomics results of Bacillus cereus grown in 4 different media:
#' - Luria Bertani (lb)
#' - modified AOAC (maoac)
#' - autoclaved soil extract from INRAE avignon France (soil)
#' - autloclaved zucchinni puree (zucchinni)
#' 
#' This romics_object was created by using the function :
#' romicsCreateObject(Example_romics_data,Example_romics_metadata, main_factor = "media")
#'
#' @format A romics_object
#'
"Example_processed_romics_object"

#' Example_romics_processed_romics_dataset
#'
#' This dataset contains proteomics results of Bacillus cereus grown in 4 different media:
#' - Luria Bertani (lb)
#' - modified AOAC (maoac)
#' - autoclaved soil extract from INRAE avignon France (soil)
#' - autloclaved zucchinni puree (zucchinni)
#' 
#' This romics_object was created by using the function :
#' romicsCreateObject(Example_romics_data,Example_romics_metadata, main_factor = "media")
#' and transformed using different romics functions.
#' to see the transformative step employed, use the following code:
#' romicsCreatePipeline(Example_processed_romics_object)
#' 
#' @format romics_object
#'
"Example_processed_romics_object"
