% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/04_Manage_Missing.R
\name{romicsDeleteRowMissing}
\alias{romicsDeleteRowMissing}
\title{romicsDeleteRowMissing()}
\usage{
romicsDeleteRowMissing(romics_object)
}
\arguments{
\item{romics_object}{has to be an romics_object created using romicsCreateObject()}
}
\value{
The function will return a filtered romics_object with the rows of the data and missing data object removed when appropriate.
}
\description{
Removes the rows that have any missing values from the romics_object
}
\details{
This function will remove rows that have any missing values in the data and missingdata layers
}
\author{
Geremy Clair
}
