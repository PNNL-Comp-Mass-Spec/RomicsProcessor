% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02_Romics_Base_Functions.R
\name{romicsExtractFactor}
\alias{romicsExtractFactor}
\title{romicsExtractFactor()}
\usage{
romicsExtractFactor(romics_object, factor = "factor")
}
\arguments{
\item{romics_object}{A romics_object created using romicsCreateObject()}

\item{factor}{A factor contained in the romics_object, the list of factors can be obtained using the function romicsFactorNames()}
}
\value{
a factor contained in an romics_object the order is the same as the columns of the romics_object$data.
}
\description{
Extract a factor from the romics_object
}
\details{
This function allows to quickly extract the content of a factor present in the romics_object.
}
\author{
Geremy Clair
}
