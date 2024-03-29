% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/16_Correlation_plot.R
\name{romicsCorrelation}
\alias{romicsCorrelation}
\title{romicsCorrelation()}
\usage{
romicsCorrelation(
  romics_object,
  by_level = FALSE,
  factor = "main",
  imputed = FALSE,
  corr_type = "pearson",
  mode = "scatter",
  ...
)
}
\arguments{
\item{romics_object}{an object of type romics_object generated with the}

\item{by_level}{Boolean indicates if the correlation plot has to be displayed for each level of the factor set with the parameter <factor = >, if multiple levels are present, one correlation plot will be generated by factor, if FALSE the global correlation plot will be generated.}

\item{factor}{Indicates what factor to use for the <by_level> plotting, by default ("main")}

\item{imputed}{Indicates if the imputed values should be used for the correlation (FALSE by default) if TRUE, the data employed will utilize both imputed and non imputed values.}

\item{corr_type}{has to be either 'pearson', 'kendall',or 'spearman', indicates the type of correlation calculated.}

\item{mode}{has to be either 'scatter' or 'colors' to indicate the type of plot to display}

\item{...}{passes argument to the function ggpairs from the package GGally, type ?ggpairs() for documentation}
}
\value{
This function will print different plot requested
}
\description{
generate one or multiple volcano plots from t.test or wilcox.test run using the functions romicsTtest() or romicsWilcoxTest() respectively
}
\details{
generate one or multiple correlation plots using the function ggpairs from the package GGally (scatter mode) or a simple color heatmap (color mode).
}
\examples{
romicsCorrelation(Example_processed_romics_object,by_level=FALSE,factor="main",imputed=T,corr_type="pearson",mode="color")
}
\author{
Geremy Clair
}
