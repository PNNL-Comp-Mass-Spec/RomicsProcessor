% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/07_Grouping.R
\name{romicsUmapPlot}
\alias{romicsUmapPlot}
\title{romicsUmapPlot()}
\usage{
romicsUmapPlot(
  romics_object,
  umap_config = umap.defaults,
  label = TRUE,
  method = c("naive", "umap-learn"),
  lock_seed = TRUE,
  seed = 42,
  color_by = c("factor", "variable"),
  factor = "main",
  variable = "variable_name",
  color_palette = viridis(n = 20),
  scale = TRUE
)
}
\arguments{
\item{romics_object}{has to be a log transformed romics_object created using romicsCreateObject() and transformed using the function log2transform() or log10transform()}

\item{umap_config}{object of class umap.config (see umap() function documentation)}

\item{method}{character, implementation. Available methods are 'naive' (an implementation written in pure R) and 'umap-learn' (requires python package 'umap-learn')}

\item{lock_seed}{has to be TRUE or FALSE to indicate if the seed is locked to enable reproducing plotting}

\item{seed}{numeric value indicating what seed to use when random seed are used}

\item{color_by}{has to be either 'factor' or 'variable to indicate if the coloring of the sample points is done using the factor set with the parameter 'factor' or with the numerical values of a variable feature set with the factor 'variable'}

\item{factor}{has to be a factor of the romics_object, the existing factors names can be found using the function romicsFactorNames()}

\item{color_palette}{has to be a R gradient color palette (viridis(n=20) by default) this gradient color palette will be used to color the points when the 'color_by' parameter is set to 'variable'}

\item{scale}{has to be TRUE or FALSE to indicate if the palette defined in color_palette has to be displayed or not.}

\item{...}{further arguments passed to or from other methods}
}
\value{
Returns either a ggplot2 showing the sample umap plot
}
\description{
Plots the samples on the umap. Can be colored either by a factor or by the values of a given variable present in the dataset. This function require the data to be full in the layer data to be run. Calculate the umap of the data layer of the romics_object using the package umap.
}
\author{
Geremy Clair
}
