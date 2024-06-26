% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scPropensity.R
\name{scPropensity}
\alias{scPropensity}
\title{scPropensity}
\usage{
scPropensity(object.file, slot, outdir, label.file = "", ...)
}
\arguments{
\item{slot}{The slot in object@meta.data containing the labels assigned to the cells}

\item{outdir}{Folder where to print the propensity matrices}

\item{object}{The input Seurat object}

\item{label.list}{external list of labels to consider (default: all labels in the metadata slot)}
}
\value{
NA
}
\description{
Wrapper function for pair propensity matrix calculation
}
