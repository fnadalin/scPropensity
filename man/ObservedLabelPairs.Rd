% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scPropensity.R
\name{ObservedLabelPairs}
\alias{ObservedLabelPairs}
\title{ObservedLabelPairs}
\usage{
ObservedLabelPairs(object, slot, label.list = NULL)
}
\arguments{
\item{object}{The input Seurat object}

\item{slot}{The slot in object@meta.data containing the labels assigned to the cells}

\item{label.list}{external list of labels to consider (default: all labels in the metadata slot)}
}
\value{
The observed matrix
}
\description{
Compute the matrix containing the observed count of paired labels
}
