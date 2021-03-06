% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{exportEmbCol}
\alias{exportEmbCol}
\title{a function to export embedding and cloumn to csv file}
\usage{
exportEmbCol(
  seurat = NULL,
  cells = NULL,
  exportCol = "seurat_clusters",
  embedding = "umap",
  outdir = "EmbCol"
)
}
\arguments{
\item{seurat}{a seurat object}

\item{cells}{cell barcode vector.}

\item{exportCol}{column names in metadata.}

\item{embedding}{reduction names in seurat object.}

\item{outdir}{path to save tables}
}
\description{
a function to export embedding and cloumn to csv file
}
