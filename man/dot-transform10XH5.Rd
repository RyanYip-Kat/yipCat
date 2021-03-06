% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cellranger.R
\name{.transform10XH5}
\alias{.transform10XH5}
\title{transform cutsom matrix into 10X h5 matrix}
\usage{
.transform10XH5(
  crgPython = "python",
  targetMatrix = "outs/filtered_feature_bc_matrix",
  CRH5 = "filtered_feature_bc_matrix.h5",
  outH5 = "filtered_feature_bc_matrix_gold.h5",
  path2CRG = "cellranger"
)
}
\arguments{
\item{crgPython}{path of  cellranger python,eg(cellranger-3.1.0/miniconda-cr-cs/4.3.21-miniconda-cr-cs-c10/bin/python)}

\item{targetMatrix}{target matrix want to transform eg(outs/filtered_feature_bc_matrix)}

\item{outH5}{h5 file to store output}

\item{path2CRG}{path of cellranger}
}
\description{
transform cutsom matrix into 10X h5 matrix
}
