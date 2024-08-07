% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/semiquantification.R
\name{calculate_con}
\alias{calculate_con}
\title{Semi-quantification of Analytes}
\usage{
calculate_con(df, sam_weight, start_col = 12)
}
\arguments{
\item{df}{A data frame containing peak areas and reference standards for each analyte.}

\item{sam_weight}{A data frame containing sample codes, weights, and final volumns.}

\item{start_col}{The starting column index of peak area columns (default is 12).}
}
\value{
A data frame with updated concentrations.
}
\description{
This function performs semi-quantification of analytes in a given data frame considering signal
correction by Naphthalene-D8. It calculates concentrations for each batch separately. For
negative concentrations, it will then use the lowest concentration of the reference standard to
do single-point correction. If sample weights are provided, the concentrations will be expressed
in samples. Otherwise, the concentrations will be basedo on instrument-measured signals only.
}
