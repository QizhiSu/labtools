% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_spectrum.R
\name{cosine_similarity}
\alias{cosine_similarity}
\title{Calculate cosine similarity between two spectra}
\usage{
cosine_similarity(spec1, spec2)
}
\arguments{
\item{spec1}{A named numeric vector representing spectrum 1 (m/z as names, intensity as values).}

\item{spec2}{A named numeric vector representing spectrum 2.}
}
\value{
A numeric value between 0 and 1 representing the cosine similarity.
If there are no common m/z values, returns 0.
}
\description{
Computes the cosine similarity between two numeric spectra, typically mass spectra,
by matching shared m/z values and evaluating vector similarity.
}
\examples{
spec1 <- update_spectrum("100:10 150:20 200:30")
spec2 <- update_spectrum("100:10 150:25 250:40")
cosine_similarity(spec1, spec2)

}
