% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/export_databases.R
\name{export4msfinder}
\alias{export4msfinder}
\title{Export Database for MS-FINDER}
\usage{
export4msfinder(data, output = "database for msfinder.txt")
}
\arguments{
\item{data}{The list of chemicals should contain at least Name, SMILES, CID,
InChIKey, Formula.}

\item{output}{The path where the structure database will be stored. It must be
in txt format. For example, "c:/data/database for msfinder.txt"}
}
\description{
\code{export4msfinder} provides a way to convert any list of chemicals into a
structure database that can be used by MS-FINDER. please refer to
\link{extract_cid} and \link{extract_meta} regarding meta data extraction.
}
