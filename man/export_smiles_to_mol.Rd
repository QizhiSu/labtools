% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chemical_structure_convertion.R
\name{export_smiles_to_mol}
\alias{export_smiles_to_mol}
\title{Export SMILES strings to MOL files}
\usage{
export_smiles_to_mol(
  df,
  id_col = "ID",
  smiles_col = "SMILES",
  output_dir = "mol_files"
)
}
\arguments{
\item{df}{A data.frame or tibble containing SMILES strings and compound identifiers.}

\item{id_col}{A string specifying the column name in \code{df} that contains unique IDs
for naming the output MOL files. Default is \code{"ID"}.}

\item{smiles_col}{A string specifying the column name in \code{df} that contains SMILES strings.
Default is \code{"SMILES"}.}

\item{output_dir}{A string specifying the directory to which MOL files will be written.
If the directory does not exist, it will be created recursively. Default is \code{"mol_files"}.}
}
\value{
This function is called for its side effect of writing MOL files to disk.
It returns a list invisibly containing processing summary with elements:
\code{success} (number of successfully processed molecules),
\code{failed} (number of failed conversions), and
\code{skipped} (number of skipped rows due to NA/NULL values).
}
\description{
This function converts SMILES strings in a data frame into individual MOL files,
which are commonly used for molecular visualization and structure-based computation.
Each row in the input data frame should contain a unique compound identifier and a valid SMILES string.
The function uses the \pkg{rcdk} backend for parsing and rendering molecular structures.
}
\details{
2D coordinates will be generated automatically for each molecule,
making the output MOL files suitable for display in chemical structure viewers.
If a SMILES string fails to parse, a warning will be issued, and that entry will be skipped.

The function uses the \pkg{rcdk} package's \code{parse.smiles()} and \code{write.molecules()}
functions to process molecular structures. The \code{generate.2d.coordinates()} function is
used to create 2D layout, which improves compatibility with structure-drawing tools.

Invalid or unparsable SMILES strings will be skipped, with a warning indicating the failed ID.
}
\examples{
\dontrun{
# Example usage
df <- data.frame(
  ID = c("Caffeine", "Aspirin"),
  SMILES = c("Cn1cnc2c1c(=O)n(c(=O)n2C)C", "CC(=O)OC1=CC=CC=C1C(=O)O"),
  stringsAsFactors = FALSE
)
export_smiles_to_mol(df, output_dir = "my_mol_files")
}

}
