#' Export SMILES strings to MOL files
#'
#' This function converts SMILES strings in a data frame into individual MOL files,
#' which are commonly used for molecular visualization and structure-based computation.
#' Each row in the input data frame should contain a unique compound identifier and a valid SMILES string.
#' The function uses the \pkg{rcdk} backend for parsing and rendering molecular structures.
#'
#' 2D coordinates will be generated automatically for each molecule,
#' making the output MOL files suitable for display in chemical structure viewers.
#' If a SMILES string fails to parse, a warning will be issued, and that entry will be skipped.
#'
#' @param df A data.frame or tibble containing SMILES strings and compound identifiers.
#' @param id_col A string specifying the column name in \code{df} that contains unique IDs
#'        for naming the output MOL files. Default is \code{"ID"}.
#' @param smiles_col A string specifying the column name in \code{df} that contains SMILES strings.
#'        Default is \code{"SMILES"}.
#' @param output_dir A string specifying the directory to which MOL files will be written.
#'        If the directory does not exist, it will be created recursively. Default is \code{"mol_files"}.
#'
#' @return This function is called for its side effect of writing MOL files to disk.
#'         It returns a list invisibly containing processing summary with elements:
#'         \code{success} (number of successfully processed molecules),
#'         \code{failed} (number of failed conversions), and
#'         \code{skipped} (number of skipped rows due to NA/NULL values).
#'
#' @details
#' The function uses the \pkg{rcdk} package's \code{parse.smiles()} and \code{write.molecules()}
#' functions to process molecular structures. The \code{generate.2d.coordinates()} function is
#' used to create 2D layout, which improves compatibility with structure-drawing tools.
#'
#' Invalid or unparsable SMILES strings will be skipped, with a warning indicating the failed ID.
#'
#' @importFrom rcdk parse.smiles generate.2d.coordinates write.molecules
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage
#' df <- data.frame(
#'   ID = c("Caffeine", "Aspirin"),
#'   SMILES = c("Cn1cnc2c1c(=O)n(c(=O)n2C)C", "CC(=O)OC1=CC=CC=C1C(=O)O"),
#'   stringsAsFactors = FALSE
#' )
#' export_smiles_to_mol(df, output_dir = "my_mol_files")
#' }
#'
export_smiles_to_mol <- function(df, id_col = "ID", smiles_col = "SMILES", output_dir = "mol_files") {
  # Input validation
  if (!is.data.frame(df)) {
    stop("df must be a data.frame or tibble")
  }

  if (nrow(df) == 0) {
    warning("Input data.frame is empty")
    return(invisible(list(success = 0, failed = 0, skipped = 0)))
  }

  if (!id_col %in% names(df)) {
    stop(sprintf("Column '%s' not found in data.frame", id_col))
  }

  if (!smiles_col %in% names(df)) {
    stop(sprintf("Column '%s' not found in data.frame", smiles_col))
  }

  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Initialize counters for summary
  success_count <- 0
  failed_count <- 0
  skipped_count <- 0

  # Loop through each row
  for (i in seq_len(nrow(df))) {
    id <- df[[id_col]][i]
    smiles <- df[[smiles_col]][i]

    # Skip rows with NA or NULL values
    if (is.na(id) || is.null(id) || is.na(smiles) || is.null(smiles)) {
      warning(sprintf("Skipping row %d: NA or NULL values found", i))
      skipped_count <- skipped_count + 1
      next
    }

    # Parse SMILES to molecule
    mol <- tryCatch(parse.smiles(smiles)[[1]], error = function(e) NULL)

    if (!is.null(mol)) {
      # Generate 2D coordinates (optional but useful for visualization)
      mol <- tryCatch(generate.2d.coordinates(mol), error = function(e) NULL)

      # Check if coordinate generation was successful
      if (!is.null(mol)) {
        # Write to MOL file
        file_path <- file.path(output_dir, paste0(id, ".mol"))
        tryCatch({
          write.molecules(mol, file_path)
          success_count <- success_count + 1
        }, error = function(e) {
          warning(sprintf("Failed to write MOL file for ID: %s. Error: %s", id, e$message))
          failed_count <- failed_count + 1
        })
      } else {
        warning(sprintf("Failed to generate 2D coordinates for ID: %s", id))
        failed_count <- failed_count + 1
      }
    } else {
      warning(sprintf("Failed to parse SMILES: %s (ID: %s)", smiles, id))
      failed_count <- failed_count + 1
    }
  }

  # Print summary
  message(sprintf("Processing complete: %d successful, %d failed, %d skipped",
                  success_count, failed_count, skipped_count))

  # Return summary invisibly
  invisible(list(success = success_count, failed = failed_count, skipped = skipped_count))
}
