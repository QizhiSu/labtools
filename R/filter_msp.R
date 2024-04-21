#' Filter an EI library by InChIKey
#'
#' \code{filter_msp} provide a simple way to filter an EI library by InChIKey, which can then be use
#' as a small *.msp file for MS-DIAL. This can be useful to targetedly find specific compounds in
#' your data.
#' @param msp Path to the library file in msp format
#' @param inchikey vector of InChIKeys to keep
#' @param output path to write the filtered library to
#'
#' @import purrr
#' @import mspcompiler
#'
#' @export
#'
filter_msp <- function(msp, inchikey, output) {
  # check if mspcompiler is installed
  if (!requireNamespace("mspcompiler", quietly = TRUE)) {
    stop("mspcompiler is not installed. Please install it first.",
         call. = FALSE)
  }

  # read in the library in msp format
  cat("Reading", basename(msp), "...\n")
  lib <- mspcompiler::read_lib(msp, type = "EI")
  cat(basename(msp), "read!\n")

  # Define a function to check for the presence of a specific InChIKey
  keep_cmp <- function(element, inchikey) { # nolint
    if (is.null(element$InChIKey)) {
      return(FALSE) # Return false for NULL InChIKey
    } else if (length(element$InChIKey) == 0) {
      return(FALSE) # Return false for empty list InChIKey
    } else if (element$InChIKey %in% inchikey) {
      return(TRUE) # Return true for an exact match
    } else {
      return(FALSE)
    }
  }

  # Filter the 'msp' list based on the presence of the desired InChIKey
  filtered_lib <- purrr::keep(lib, ~ keep_cmp(., inchikey))

  # determine if all InChIKeys were found
  filtered_inchikeys <- sapply(filtered_lib, function(x) x$InChIKey)
  filtered_inchikeys <- unique(filtered_inchikeys)
  if (all(inchikey %in% filtered_inchikeys)) {
    cat("All InChIKeys found!\n")
  } else {
    missing_inchikeys <- inchikey[!inchikey %in% filtered_inchikeys]
    which_missing <- which(inchikey %in% missing_inchikeys)
    cat(
      "Warning: The following compounds were not found:\n",
      paste(which_missing, collapse = "\n")
    )
  }

  # write the filtered library to file
  mspcompiler::write_EI_msp(filtered_lib, output)
  cat("Filtered library written to", output, "!\n")
}