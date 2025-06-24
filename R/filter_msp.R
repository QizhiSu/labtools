#' Filter an EI library by InChIKey
#'
#' \code{filter_msp} provide a simple way to filter an EI library by InChIKey, which can then be use
#' as a small *.msp file for MS-DIAL. This can be useful to targetedly find specific compounds in
#' your data.
#' @param msp Path to the library file in msp format
#' @param cmp_list A data.frame contains at least Name and InChIKey to be kept
#' @param keep_napd8 If TRUE, add Naphthalene-D8 to the list of compounds to keep
#' @param output path to write the filtered library to
#'
#' @import purrr
#' @import dplyr
#' @import mspcompiler
#'
#' @export
#'
filter_msp <- function(msp, cmp_list, keep_napd8 = TRUE, output) {
  # check if mspcompiler is installed
  if (!requireNamespace("mspcompiler", quietly = TRUE)) {
    stop("mspcompiler is not installed. Please install it first.",
      call. = FALSE
    )
  }

  # check if file exists
  if (!file.exists(msp)) {
    stop(msp, " file does not exist.", call. = FALSE)
  }

  if (!file.exists(dirname(output))) {
    stop(dirname(output), " path does not exist.", call. = FALSE)
  }

  # read in the library in msp format
  cat("Reading", basename(msp), "...\n")
  lib <- mspcompiler::read_lib(msp, type = "EI", remove_ri = FALSE)
  cat(basename(msp), "read!\n")

  # define the InChIKeys to keep
  if (keep_napd8) {
    cmp_list <- dplyr::bind_rows(
      cmp_list, data.frame(Name = "Naphthalene-D8", InChIKey = "UFWIBTONFRDIAS-LUPYIRNKSA-N")
    )
  }
  inchikey <- cmp_list$InChIKey

  # Define a function to check for the presence of a specific InChIKey
  keep_cmp <- function(element, inchikey) {
    !is.null(element$InChIKey) &&
      length(element$InChIKey) > 0 &&
      element$InChIKey %in% inchikey
  }

  # Filter the 'msp' list based on the presence of the desired InChIKey
  filtered_lib <- purrr::keep(lib, ~ keep_cmp(., inchikey))

  # determine if all InChIKeys were found
  filtered_inchikey <- sapply(filtered_lib, function(x) x$InChIKey)
  filtered_inchikey <- unique(filtered_inchikey)
  if (all(inchikey %in% filtered_inchikey)) {
    cat("All InChIKeys found!\n")
  } else {
    missing_inchikey <- inchikey[!inchikey %in% filtered_inchikey]
    which_missing <- which(inchikey %in% missing_inchikey)
    cat(
      "Warning: The following compounds were not found:\n",
      paste(which_missing, collapse = ", ")
    )
  }

  # change the Name to the one provided in the cmp_list
  change_name <- function(lib, cmp_list) {
    lib$Name <- cmp_list$Name[match(lib$InChIKey, cmp_list$InChIKey)]
    return(lib)
  }
  filtered_lib <- lapply(filtered_lib, change_name, cmp_list = cmp_list)

  # write the filtered library to file
  mspcompiler::write_EI_msp(filtered_lib, output)
  cat("\n", "Filtered library written to", output, "!\n")
}