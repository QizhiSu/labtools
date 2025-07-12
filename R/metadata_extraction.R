#' Parse and clean CAS number strings
#'
#' This function processes strings that may contain multiple CAS numbers
#' separated by arbitrary characters and returns a cleaned character vector
#' of valid CAS numbers in the format "xxx-xx-x". Leading zeros are removed.
#'
#' The function handles various input formats including:
#' - Multiple CAS numbers separated by semicolons, commas, or other delimiters
#' - CAS numbers with leading zeros
#' - CAS numbers enclosed in brackets or other characters
#'
#' @param cas_string A character string containing one or more CAS numbers.
#'   Can be NA or empty string.
#'
#' @return A character vector of cleaned and validated CAS numbers in the format
#'   "xxx-xx-x" where x represents digits. If no valid CAS is found, returns
#'   \code{NA_character_} and emits a warning message.
#'
#' @details
#' The function uses regular expressions to:
#' \itemize{
#'   \item Remove all non-digit and non-hyphen characters
#'   \item Split the string into potential CAS number components
#'   \item Remove leading zeros from the first block of digits
#'   \item Validate against the standard CAS format (2-7 digits, 2 digits, 1 digit)
#' }
#'
#' @examples
#' # Multiple CAS numbers with various separators
#' parse_cas_clean("00064-17-5; [50-00-0]\n00110-54-3")
#' # Returns: c("64-17-5", "50-00-0", "110-54-3")
#'
#' # Invalid input
#' parse_cas_clean("abc;1234") # Returns NA with a message
#'
#' # Empty or NA input
#' parse_cas_clean(NA) # Returns NA
#' parse_cas_clean("") # Returns NA
#'
#' @seealso \code{\link{extract_cid}} for extracting PubChem CIDs using CAS numbers
#' @export
parse_cas_clean <- function(cas_string) {
  # Input validation
  if (is.na(cas_string) || cas_string == "" || is.null(cas_string)) {
    return(NA_character_)
  }

  # Ensure input is character
  cas_string <- as.character(cas_string)

  # Replace all non-digit and non-hyphen characters with space
  cas_string <- gsub("[^0-9\\-]+", " ", cas_string)

  # Split by whitespace and remove empty strings
  cas_parts <- unlist(strsplit(cas_string, "\\s+"))
  cas_parts <- cas_parts[nzchar(cas_parts)]

  # Return NA if no parts found
  if (length(cas_parts) == 0) {
    return(NA_character_)
  }

  # Remove leading zeros in the first block (if full format matches)
  cas_parts <- gsub("^0*(\\d{2,7}-\\d{2}-\\d{1})$", "\\1", cas_parts)

  # Keep only valid CAS number format (2-7 digits, 2 digits, 1 digit)
  cas_pattern <- "^\\d{2,7}-\\d{2}-\\d{1}$"
  cas_parts <- trimws(cas_parts)
  cas_valid <- unique(cas_parts[grepl(cas_pattern, cas_parts)])

  if (length(cas_valid) == 0) {
    message("WARNING: No valid CAS numbers found in: ", cas_string)
    return(NA_character_)
  }

  cas_valid
}

#' Extract PubChem CIDs via InChIKey, CAS, or Name (wrapper of webchem::get_cid)
#'
#' A convenience wrapper around \code{\link[webchem]{get_cid}} that looks up PubChem Compound IDs (CIDs)
#' using InChIKey, CAS number(s), or chemical names - in that order of priority. It supports flexible column
#' specifications (either by name or index), parses CAS strings containing multiple entries, and gives informative
#' messages if specified columns are missing.
#'
#' @param data A data.frame containing chemical information.
#' @param name_col Name or index of the column containing substance names (optional).
#'   Can be a character string (column name) or integer (column index). Default is "Name".
#' @param cas_col Name or index of the column containing CAS numbers (optional).
#'   Can be a character string (column name) or integer (column index). Default is "CAS".
#' @param inchikey_col Name or index of the column containing InChIKeys (optional).
#'   Can be a character string (column name) or integer (column index). Default is "InChIKey".
#' @param verbose Logical. Whether to print progress messages. Default is \code{TRUE}.
#'
#' @return A data.frame with a new or updated \code{CID} column containing matched PubChem CIDs.
#'   Existing CID values are preserved and only missing values are filled.
#'
#' @details
#' The function searches for CIDs in the following priority order:
#' \enumerate{
#'   \item InChIKey (most reliable)
#'   \item CAS numbers (parsed using \code{\link{parse_cas_clean}})
#'   \item Chemical names (least reliable due to naming variations)
#' }
#'
#' The function handles network errors gracefully and provides detailed progress information.
#' Multiple CAS numbers in a single cell are supported and parsed automatically.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage with name and CAS columns
#' df <- data.frame(
#'   Name = c("benzene", "acetic acid"),
#'   CAS = c("71-43-2", "64-19-7")
#' )
#' result <- extract_cid(df)
#'
#' # Using column indices instead of names
#' extract_cid(df, name_col = 1, cas_col = 2)
#'
#' # With InChIKey column
#' df_with_inchi <- data.frame(
#'   Name = "benzene",
#'   InChIKey = "UHOVQNZJYSORNB-UHFFFAOYSA-N"
#' )
#' extract_cid(df_with_inchi)
#' }
#'
#' @seealso \code{\link{parse_cas_clean}} for CAS number parsing,
#'   \code{\link[webchem]{get_cid}} for the underlying PubChem query function
extract_cid <- function(data,
                        name_col = "Name",
                        cas_col = "CAS",
                        inchikey_col = "InChIKey",
                        verbose = TRUE) {
  # Ensure CID column exists
  if (!"CID" %in% names(data)) data$CID <- NA_integer_

  # Helper: Convert column index to name
  col_to_name <- function(col, df) {
    if (is.null(col)) {
      return(NULL)
    }
    if (is.numeric(col)) names(df)[col] else as.character(col)
  }

  # Convert input column identifiers to actual column names
  name_col <- col_to_name(name_col, data)
  cas_col <- col_to_name(cas_col, data)
  inchikey_col <- col_to_name(inchikey_col, data)

  # Check column existence and provide informative messages
  all_cols <- list(name = name_col, cas = cas_col, inchikey = inchikey_col)
  present <- vapply(all_cols, function(x) !is.null(x) && x %in% names(data), logical(1))
  missing_cols <- unlist(all_cols[!present])

  if (all(!present)) {
    stop(
      "ERROR: None of the specified columns exist in the input data: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (any(!present)) {
    message(
      "WARNING: The following columns are not found and will be skipped: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  # Helper: Create lookup named vector from webchem result
  make_lookup <- function(res_df) {
    if (is.null(res_df) || nrow(res_df) == 0) {
      return(list())
    }

    # Check for expected columns in webchem result
    if (!"cid" %in% names(res_df)) {
      warning("Expected 'cid' column not found in webchem result")
      return(list())
    }

    # Use 'query' column if available, otherwise use row names or indices
    query_col <- if ("query" %in% names(res_df)) {
      res_df$query
    } else {
      rownames(res_df)
    }

    valid <- !is.na(res_df$cid) & res_df$cid != 0
    if (sum(valid) == 0) {
      return(list())
    }

    setNames(res_df$cid[valid], query_col[valid])
  }

  # --- Track results ---
  n_start <- sum(!is.na(data$CID))
  n_total <- nrow(data)
  n_found_inchikey <- 0
  n_found_cas <- 0
  n_found_name <- 0

  # Step 1: InChIKey
  if (present["inchikey"]) {
    idx <- which(is.na(data$CID))
    ikeys <- unique(na.omit(data[[inchikey_col]][idx]))
    ikeys <- ikeys[nzchar(ikeys)]
    if (length(ikeys) > 0) {
      if (verbose) message("Extracting CID from InChIKey...")
      res <- tryCatch(
        webchem::get_cid(ikeys, from = "inchikey", match = "first", verbose = verbose),
        error = function(e) NULL
      )
      lookup <- make_lookup(res)
      for (i in which(is.na(data$CID))) {
        ikey <- data[[inchikey_col]][i]
        if (!is.na(ikey) && ikey %in% names(lookup)) {
          data$CID[i] <- as.integer(lookup[[ikey]])
          n_found_inchikey <- n_found_inchikey + 1
        }
      }
    }
  }

  # Step 2: CAS
  if (present["cas"]) {
    idx <- which(is.na(data$CID))
    cas_lists <- lapply(data[[cas_col]][idx], parse_cas_clean)
    all_cas <- unique(unlist(cas_lists))
    all_cas <- all_cas[!is.na(all_cas)]
    if (length(all_cas) > 0) {
      if (verbose) message("Extracting CID from CAS...")
      res <- tryCatch(
        webchem::get_cid(all_cas, match = "first", verbose = verbose),
        error = function(e) NULL
      )
      lookup <- make_lookup(res)
      for (j in seq_along(idx)) {
        row <- idx[j]
        cas_vec <- cas_lists[[j]]
        if (length(cas_vec) == 0 || is.na(cas_vec[1])) next
        for (cas in cas_vec) {
          if (cas %in% names(lookup)) {
            data$CID[row] <- as.integer(lookup[[cas]])
            n_found_cas <- n_found_cas + 1
            break
          }
        }
      }
    }
  }

  # Step 3: Name
  if (present["name"]) {
    idx <- which(is.na(data$CID))
    name_vec <- unique(na.omit(data[[name_col]][idx]))
    name_vec <- name_vec[nzchar(name_vec)]
    if (length(name_vec) > 0) {
      if (verbose) message("Extracting CID from Name...")
      res <- tryCatch(
        webchem::get_cid(name_vec, from = "name", match = "first", verbose = verbose),
        error = function(e) NULL
      )
      lookup <- make_lookup(res)
      for (i in which(is.na(data$CID))) {
        n_val <- data[[name_col]][i]
        if (!is.na(n_val) && n_val %in% names(lookup)) {
          data$CID[i] <- as.integer(lookup[[n_val]])
          n_found_name <- n_found_name + 1
        }
      }
    }
  }

  # Final summary
  n_end <- sum(!is.na(data$CID))
  message("CID extraction summary:")
  if (present["inchikey"]) message("  From InChIKey: ", n_found_inchikey)
  if (present["cas"]) message("  From CAS: ", n_found_cas)
  if (present["name"]) message("  From Name: ", n_found_name)
  message(
    "Total found: ", n_end, " / ", n_total,
    " (newly found: ", n_end - n_start, ", not found: ", n_total - n_end, ")"
  )

  return(data)
}


#-------------------------------------------------------------------------------
#' Extract metadata from PubChem for compounds with known CIDs
#'
#' This function enriches a data frame with structural and contextual metadata from PubChem
#' using available CIDs, including molecular properties, CAS numbers, synonyms, Flavornet descriptors,
#' and uses information.
#'
#' @param data A data.frame containing a 'CID' column (PubChem compound IDs)
#' @param cas Logical. If TRUE, attempt to extract CAS numbers from pc_sect and synonyms.
#' @param flavornet Logical. If TRUE, attempt to extract Flavornet sensory descriptors.
#' @param synonyms Logical. If TRUE, attempt to extract synonym list from PubChem.
#' @param uses Logical. If TRUE, attempt to extract compound use information from PubChem.
#' @param verbose Logical. If TRUE, print progress messages.
#' @param checkpoint_dir Directory to save checkpoint files (optional, currently not used).
#'
#' @return A data.frame with new metadata columns added.
#' @importFrom stats aggregate
#' @export
extract_meta <- function(data,
                         cas = FALSE,
                         flavornet = FALSE,
                         synonyms = FALSE,
                         uses = FALSE,
                         verbose = TRUE,
                         checkpoint_dir = ".") {
  if (!"CID" %in% names(data)) stop("Missing CID column. Please run extract_cid() first.")
  data$CID <- as.integer(data$CID)
  cid_vector <- na.omit(unique(data$CID))

  # --- Structural metadata ---
  if (verbose) message("\nExtracting core PubChem metadata...")
  props <- c("MolecularFormula", "MolecularWeight", "SMILES", "InChI", "InChIKey", "IUPACName", "ExactMass")
  cid_info <- tryCatch(
    webchem::pc_prop(cid_vector, properties = props, verbose = verbose),
    error = function(e) NULL
  )

  if (!is.null(cid_info) && nrow(cid_info) > 0) {
    col_map <- c(
      MolecularFormula = "Formula", MolecularWeight = "MW", SMILES = "SMILES",
      InChI = "InChI", InChIKey = "InChIKey", IUPACName = "IUPAC_Name", ExactMass = "Exact_Mass"
    )
    for (old in names(col_map)) {
      if (old %in% names(cid_info)) names(cid_info)[names(cid_info) == old] <- col_map[[old]]
    }
    idx <- match(data$CID, cid_info$CID)
    for (col in names(cid_info)) {
      if (col != "CID") {
        new_col <- col
        if (col %in% names(data)) {
          suffix_idx <- 1
          while (paste0(col, "_", suffix_idx) %in% names(data)) {
            suffix_idx <- suffix_idx + 1
          }
          new_col <- paste0(col, "_", suffix_idx)
          if (verbose) message(sprintf("!! Column '%s' already exists, writing to '%s' instead.", col, new_col))
        }
        data[[new_col]] <- cid_info[[col]][idx]
      }
    }
  }

  # --- CAS numbers ---
  if (cas) {
    if (verbose) message("\nExtracting CAS numbers from pc_sect and synonyms...")
    data$CAS <- NA_character_
    cas_sect <- tryCatch(
      webchem::pc_sect(cid_vector, section = "cas", verbose = verbose),
      error = function(e) NULL
    )
    if (!is.null(cas_sect) && "CID" %in% names(cas_sect)) {
      cas_map <- setNames(cas_sect$Result, cas_sect$CID)
      idx <- match(data$CID, names(cas_map))
      data$CAS[!is.na(idx)] <- cas_map[as.character(data$CID[!is.na(idx)])]
      if (verbose) message(sprintf("CAS found via pc_sect: %d", sum(!is.na(data$CAS))))
    }

    # Try synonyms if still missing
    missing_cas_idx <- which(is.na(data$CAS) & !is.na(data$CID))
    if (length(missing_cas_idx) > 0) {
      syn_df <- tryCatch(webchem::pc_synonyms(data$CID[missing_cas_idx], from = "cid", verbose = verbose),
        error = function(e) NULL
      )
      if (is.data.frame(syn_df) && "CID" %in% names(syn_df)) {
        for (i in seq_len(nrow(syn_df))) {
          cid <- syn_df$CID[i]
          syns <- unlist(syn_df$Synonym[i])
          cas_like <- grep("^\\d{2,7}-\\d{2}-\\d$", syns, value = TRUE)
          if (length(cas_like) > 0) {
            row_i <- which(data$CID == cid & is.na(data$CAS))
            data$CAS[row_i] <- cas_like[1]
            if (verbose) message(sprintf("CAS found via synonym for CID %s: %s", cid, cas_like[1]))
          }
        }
      }
    }
  }

  # --- Flavornet ---
  if (flavornet) {
    if (verbose) message("\nExtracting Flavornet sensory descriptors...")
    data$Flavornet <- NA_character_
    cas_vals <- na.omit(data$CAS)
    if (length(cas_vals) > 0) {
      for (i in seq_len(nrow(data))) {
        cas_i <- data$CAS[i]
        if (!is.na(cas_i)) {
          val <- tryCatch(webchem::fn_percept(cas_i), error = function(e) NA_character_)
          if (!is.na(val)) {
            data$Flavornet[i] <- val
            if (verbose) message(sprintf("Flavornet found for CAS %s: %s", cas_i, val))
          } else {
            if (verbose) message(sprintf("No Flavornet data for CAS %s", cas_i))
          }
        }
      }
    }
  }

  # --- Synonyms ---
  if (synonyms) {
    if (verbose) message("\nExtracting synonyms from PubChem...")

    cid_vector <- na.omit(unique(data$CID))
    syn_list <- tryCatch(
      webchem::pc_synonyms(cid_vector, from = "cid", verbose = verbose),
      error = function(e) {
        if (verbose) message("Error in pc_synonyms: ", e$message)
        NULL
      }
    )

    if (is.list(syn_list) && length(syn_list) > 0) {
      if (verbose) message(sprintf("Retrieved synonyms for %d compounds", length(syn_list)))

      # Collapse each CID's synonyms into a string
      syn_map <- vapply(
        syn_list,
        function(syns) {
          if (length(syns) > 0 && any(nzchar(syns))) {
            paste(unique(syns[!is.na(syns) & nzchar(syns)]), collapse = "; ")
          } else {
            NA_character_
          }
        },
        character(1)
      )

      # Create column if not already
      if (!"Synonyms" %in% names(data)) data$Synonyms <- NA_character_

      # Assign to data
      assigned <- 0
      for (i in seq_len(nrow(data))) {
        cid <- as.character(data$CID[i])
        if (!is.na(cid) && cid %in% names(syn_map)) {
          data$Synonyms[i] <- syn_map[[cid]]
          assigned <- assigned + 1
        }
      }

      if (verbose) message(sprintf("Synonyms assigned to %d rows", assigned))
    } else {
      if (verbose) message("WARNING: No synonym data retrieved or invalid format")
    }
  }

  # --- Uses ---
  if (uses) {
    if (verbose) message("\nExtracting 'Uses' information from PubChem...")
    uses_df <- tryCatch(webchem::pc_sect(cid_vector, section = "uses", verbose = verbose),
      error = function(e) NULL
    )
    data$Uses <- NA_character_
    if (!is.null(uses_df) && "CID" %in% names(uses_df)) {
      uses_agg <- aggregate(Result ~ CID, data = uses_df, FUN = function(x) paste(unique(x), collapse = "; "))
      uses_map <- setNames(uses_agg$Result, uses_agg$CID)
      idx <- match(data$CID, names(uses_map))
      data$Uses[!is.na(idx)] <- uses_map[as.character(data$CID[!is.na(idx)])]
      if (verbose) {
        message(sprintf(
          "Uses extraction done: %d found, %d not found.",
          sum(!is.na(data$Uses)), nrow(data) - sum(!is.na(data$Uses))
        ))
      }
    } else {
      if (verbose) message("No 'Uses' data found for any CID.")
    }
  }

  if (requireNamespace("beepr", quietly = TRUE)) beepr::beep(1)
  return(data)
}


#-------------------------------------------------------------------------------
#' Extract chemical classification via ClassyFire
#'
#' This function queries the [ClassyFire](http://classyfire.wishartlab.com/) database
#' using InChIKey for each unique compound in a dataset and retrieves hierarchical
#' chemical classifications (e.g., Kingdom, Superclass, Class, Subclass).
#'
#' It displays progress and messages for each compound during processing. If a classification
#' is found, it is merged back to the original dataset by InChIKey. Classification data
#' is optionally repositioned after known reference columns (e.g., "ExactMass").
#'
#' @param data A data.frame or tibble containing at least one column with InChIKeys.
#' @param inchikey_col A string specifying the name of the column that contains InChIKeys. Default is `"InChIKey"`.
#' @param name_col A string specifying the name of the column that contains compound names (for message display). Default is `"Name"`.
#'
#' @return A data.frame with appended classification columns such as Kingdom, Superclass, Class, etc.
#'
#' @details
#' This function only uses InChIKey to retrieve data via `classyfireR::get_classification()`.
#' The fallback using SMILES/InChI via `submit_query()` is **not used**, as that API is currently unavailable.
#'
#' @examples
#' \dontrun{
#' library(classyfireR)
#' df <- data.frame(
#'   InChIKey = c("XMGQYMWWDOXHJM-UHFFFAOYSA-N"),
#'   Name = c("Limonene")
#' )
#' extract_classyfire(df)
#' }
#'
#' @importFrom dplyr left_join relocate all_of bind_rows
#' @export
#'
extract_classyfire <- function(data, inchikey_col = "InChIKey", name_col = "Name") {
  if (!requireNamespace("classyfireR", quietly = TRUE)) {
    stop("Please install the 'classyfireR' package.")
  }

  if (!(inchikey_col %in% colnames(data))) {
    stop("InChIKey column not found in the data.")
  }

  # Internal function to extract classification from classification() result
  extract_cla <- function(x) {
    if (is.null(x) || !is.data.frame(x) || nrow(x) == 0) {
      return(NULL)
    }

    # Create a data.frame with classification levels as column names
    # and classification values as the single row
    if ("Classification" %in% names(x) && "Level" %in% names(x)) {
      cla <- data.frame(t(x$Classification), stringsAsFactors = FALSE)
      colnames(cla) <- x$Level
      cla
    } else {
      NULL
    }
  }

  inchikeys <- unique(data[[inchikey_col]])
  classy_list <- list()

  for (i in seq_along(inchikeys)) {
    inchikey <- inchikeys[i]
    message(sprintf("%d/%d", i, length(inchikeys)))

    compound_name <- if (name_col %in% colnames(data)) {
      unique(data[data[[inchikey_col]] == inchikey, name_col])[1]
    } else {
      inchikey
    }

    result <- suppressMessages(tryCatch(
      {
        classyfireR::get_classification(inchikey)
      },
      error = function(e) NULL
    ))

    if (!is.null(result)) {
      class_info <- tryCatch(
        {
          classyfireR::classification(result)
        },
        error = function(e) NULL
      )

      class_df <- tryCatch(
        {
          extract_cla(class_info)
        },
        error = function(e) NULL
      )

      if (!is.null(class_df)) {
        message(sprintf("Classification found via InChIKey for %s", compound_name))
        class_df[[inchikey_col]] <- inchikey
        classy_list[[inchikey]] <- class_df
      } else {
        message(sprintf("XX Unclassified for %s", compound_name))
      }
    } else {
      message(sprintf("XX Unclassified for %s", compound_name))
    }
  }

  if (length(classy_list) == 0) {
    warning("No classification found for any compound.")
    return(data)
  }

  # bind results
  classy_df <- dplyr::bind_rows(classy_list)
  data <- dplyr::left_join(data, classy_df, by = inchikey_col)

  # repositioning classification columns if needed
  known_cols <- c("Flavornet", "CAS_retrieved", "ExactMass")
  target_cols <- colnames(classy_df)[!colnames(classy_df) %in% inchikey_col]

  for (ref in known_cols) {
    if (ref %in% colnames(data)) {
      data <- dplyr::relocate(data, dplyr::all_of(target_cols), .after = ref)
      break
    }
  }

  if (requireNamespace("beepr", quietly = TRUE)) beepr::beep(1)
  return(data)
}


#-------------------------------------------------------------------------------
#' Assign structural metadata (SMILES, InChIKey, InChI) to your dataset
#'
#' If some compounds are not found in PubChem, \code{extract_meta()} will not
#' retrieve their structural information such as SMILES, InChIKey, or InChI.
#' This function provides a way to manually assign these values based on name matching.
#'
#' A metadata file (e.g., CSV or TXT) containing at least compound names and SMILES is required.
#' Column names in the file are case-insensitive and must include "Name" and "SMILES".
#' If available, "InChIKey" and "InChI" will also be assigned.
#'
#' There are two ways to prepare this metadata file:
#' \enumerate{
#'   \item Manually create a file with at least "Name" and "SMILES" columns.
#'   \item Convert MOL files into a structured table using
#'         \code{combine_mol2sdf()} and \code{extract_structure()} from the
#'         \pkg{mspcompiler} package:
#'         \url{https://github.com/QizhiSu/mspcompiler}.
#' }
#'
#' @param data A data.frame or tibble, typically the result after running \code{extract_meta()}.
#'   Must contain a "Name" column for matching.
#' @param meta_file A text-based file (e.g., *.csv or *.txt) containing at least "Name" and "SMILES",
#'   and optionally "InChIKey" and "InChI". File can be created manually or generated from MOL files.
#'   Supported formats include CSV, TSV, Excel, and other formats supported by \code{\link[rio]{import}}.
#'
#' @return A data.frame or tibble with missing SMILES, InChIKey, and InChI values filled in where available.
#'   Only missing values are updated; existing values are preserved.
#'
#' @details
#' The function performs case-insensitive, whitespace-trimmed matching between compound names
#' in the main dataset and the metadata file. This ensures robust matching even with minor
#' formatting differences.
#'
#' Column name matching is also case-insensitive for flexibility in input file formats.
#'
#' @note
#' Compound names in your metadata file must match those in your main dataset for successful
#' assignment. Consider standardizing compound names before using this function.
#'
#' @examples
#' \dontrun{
#' # Assuming you have a CSV file with compound metadata
#' data_with_meta <- assign_meta(my_data, "compound_metadata.csv")
#'
#' # The metadata file should have columns like:
#' # Name, SMILES, InChIKey, InChI
#' # "benzene", "c1ccccc1", "UHOVQNZJYSORNB-UHFFFAOYSA-N", "InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H"
#' }
#'
#' @importFrom rio import
#' @export
#'
assign_meta <- function(data, meta_file) {
  # Input validation
  if (!file.exists(meta_file)) {
    stop("Metadata file not found: ", meta_file)
  }

  if (!"Name" %in% names(data)) {
    stop("Input data must contain a 'Name' column for matching")
  }

  # Import metadata file
  meta_data <- tryCatch(
    {
      rio::import(meta_file)
    },
    error = function(e) {
      stop("Failed to read metadata file: ", e$message)
    }
  )

  # Validate metadata file structure
  if (!"Name" %in% names(meta_data)) {
    # Try to find name column with case-insensitive matching
    name_cols <- grep("name", names(meta_data), ignore.case = TRUE, value = TRUE)
    if (length(name_cols) == 0) {
      stop("Metadata file must contain a 'Name' column")
    }
    names(meta_data)[names(meta_data) == name_cols[1]] <- "Name"
  }

  # Standardize column names (case-insensitive matching)
  col_mapping <- list(
    "smiles" = "SMILES",
    "inchikey" = "InChIKey",
    "^inchi$" = "InChI"
  )

  for (pattern in names(col_mapping)) {
    target_col <- col_mapping[[pattern]]
    matching_cols <- grep(pattern, names(meta_data), ignore.case = TRUE, value = TRUE)
    if (length(matching_cols) > 0) {
      names(meta_data)[names(meta_data) == matching_cols[1]] <- target_col
    }
  }

  # Check for required SMILES column
  if (!"SMILES" %in% names(meta_data)) {
    stop("Metadata file must contain a 'SMILES' column")
  }

  # Standardize matching keys (case-insensitive, whitespace-trimmed)
  key_data <- trimws(tolower(data$Name))
  key_meta <- trimws(tolower(meta_data$Name))
  idx <- match(key_data, key_meta)

  # Count matches for reporting
  n_matches <- sum(!is.na(idx))
  message(sprintf(
    "Found %d matches out of %d compounds (%.1f%%)",
    n_matches, nrow(data), 100 * n_matches / nrow(data)
  ))

  # Assign SMILES
  if (!"SMILES" %in% names(data)) data$SMILES <- NA_character_
  is_missing_smiles <- is.na(data$SMILES)
  n_smiles_before <- sum(!is.na(data$SMILES))
  data$SMILES[is_missing_smiles] <- meta_data$SMILES[idx[is_missing_smiles]]
  n_smiles_after <- sum(!is.na(data$SMILES))
  message(sprintf(
    "SMILES: %d assigned (%d -> %d)",
    n_smiles_after - n_smiles_before, n_smiles_before, n_smiles_after
  ))

  # Assign InChIKey if available
  if ("InChIKey" %in% names(meta_data)) {
    if (!"InChIKey" %in% names(data)) data$InChIKey <- NA_character_
    is_missing_inchikey <- is.na(data$InChIKey)
    n_inchikey_before <- sum(!is.na(data$InChIKey))
    data$InChIKey[is_missing_inchikey] <- meta_data$InChIKey[idx[is_missing_inchikey]]
    n_inchikey_after <- sum(!is.na(data$InChIKey))
    message(sprintf(
      "InChIKey: %d assigned (%d -> %d)",
      n_inchikey_after - n_inchikey_before, n_inchikey_before, n_inchikey_after
    ))
  }

  # Assign InChI if available
  if ("InChI" %in% names(meta_data)) {
    if (!"InChI" %in% names(data)) data$InChI <- NA_character_
    is_missing_inchi <- is.na(data$InChI)
    n_inchi_before <- sum(!is.na(data$InChI))
    data$InChI[is_missing_inchi] <- meta_data$InChI[idx[is_missing_inchi]]
    n_inchi_after <- sum(!is.na(data$InChI))
    message(sprintf(
      "InChI: %d assigned (%d -> %d)",
      n_inchi_after - n_inchi_before, n_inchi_before, n_inchi_after
    ))
  }

  data
}
