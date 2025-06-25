#' Extract PubChem Compound ID (CID) using various identifiers
#'
#' \code{extract_cid()} retrieves PubChem CIDs based on compound identifiers
#' such as InChIKey, CAS number, or compound name. It attempts each identifier
#' in order and records the CID in a new column. Lookup results are saved to a
#' checkpoint file to support resuming interrupted processes.
#'
#' The function uses \code{webchem::get_cid()} under the hood and handles
#' timeouts or web request errors gracefully. If a compound cannot be matched
#' to a CID, it is reported clearly in the console output. For efficiency and
#' fault tolerance, results are saved to disk every 10 records using
#' \code{saveRDS()}. If a matching checkpoint is found, the function resumes
#' from it.
#'
#' @param data A \code{data.frame} or \code{tibble} containing compound
#' information. A column named \code{"CID"} will be created if not already present.
#' @param name_col The column name or index containing compound names to be used for CID lookup. Defaults to `"Name"``.
#' @param cas_col The column name or index containing CAS numbers to be used for CID lookup. Defaults to `"CAS"``.
#' @param inchikey_col The column name or index containing InChIKeys to be used for CID lookup. Defaults to `"InChIKey"`.
#' @param timeout Timeout (in seconds) for each PubChem request. Defaults to 120.
#' @param verbose Logical. If \code{TRUE}, print progress and status messages. Defaults to \code{TRUE}.
#' @param checkpoint_file File path for storing and resuming from checkpoint. Defaults to \code{"cid_checkpoint.rds"}.
#' @param use_checkpoint Logical. Whether to use and save checkpoint file. Defaults to \code{TRUE}.
#'
#' @return The input \code{data.frame} with a \code{CID} column filled in. Any found CIDs
#' will be shown with a success message in the console indicating the identifier type used.
#'
#' @details
#' CID lookup is attempted in the following order (if columns are provided): InChIKey, CAS, Name.
#' Timeouts and unexpected errors during PubChem queries are caught silently and do not interrupt the process.
#' The function prints a summary when finished and optionally plays a notification sound if the \code{beepr} package is installed.
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   Name = c("bisphenol A", "acetone"),
#'   InChIKey = c("JUZNXVGZRYFYFW-UHFFFAOYSA-N", NA),
#'   CAS = c("80-05-7", NA),
#'   stringsAsFactors = FALSE
#' )
#' df_with_cid <- extract_cid(df, name_col = "Name", cas_col = "CAS", inchikey_col = "InChIKey")
#' }
#'
#' @seealso \code{\link[webchem]{get_cid}}, \code{\link[beepr]{beep}}
#'
#' @importFrom webchem get_cid
#' @importFrom R.utils withTimeout
#' @importFrom digest digest
#' @importFrom beepr beep
#'
#' @export

extract_cid <- function(data,
                        name_col = "Name",
                        cas_col = "CAS",
                        inchikey_col = "InChIKey",
                        timeout = 180,
                        verbose = TRUE,
                        checkpoint_file = "cid_checkpoint.rds",
                        use_checkpoint = TRUE) {

  # 内部函数：单个值查找 CID
  lookup_cid <- function(identifier, from, timeout, verbose, max_retries = 3, wait = 2) {
    if (is.na(identifier) || is.null(identifier) || identifier == "") {
      return(NA_character_)
    }

    for (attempt in seq_len(max_retries)) {
      result <- tryCatch(
        {
          R.utils::withTimeout({
            res <- webchem::get_cid(identifier, match = "first", from = from)
            cid <- res$cid
            if (is.null(cid) || length(cid) == 0 || cid == 0) {
              return(NA_character_)  # 查不到 => 不重试
            } else {
              return(as.character(cid))
            }
          }, timeout = timeout)
        },
        TimeoutException = function(e) NA_character_,
        error = function(e) NA_character_
      )

      # 成功或查不到 => 不重试
      if (!is.na(result)) return(result)

      # 否则等待重试
      if (attempt < max_retries) Sys.sleep(wait)
    }

    return(NA_character_) # 所有尝试都失败
  }

  # 初始化 CID 列
  if (!"CID" %in% colnames(data)) {
    data$CID <- NA_character_
  }

  # checkpoint 比对用 hash
  data_hash <- digest::digest(data)

  # 尝试从 checkpoint 恢复
  if (use_checkpoint && file.exists(checkpoint_file)) {
    checkpoint <- readRDS(checkpoint_file)
    if (!is.null(attr(checkpoint, "data_hash")) && attr(checkpoint, "data_hash") == data_hash) {
      if (verbose) message("✅ Resuming from checkpoint: ", checkpoint_file)
      data <- checkpoint
    } else if (verbose) {
      message("⚠️ Checkpoint file exists but doesn't match current data. Starting fresh.")
    }
  }

  id_cols <- list(
    inchikey = inchikey_col,
    cas = cas_col,
    name = name_col
  )

  n <- nrow(data)
  for (i in seq_len(n)) {
    if (!is.na(data$CID[i])) next
    if (verbose) message(i, "/", n)

    found_cid <- NA_character_
    used_method <- NULL

    for (id_type in names(id_cols)) {
      col <- id_cols[[id_type]]
      if (!is.null(col)) {
        id_value <- data[[col]][i]
        cid <- lookup_cid(id_value, from = id_type, timeout = timeout, verbose = verbose)
        if (!is.na(cid)) {
          found_cid <- cid
          used_method <- id_type
          break
        }
      }
    }

    data$CID[i] <- found_cid

    compound_name <- if (!is.null(name_col)) data[[name_col]][i] else paste0("row ", i)

    if (!is.na(found_cid)) {
      if (verbose) message("✅ CID ", found_cid, " found via ", toupper(used_method), " for: \"", compound_name, "\"")
    } else {
      if (verbose) message("❌ No CID found for: \"", compound_name, "\"")
    }

    if (i %% 10 == 0) {
      attr(data, "data_hash") <- data_hash
      saveRDS(data, checkpoint_file)
      if (verbose) message("💾 Checkpoint saved at row ", i)
    }
  }

  # 最后一次保存
  attr(data, "data_hash") <- data_hash
  saveRDS(data, checkpoint_file)

  data$CID <- as.integer(data$CID)

  found <- sum(!is.na(data$CID))
  not_found <- sum(is.na(data$CID))
  message("✅ CID extraction completed: ", found, " found, ", not_found, " not found.")

  # 提示音（需 beepr 包）
  if (requireNamespace("beepr", quietly = TRUE)) {
    beepr::beep(1)
  }

  return(data)
}


#-------------------------------------------------------------------------------
#' Extract metadata from PubChem for compounds with CIDs
#'
#' \code{extract_meta()} extracts chemical metadata from PubChem based on the
#' CID (Compound ID) column in the input data.frame or tibble. It retrieves
#' structural properties such as molecular formula, molecular weight, SMILES,
#' InChI, InChIKey, IUPAC name, and exact mass. Additionally, it can extract CAS
#' numbers, Flavornet sensory data, and synonyms from PubChem when requested.
#'
#' The function supports checkpointing to save intermediate results, allowing
#' long-running queries to resume from the last saved state without repeating
#' completed steps.
#'
#' @param data A data.frame or tibble containing at least a \code{CID} column
#'   (integer or numeric). Typically produced by \code{\link{extract_cid}}.
#' @param cas Logical, whether to extract CAS registry numbers. If \code{TRUE},
#'   CAS extraction is performed with fallback from PubChem sections to synonyms.
#' @param flavornet Logical, whether to extract Flavornet sensory descriptors.
#'   Requires \code{cas = TRUE} because it depends on CAS numbers.
#' @param synonyms Logical, whether to extract compound synonyms from PubChem.
#' @param verbose Logical, whether to print progress and status messages.
#' @param checkpoint_dir Character, directory path where checkpoint files are
#'   stored and loaded. Defaults to current directory \code{"."}.
#'
#' @return A data.frame or tibble identical to the input \code{data}, but with
#'   added columns for the extracted metadata:
#'   \itemize{
#'     \item Core PubChem properties: \code{Formula}, \code{MW}, \code{SMILES},
#'       \code{InChI}, \code{InChIKey}, \code{IUPAC_Name}, \code{Exact_Mass}.
#'     \item CAS number (in a new column \code{CAS_1} if \code{CAS} already exists).
#'     \item Flavornet sensory data in \code{Flavornet}.
#'     \item Compound synonyms in \code{Synonyms}.
#'   }
#'
#'   Intermediate results are saved as checkpoint files in \code{checkpoint_dir}.
#'   On re-running with the same input data, the function will resume from the
#'   checkpoint to avoid repeating completed queries.
#'
#' @details
#' This function requires internet access to query the PubChem and Flavornet
#' databases via the \pkg{webchem} package.
#'
#' The CAS number extraction tries first from PubChem's "section" data, and if
#' unsuccessful, falls back to searching the synonyms list for CAS-like strings.
#'
#' Checkpoints are matched via a hash of the input \code{data} to ensure
#' correspondence between saved intermediate results and current inputs.
#'
#' @seealso
#' \code{\link{extract_cid}}, \pkg{webchem}
#'
#' @examples
#' \dontrun{
#' library(webchem)
#' # Suppose df_cid is your data with a CID column already obtained by extract_cid()
#' df_meta <- extract_meta(df_cid, cas = TRUE, flavornet = TRUE, synonyms = TRUE)
#' }
#'
#' @importFrom digest digest
#' @importFrom webchem pc_prop pc_sect pc_synonyms fn_percept
#' @export
#'
extract_meta <- function(data,
                         cas = FALSE,
                         flavornet = FALSE,
                         synonyms = FALSE,
                         verbose = TRUE,
                         checkpoint_dir = ".") {

  if (!("CID" %in% names(data))) stop("Missing CID column. Please run extract_cid() first.")
  if (flavornet && !cas) stop("To extract Flavornet, please also set cas = TRUE.")

  message("🧪 Input data: ", nrow(data), " rows total.")
  message("📌 CID status: ", sum(!is.na(data$CID)), " non-empty, ", sum(is.na(data$CID)), " missing.")

  data$CID <- as.integer(data$CID)
  data_hash <- digest::digest(data)

  # --- Structural metadata ---
  if (verbose) message("\n🔍 Extracting core PubChem metadata...")
  props <- c("MolecularFormula", "MolecularWeight", "IsomericSMILES", "InChI", "InChIKey", "IUPACName", "ExactMass")
  cid_list <- data$CID
  cid_info <- tryCatch(
    webchem::pc_prop(na.omit(unique(cid_list)), properties = props),
    error = function(e) NULL
  )

  if (!is.null(cid_info) && nrow(cid_info) > 0) {
  col_map <- c(
    MolecularFormula = "Formula", MolecularWeight = "MW", IsomericSMILES = "SMILES",
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
        # 如果已经有这个列名，就加后缀 _1
        suffix_idx <- 1
        while (paste0(col, "_", suffix_idx) %in% names(data)) {
          suffix_idx <- suffix_idx + 1
        }
        new_col <- paste0(col, "_", suffix_idx)
        if (verbose) message(sprintf("⚠️ Column '%s' already exists, writing to '%s' instead.", col, new_col))
      }
      data[[new_col]] <- cid_info[[col]][idx]
    }
  }
}

  # --- CAS ---
  if (cas) {
    if (verbose) message("\n🔎 Extracting CAS (section + synonyms fallback).")
    cas_col <- if ("CAS" %in% names(data)) "CAS_1" else "CAS"
    if (!(cas_col %in% names(data))) data[[cas_col]] <- NA_character_
    cas_cp <- file.path(checkpoint_dir, "cas_checkpoint.rds")
    cas_source <- rep(NA_character_, nrow(data))

    checkpoint_matched <- FALSE
    if (file.exists(cas_cp)) {
      cp <- readRDS(cas_cp)
      if (!is.null(attr(cp, "data_hash")) && attr(cp, "data_hash") == data_hash) {
        if (verbose) message("✅ CAS checkpoint matched.")
        data[[cas_col]] <- cp$value
        cas_source <- cp$source
        checkpoint_matched <- TRUE
      } else if (verbose) message("⚠️ CAS checkpoint mismatch. Ignoring.")
    }

    found <- sum(!is.na(data[[cas_col]]))
    total <- nrow(data)
    for (i in seq_len(total)) {
      if (checkpoint_matched && !is.na(data[[cas_col]][i])) {
        # suppress "already present" messages when checkpoint matched
        next
      }
      if (checkpoint_matched && is.na(data$CID[i])) {
        # suppress messages when checkpoint matched and no CID
        next
      }
      if (!checkpoint_matched && verbose) message(sprintf("CAS: %d/%d", i, total))

      cid <- data$CID[i]
      compound_name <- if ("Name" %in% names(data)) data$Name[i] else paste0("row ", i)

      if (is.na(cid)) {
        if (!checkpoint_matched && verbose) message(sprintf("❌ No CID for row %d (%s), skipping CAS extraction.", i, compound_name))
        next
      }

      if (!is.na(data[[cas_col]][i])) {
        if (!checkpoint_matched && verbose) message(sprintf("⏭️ CAS already present for CID %s, skipping.", cid))
        next
      }

      val <- NA_character_
      source_type <- NA_character_

      sect_raw <- tryCatch(webchem::pc_sect(cid, "cas"), error = function(e) NULL)
      sect <- if (!is.null(sect_raw) && "Result" %in% names(sect_raw)) sect_raw$Result else NULL
      if (!is.null(sect) && length(sect) > 0 && nzchar(sect[[1]])) {
        val <- sect[[1]]
        source_type <- "section"
        found <- found + 1
        if (!checkpoint_matched && verbose) message(sprintf("✅ CAS %s found for CID %s (section).", val, cid))
      } else {
        syns <- tryCatch(webchem::pc_synonyms(cid, from = "cid"), error = function(e) NULL)
        if (is.list(syns)) {
          cas_like <- grep("^\\d{2,7}-\\d{2}-\\d$", unlist(syns), value = TRUE)
          if (length(cas_like) > 0) {
            val <- cas_like[1]
            source_type <- "synonym"
            found <- found + 1
            if (!checkpoint_matched && verbose) message(sprintf("✅ CAS %s found for CID %s (synonym).", val, cid))
          }
        }
      }

      data[[cas_col]][i] <- val
      cas_source[i] <- source_type
      if (is.na(val) && !checkpoint_matched && verbose) message(sprintf("❌ No CAS found for CID %s.", cid))

      saveRDS(structure(list(value = data[[cas_col]], source = cas_source), data_hash = data_hash), cas_cp)
    }

    message(sprintf("✅ CAS extraction done, Found: %d, Not found: %d", found, total - found))
  }

  # --- Flavornet ---
  if (flavornet) {
    if (verbose) message("\n🌸 Extracting Flavornet sensory data.")
    if (!"Flavornet" %in% names(data)) data$Flavornet <- NA_character_
    flav_cp <- file.path(checkpoint_dir, "flavornet_checkpoint.rds")

    checkpoint_matched <- FALSE
    if (file.exists(flav_cp)) {
      cp <- readRDS(flav_cp)
      if (!is.null(attr(cp, "data_hash")) && attr(cp, "data_hash") == data_hash) {
        if (verbose) message("✅ Flavornet checkpoint matched.")
        data$Flavornet <- cp$value
        checkpoint_matched <- TRUE
      } else if (verbose) message("⚠️ Flavornet checkpoint mismatch. Ignoring.")
    }

    cas_field <- if ("CAS_1" %in% names(data)) "CAS_1" else "CAS"
    found <- sum(!is.na(data$Flavornet))
    total <- nrow(data)
    for (i in seq_len(total)) {
      if (checkpoint_matched && !is.na(data$Flavornet[i])) {
        # suppress "already present" messages when checkpoint matched
        next
      }
      if (checkpoint_matched && is.na(data[[cas_field]][i])) {
        # suppress messages when checkpoint matched and no CAS
        next
      }
      if (!checkpoint_matched && verbose) message(sprintf("Flavornet: %d/%d", i, total))

      cas_val <- data[[cas_field]][i]
      if (is.na(cas_val)) {
        if (!checkpoint_matched && verbose) message(sprintf("❌ No CAS for row %d, skipping Flavornet extraction.", i))
        next
      }

      if (!is.na(data$Flavornet[i])) {
        if (!checkpoint_matched && verbose) message(sprintf("⏭️ Flavornet already present for CAS %s, skipping.", cas_val))
        next
      }

      val <- tryCatch(webchem::fn_percept(cas_val), error = function(e) NA_character_)
      if (!is.na(val)) {
        found <- found + 1
        if (!checkpoint_matched && verbose) message(sprintf("✅ Flavornet found for CAS %s: %s", cas_val, val))
      } else {
        if (!checkpoint_matched && verbose) message(sprintf("❌ No Flavornet found for CAS %s.", cas_val))
      }
      data$Flavornet[i] <- val

      saveRDS(structure(list(value = data$Flavornet), data_hash = data_hash), flav_cp)
    }
    message(sprintf("✅ Flavornet extraction done, Found: %d, Not found: %d", found, total - found))
  }

  # --- Synonyms ---
  if (synonyms) {
    if (verbose) message("\n📚 Extracting synonyms from PubChem.")
    if (!"Synonyms" %in% names(data)) data$Synonyms <- NA_character_
    syn_cp <- file.path(checkpoint_dir, "syn_checkpoint.rds")

    checkpoint_matched <- FALSE
    if (file.exists(syn_cp)) {
      cp <- readRDS(syn_cp)
      if (!is.null(attr(cp, "data_hash")) && attr(cp, "data_hash") == data_hash) {
        if (verbose) message("✅ Synonym checkpoint matched.")
        data$Synonyms <- cp$value
        checkpoint_matched <- TRUE
      } else if (verbose) message("⚠️ Synonym checkpoint mismatch. Ignoring.")
    }

    found <- sum(!is.na(data$Synonyms))
    total <- nrow(data)
    for (i in seq_len(total)) {
      if (checkpoint_matched && !is.na(data$Synonyms[i])) {
        # suppress "already present" messages when checkpoint matched
        next
      }
      if (checkpoint_matched && is.na(data$CID[i])) {
        # suppress messages when checkpoint matched and no CID
        next
      }
      if (!checkpoint_matched && verbose) message(sprintf("Synonyms: %d/%d", i, total))

      cid <- data$CID[i]
      if (is.na(cid)) {
        if (!checkpoint_matched && verbose) message(sprintf("❌ No CID for row %d, skipping synonyms extraction.", i))
        next
      }

      if (!is.na(data$Synonyms[i])) {
        if (!checkpoint_matched && verbose) message(sprintf("⏭️ Synonyms already present for CID %s, skipping.", cid))
        next
      }

      syns <- tryCatch(webchem::pc_synonyms(cid, from = "cid"), error = function(e) NULL)
      if (is.list(syns)) {
        syn_str <- paste(unique(unlist(syns)), collapse = "; ")
        if (nzchar(syn_str)) {
          data$Synonyms[i] <- syn_str
          found <- found + 1
          if (!checkpoint_matched && verbose) message(sprintf("✅ Synonyms found for CID %s.", cid))
          saveRDS(structure(list(value = data$Synonyms), data_hash = data_hash), syn_cp)
          next
        }
      }
      if (!checkpoint_matched && verbose) message(sprintf("❌ No synonyms found for CID %s.", cid))
      saveRDS(structure(list(value = data$Synonyms), data_hash = data_hash), syn_cp)
    }
    message(sprintf("✅ Synonyms extraction done, Found: %d, Not found: %d", found, total - found))
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
#' df <- data.frame(InChIKey = c("XMGQYMWWDOXHJM-UHFFFAOYSA-N"),
#'                  Name = c("Limonene"))
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
    cla <- data.frame(matrix(ncol = nrow(x), nrow = 0)) # tibble does not work
    cla <- rbind(cla, x$Classification)
    colnames(cla) <- x$Level
    return(cla)
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

    result <- suppressMessages(tryCatch({
      classyfireR::get_classification(inchikey)
    }, error = function(e) NULL))

    if (!is.null(result)) {
      class_info <- tryCatch({
        classyfireR::classification(result)
      }, error = function(e) NULL)

      class_df <- tryCatch({
        extract_cla(class_info)
      }, error = function(e) NULL)

      if (!is.null(class_df)) {
        message(sprintf("✅ Classification found via InChIKey for %s", compound_name))
        class_df[[inchikey_col]] <- inchikey
        classy_list[[inchikey]] <- class_df
      } else {
        message(sprintf("❌ Unclassified for %s", compound_name))
      }
    } else {
      message(sprintf("❌ Unclassified for %s", compound_name))
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
#' Note: compound names in your metadata file must exactly match those in your main dataset
#' (case-insensitive, whitespace-trimmed) for successful assignment.
#'
#' @param data A data.frame or tibble, typically the result after running \code{extract_meta()}.
#' @param meta_file A text-based file (e.g., *.csv or *.txt) containing at least "Name" and "SMILES",
#' and optionally "InChIKey" and "InChI". File can be created manually or generated from MOL files.
#'
#' @return A data.frame or tibble with missing SMILES, InChIKey, and InChI values filled in where available.
#'
#' @export
#'
assign_meta <- function(data, meta_file) {
  meta_data <- rio::import(meta_file)

  # 标准化列名
  names(data)[grepl("name", names(data), ignore.case = TRUE)] <- "Name"
  names(meta_data)[grepl("name", names(meta_data), ignore.case = TRUE)] <- "Name"
  names(meta_data)[grepl("smiles", names(meta_data), ignore.case = TRUE)] <- "SMILES"
  names(meta_data)[grepl("inchikey", names(meta_data), ignore.case = TRUE)] <- "InChIKey"
  names(meta_data)[grepl("^inchi$", names(meta_data), ignore.case = TRUE)] <- "InChI"

  # 标准化匹配键
  key_data <- trimws(tolower(data$Name))
  key_meta <- trimws(tolower(meta_data$Name))
  idx <- match(key_data, key_meta)

  # 如果没有 SMILES 列，创建一个
  if (!"SMILES" %in% names(data)) data$SMILES <- NA
  is_missing <- is.na(data$SMILES)
  data$SMILES[is_missing] <- meta_data$SMILES[idx[is_missing]]

  # 赋值 InChIKey
  if ("InChIKey" %in% names(meta_data)) {
    if (!"InChIKey" %in% names(data)) data$InChIKey <- NA
    is_missing <- is.na(data$InChIKey)
    data$InChIKey[is_missing] <- meta_data$InChIKey[idx[is_missing]]
  }

  # 赋值 InChI
  if ("InChI" %in% names(meta_data)) {
    if (!"InChI" %in% names(data)) data$InChI <- NA
    is_missing <- is.na(data$InChI)
    data$InChI[is_missing] <- meta_data$InChI[idx[is_missing]]
  }

  return(data)
}
