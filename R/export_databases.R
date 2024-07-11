#' Export Database for MS-FINDER
#'
#' \code{export4msfinder} provides a way to convert any list of chemicals into a
#' structure database that can be used by MS-FINDER. please refer to
#'  \link{extract_cid} and \link{extract_meta} regarding meta data extraction.
#'
#' @param data The list of chemicals should contain at least Name, SMILES, CID,
#' InChIKey, Formula.
#' @param output The path where the structure database will be stored. It must be
#' in txt format. For example, "c:/data/database for msfinder.txt"
#'
#' @export
#'
#' @import dplyr
#' @importFrom stats setNames
#' @importFrom utils write.table
#'
export4msfinder <- function(data,
                            output = "database for msfinder.txt") {
  column_names <- c("Name", "InChIKey", "CID", "ExactMass", "Formula", "SMILES")
  if (all(column_names %in% colnames(data)) == FALSE) {
    stop(
      "Your data does not have Name, InChIKey, CID, ExactMass, Formula, and/or SMILES ",
      "columns, please double-check your data. You can change them accordingly ",
      "or refer to the extract_cid and extract_meta functions."
    )
  }

  data <-
    data %>%
    dplyr::mutate(
      Short_InChIKey = sub(".{13}$", "", InChIKey),
      Database_ID = dplyr::row_number()
    ) %>%
    dplyr::select(
      Name,
      InChIKey,
      Short_InChIKey,
      CID,
      ExactMass,
      Formula,
      SMILES,
      Database_ID
    ) %>%
    setNames(c(
      "Title",
      "InChIKey",
      "Short InChIKey",
      "Pubchem CID",
      "Exact mass",
      "Formula",
      "SMILES",
      "Database ID"
    ))

  # write txt
  write.table(data, output, sep = "\t", row.names = FALSE)
}

#-------------------------------------------------------------------------------
#' Export database for MS-DIAL
#'
#' \code{export4msdial} provides a way to convert any list of chemicals into a
#' structure database that can be used by MS-DIAL for post-identification.
#' please refer to \link{extract_cid} and \link{extract_meta} regarding meta
#' data extraction.

#' @param data The list of chemicals should contain at least Name and ExactMass
#' @param polarity The ESI polarity applied in your dat. Either "pos" or "neg".
#' @param output The path where the structure database will be stored. It must be
#' in txt format. For example, "c:/data/database for msdial.txt"
#'
#' @export
#'
#' @import dplyr
#' @importFrom utils write.table
#'
export4msdial <- function(data,
                          polarity = "pos",
                          output = "database for msdial.txt") {
  column_names <- c("Name", "ExactMass", "SMILES", "InChIKey")
  if (all(column_names %in% colnames(data)) == FALSE) {
    stop(
      "Your data does not have Name, ExactMass, InChIKey, Formula, SMILES, and/or SMILES columns",
      "please double-check your data. You can change them accordingly ",
      "or refer to the extract_cid and extract_meta functions."
    )
  }

  # subset data
  data <- data[, c("Name", "ExactMass", "InChIKey", "Formula", "SMILES")]
  colnames(data) <- c("Metabolite", "MZ", "InChIKey", "Formula", "SMILES")
  data$RT <- NA

  # keep number of rows for later use
  n_row <- nrow(data)

  # assign adduct
  if (polarity == "pos") {
    # repeat rows 8 times
    data <- data[rep(seq_len(nrow(data)), times = 8), ]
    rownames(data) <- NULL

    # define adducts
    adducts <- c(
      "[M+H]+", "[M+Na]+", "[M+NH4]+", "[M+K]+",
      "[2M+H]+", "[2M+Na]+", "[2M+NH4]+", "[2M+K]+"
    )
    data$Adduct <- rep(adducts, each = n_row)

    data <-
      data %>%
      mutate(MZ = as.numeric(MZ)) %>%
      mutate(
        MZ = ifelse(Adduct == "[M+H]+", MZ + 1.007276, MZ),
        MZ = ifelse(Adduct == "[M+Na]+", MZ + 22.989218, MZ),
        MZ = ifelse(Adduct == "[M+NH4]+", MZ + 18.033823, MZ),
        MZ = ifelse(Adduct == "[M+K]+", MZ + 38.963158, MZ),
        MZ = ifelse(Adduct == "[2M+H]+", 2 * MZ + 1.007276, MZ),
        MZ = ifelse(Adduct == "[2M+Na]+", 2 * MZ + 22.989218, MZ),
        MZ = ifelse(Adduct == "[2M+NH4]+", 2 * MZ + 18.033823, MZ),
        MZ = ifelse(Adduct == "[2M+K]+", 2 * MZ + 38.963158, MZ)
      )
  } else {
    # for negative mode
    # repeat rows 8 times
    data <- data[rep(seq_len(nrow(data)), times = 5), ]
    rownames(data) <- NULL

    # define adducts
    adducts <- c("[M-H]-", "[M-H2O-H]-", "[M+FA-H]-", "[2M-H]-", "[2M+FA-H]-")
    data$Adduct <- rep(adducts, each = n_row)

    data <-
      data %>%
      mutate(MZ = as.numeric(MZ)) %>%
      mutate(
        MZ = ifelse(Adduct == "[M-H]-", MZ - 1.007276, MZ),
        MZ = ifelse(Adduct == "[M-H2O-H]-", MZ - 19.01839, MZ),
        MZ = ifelse(Adduct == "[M+FA-H]-", MZ + 44.998201, MZ),
        MZ = ifelse(Adduct == "[2M-H]-", 2 * MZ - 1.007276, MZ),
        MZ = ifelse(Adduct == "[2M+FA-H]-", 2 * MZ + 44.998201, MZ)
      )
  }

  data <- relocate(data, RT:Adduct, .after = MZ)
  # write txt
  write.table(data, output, sep = "\t", row.names = FALSE)
}
