#' Export structure database for MS-FINDER
#'
#' \code{export4msfinder} provides a way to convert any list of chemicals into a
#' structure database that can be used by MS-FINDER.
#'
#' @param input The path of the list of chemicals. The file can be either in txt,
#' csv, xlsl, or xls form, but the chemical name should be in the first column. #'
#' for example, "c:/data/list_of_chemicals.txt"
#' @param cas_col If you have CAS number for each chemical, then specify in which
#' column the CAS number locates. It will be used for extracting chemical information
#' from Pubchem. This is recommended. If no CAS number is available, just levea it,
#' and it will use the chemical name for extraction.
#' @param output The path where the structure database will be stored. It must be
#' in txt format. For example, "c:/data/structure_database_for_msfinder.txt"
#'
#' @export
#'
#' @importFrom rio import export
#' @importFrom webchem get_cid pc_prop
#' @importFrom stats setNames
export4msfinder <- function(input,
                            cas_col = FALSE,
                            output = "structure data for msfinder.txt") {
  data <- rio::import(input)
  if (!"CID" %in% colnames(data)) data$CID <- NA_character_

  message("Extracting chemical information from Pubchem.\n",
          "Please waite...\n ",
          "...")
  # cas is optional, but name is necessary
  if (cas_col != FALSE) {
    data <-
      data %>%
      mutate(CID = case_when(
        is.na(CID) ~
          webchem::get_cid(.[, cas_col], match = "first")$cid %>% as.character(),
        TRUE ~ CID
      ))
  }

  data <-
    data %>%
    mutate(CID = case_when(
      is.na(CID) ~
        webchem::get_cid(.[, 1], from = "name", match = "first")$cid %>%
        as.character(),
      TRUE ~ CID
    ))

  # extract meta data
  data <-
    data %>%
    mutate(CID = as.integer(CID)) %>%
    select(1, CID)

  data <-
    data %>%
    filter(!is.na(CID)) %>%
    distinct(CID) %>%
    webchem::pc_prop(properties = c("IsomericSMILES",
                           "InChIKey",
                           "ExactMass",
                           "MolecularFormula")) %>%
    left_join(select(data, 1, CID), ., by = "CID") %>%
    mutate(Short_InChIKey = sub(".{13}$", "", InChIKey),
           Database_ID = dplyr::row_number()) %>%
    relocate(InChIKey,
             Short_InChIKey,
             CID,
             ExactMass,
             MolecularFormula,
             IsomericSMILES,
             .before = Database_ID) %>%
    setNames(c("Title",
             "InChIKey",
             "Short InChIKey",
             "Pubchem CID",
             "Exact mass",
             "Formula",
             "SMILES",
             "Database ID"))
  message("Extraction finished")

  rio::export(data, output)
}
