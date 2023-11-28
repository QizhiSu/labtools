#' Read in MS-DIAL exported *.txt file and clean it
#'
#' @param file The *.txt file exported from MS-DIAL.
#' @param type Type of the data analyed by MS-DIAL, can be either gcms or lcms.
#' lcms is not currently implemented.
#'
#' @return A cleaned table with most important information reserved.
#' @export
#'
#' @import dplyr
#' @importFrom janitor row_to_names clean_names
#' @importFrom tibble rowid_to_column
read_msdial <- function(file,
                        type = "gcms") {
  if(type == "gcms") {
    data <- rio::import(file, skip = 3, header = TRUE) %>%
      janitor::row_to_names(1) %>%
      suppressWarnings() %>%
      janitor::clean_names(case = "none")

    data <- data %>%
      select(2:5, 8:9, 11:12, 16, 19, 29:ncol(data)) %>%
      select(-contains(c("BK", "bk", "blank", "Blank", "BLANK"))) %>%
      select(!matches("^X\\d")) %>%
      rename(RT = `Average_Rt_min`,
             RI = `Average_RI`,
             Reference_RI = `Reference_RI`,
             Quant_mass = `Quant_mass`,
             Name = `Metabolite_name`,
             Score = `Total_score`) %>%
      relocate(Reference_RI, .after = RI) %>%
      relocate(Score, .after = Name) %>%
      relocate(Formula, .after = SMILES) %>%
      filter(Name != "Unknown") %>%
      mutate(RT = as.numeric(RT) %>% round(digits = 2),
             RI = as.numeric(RI) %>% round(),
             Score = as.numeric(Score) %>% round(),
             Delta_RI = (RI - as.numeric(Reference_RI)) %>% round(),
             Quant_mass = as.numeric(Quant_mass) %>% round(),
             across(11:ncol(.), as.numeric),
             across(11:ncol(.), round)) %>%
      relocate(Delta_RI, .after = Reference_RI) %>%
      tibble::rowid_to_column() %>%
      rename(ID = rowid)
  }

  return(data)
}



#' Calculate the detection frequency and detection rate
#'
#' The frequency and rate of samples that have a compound detected in total number of
#' samples. It is calculated in terms of samples excluding replicates.
#'
#' @param data The data returned by the \code{read_msdial()} function.
#' @param num_sample Total number of samples.
#' @param sep The separator between samples in the Comment field in MS-DIAL.
#' For example, if you put "1, 2, 3, 4, 5, 6" in the Comment field in MS-DIAL,
#' where each number denoted as sample code. Then the separator is ",".
#'
#' @return A table with detection frequency and rate added after Comment column.
#' @export
#'
#' @import dplyr
#'
calculate_freq <- function(data,
                          num_sample,
                          sep = ","){
  # Define a function for calculate pct for a single string.
  get_freq <- function(x) {
    # str_split returns a list and has to be unlisted into a vector
    freq <- length(unlist(strsplit(x, sep)) %>% trimws())

    return(freq)
  }

  # Apply to the whole data set.
  data <- data %>%
    mutate(Detection_frequency = lapply(.$Comment, get_freq) %>% unlist(),
           Detection_rate = round(Detection_frequency/num_sample*100, 1)) %>%
    relocate(Detection_frequency:Detection_rate, .after = Comment)

  return(data)
}

