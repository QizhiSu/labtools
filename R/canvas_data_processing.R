#' Manually align Canvas data
#'
#' \code{read_canvas} Facilitate the alignment of canvas data by directly
#' reading .*txt files export from Canvas.
#'
#' Canvas is currently only capable of processing data one by one, which is error-prone
#' and tedious to align and combine data from different samples. This function
#' aims to facilitate this process and help us find conflicting identification
#' among different samples. The first step is to manually identify and mark peaks
#' of interest in Canvas. Then, exported the marked peaks as .txt from all samples
#' into a folder. The function will then read all .txt files and combine them into
#' a single data by matching the chemical name. Based on the defined tolerance of
#' retention index and 2nd dimension retention time, the function will evaluate
#' every chemical compound and then prompt up message about in which samples the
#' compounds are detected and which samples have the maximum and minimum RI, with
#' the aim to help identify which samples could have mis-identification results.
#' In order to be reproducible and tracible, mis-identification results should be
#' modified in Canvas and re-exported into .txt replacing the previous ones.
#' Manual modification of the .txt or the combined data is not recommended.
#'
#'
#'
#' @param path Path of the folder contains .*txt files originally exported from
#' Canvas
#' @param ri_align_tolerance RI tolerance acceptable for alignment
#' @param rt_2d_tolerance 2nd dimension RT acceptable for alignment
#' @param keep whether 'area', 'height', or 'both', telling which information to
#' be kept. If keep = 'area', then only peak areas of the compounds in detected
#' samples will be kept; if keep = 'height', then only peak height of the compounds
#' in detected samples will be kept; if keep = 'both', then both peak areas and
#' peak height of the compounds in detected samples will be kept.
#'
#' @return A \code{data.frame} with compound name, average RT, RI, 2d RT, CAS, as
#' well as peak are/height in each sample.
#' @export
#'
#' @import dplyr
#' @import purrr
#' @importFrom glue glue
#' @importFrom utils capture.output read.table
read_canvas <- function(path,
                        ri_align_tolerance = 50,
                        rt_2d_tolerance = 0.05,
                        keep = "area") { # can be 'area', 'height', or 'both'
  ## A function to read and organize single canvas data
  read_single <- function(file) {
    data <- read.table(
      file,
      header = TRUE,
      sep = '\t',
      fileEncoding = 'gbk',
      col.names = c(
        "Number",
        "Name",
        "RT",
        "RT_2D",
        "Area",
        "Height",
        "Percentage",
        "Peak_width",
        "CAS",
        "Formula",
        "Match",
        "R_match",
        "Possibility",
        "RI",
        "NIST_RI",
        "MW",
        "Note"
      )
    ) %>%
      dplyr::select(2:6, 9:16) %>%
      filter(Name != "") %>%
      mutate(Score = round(rowMeans(select(., Match, R_match))),
             RT = round(as.numeric(RT), 3),
             RT_2D = round(as.numeric(RT_2D), 3),
             Area = round(as.numeric(Area)),
             Height = round(as.numeric(Height)),
             CAS = ifelse(CAS == '0-0-0', NA, CAS),
             Match = as.numeric(Match),
             R_match = as.numeric(R_match),
             Possibility = round(as.numeric(Possibility)),
             RI = as.numeric(RI),
             MW = as.numeric(MW)) %>%
      relocate(Score, .after = R_match) %>%
      arrange(RT)

    return(data)
  }


#-------------------------------------------------------------------------------

  ## read all files into a list
  # get file paths of all files in the path
  files_path <- list.files(path, pattern = '.*txt', full.names = TRUE)
  file_name <- sub("\\.txt", "", basename(files_path)) # get the file name
  file_name <- gsub("-", '_', file_name) # replace "-" with '_'
  files <- lapply(files_path, read_single) # read all files into a list
  names(files) <- file_name # rename list elements with their correspondent file name
  # add file name to column names of each list element
  files <-
    purrr::imap(files, ~ rename_with(.x, function(z)
      paste(z, .y, sep = "_#"),-c("Name")))

#-------------------------------------------------------------------------------

  ## Combine data from all samples.
  # At present by 'Name', can be further improved in the future
  cb_data <-
    files %>%
    purrr::reduce(dplyr::full_join, by = "Name") %>% # combine into a dataframe
    mutate(RI_SD = select(., matches("^RI")) %>% # sample name should not have "RI" at the beginning
             purrr::pmap(~ round(sd(c(...), na.rm = TRUE))) %>%
             unlist(),
           RT_2D_SD = select(., matches("^RT_2D")) %>%
             purrr::pmap(~ round(sd(c(...), na.rm = TRUE), 3)) %>%
             unlist())

#-------------------------------------------------------------------------------

   ## A function to message or warnings a dataframe
  print_and_capture <- function(x)
  {
    paste(capture.output(print(x)), collapse = "\n")
  }

#-------------------------------------------------------------------------------

  ## check alignment and warnings probable misalignment
  # evaluate compound by compound to see if different sample could have significantly
  # different RI or RT_2D
  for (i in 1:nrow(cb_data)) {
    if(!is.na(cb_data$RI_SD[i])) { # if RI_SD is.na, it suggests that only detected in one sample
      if(cb_data$RI_SD[i] > ri_align_tolerance) {
        # if RI_SD is smaller than the defined value, message which samples could have problem
        samples <-
          cb_data %>%
          select(matches('^RI') & -contains('SD')) %>% # subset RI columns
          slice(i) %>% # keep only the i row
          select(which(!is.na(.))) %>% # keep non na columns
          rename_with(~gsub('RI_#?', '', .))

          sample_name <- # sample names into a vector
            samples %>%
            colnames() %>%
            paste(collapse = ', ')
          sample_max <- colnames(samples)[which.max(samples)]
          sample_min <- colnames(samples)[which.min(samples)]

          message( # message all samples that have detected the compound
            # and which sample have maximum and minimum value respectively
            glue::glue(
              'Compound {cb_data$Name[i]} was detected in {ncol(samples)} ',
              'samples, but the retention index varied significantly across samples ',
              'as detailed below:\n',
              '{print_and_capture(samples)}\n',
              'The average RI was {rowMeans(samples)}. Sample {sample_max} had ',
              'maximum RI, while {sample_min} had minimum.\n',
              '**************************************************************\n\n'
            )
          )
      }else{
          if(cb_data$RT_2D_SD[i] > rt_2d_tolerance) {
            samples <-
              cb_data %>%
              select(matches('^RT_2D') & -contains('SD')) %>% # subset RI columns
              slice(i) %>% # keep only the i row
              select(which(!is.na(.))) %>% # keep non na columns
              rename_with(~gsub('RT_2D_#?', '', .))

            sample_name <- # sample names into a vector
              samples %>%
              colnames() %>%
              paste(collapse = ', ')
            sample_max <- colnames(samples)[which.max(samples)]
            sample_min <- colnames(samples)[which.min(samples)]

            message(
              glue::glue(
                'Compound {cb_data$Name[i]} was detected in {ncol(samples)} ',
                'samples, but the 2nd dimesional retention time varied significantly ',
                'across samples as detailed below:\n',
                '{print_and_capture(samples)}\n',
                'The average RI was {rowMeans(samples)}. Sample {sample_max} had ',
                'maximum RI, while {sample_min} had minimum.\n',
                '**************************************************************\n\n'
              )
            )
          }
        }
      }
  }

#-------------------------------------------------------------------------------

  ## simplify the combined data
  sim_data <-
    cb_data %>%
    mutate(RT = select(., contains("RT_") & -contains("RT_2D")) %>%
             rowMeans(na.rm = TRUE) %>%
             round(3),
           RT_2D = select(., contains("RT_2D")) %>%
             rowMeans(na.rm = TRUE) %>%
             round(3),
           RI = select(., matches("^RI") & -contains("SD")) %>%
             rowMeans(., na.rm = TRUE) %>%
             round(0),
           NIST_RI = coalesce(!!!select(., contains("NIST_RI"))),
           Match = coalesce(!!!select(., contains("Match"))),
           R_match = coalesce(!!!select(., contains("R_match"))),
           CAS = coalesce(!!!select(., contains("CAS"))),
           Score = pmax(!!!select(., contains("Score")), na.rm = TRUE),
           Count = paste0(rowSums(!is.na(select(., contains("Area")))),
                          "/", length(samples)),
           Ratio = round(rowSums(!is.na(select(., contains("Area"))))/
                           length(samples) * 100, digits = 1)
    ) %>%
    select(-contains("#"),
           contains(c("Area", "Height")),
           -c('RI_SD', 'RT_2D_SD')) %>%
    relocate(CAS, .before = Match)


  ## which information to be kept
  if(keep == 'area') {
    sim_data <-
      sim_data %>%
      select(-contains('#'), contains('Area')) %>%
      rename_with(~sub('Area_#', '', .))
  }
  if(keep == 'height') {
    sim_data <-
      sim_data %>%
      select(-contains('#'), contains('Height')) %>%
      rename_with(~sub('Height_#', '', .))
  }

  return(sim_data)
  }


