#' Read in MS-DIAL exported *.txt file and clean it
#'
#' @param file The *.txt file exported from MS-DIAL.
#' @param type Type of the data analyed by MS-DIAL, can be either gcms or lcms.
#' lcms is not currently implemented.
#' @param keep_unknown If TRUE, keep unknown compounds.
#' @param keep_spectrum If TRUE, keep EI_spectrum column.
#' @param keep_mean_sd If TRUE, keep mean and sd columns
#'
#' @return A cleaned table with most important information reserved.
#' @export
#'
#' @import dplyr
#' @importFrom janitor row_to_names clean_names
#' @importFrom tibble rowid_to_column
read_msdial <- function(file,
                        type = "gcms",
                        keep_unknown = FALSE,
                        keep_spectrum = FALSE,
                        keep_mean_sd = FALSE) {
  if (type == "gcms") {
    # read in the first 4 rows for compiling colnames
    tmp <- rio::import(file, nrows = 4, header = TRUE)
    col_names <- paste(tmp[4, ], tmp[3, ], sep = "_")
    col_names <- sub("_1$", "", col_names) # ensure remove only the last _1
    col_names <- sub("_Batch ID", "", col_names)
    col_names <- sub("_$", "", col_names)
    col_names <- sub("_Average", "_Mean", col_names)
    col_names <- sub("_Stdev", "_SD", col_names)
    col_names <- gsub(" |/|-|\\(|\\.", "_", col_names)
    col_names <- gsub("\\)|_%", "", col_names)
    col_names <- sub("__", "_", col_names)

    # read data
    data <- rio::import(file, skip = 4, header = TRUE)
    colnames(data) <- col_names # assign colnames
    data <- janitor::clean_names(data, case = "none")

    data <- data %>%
      select(2:5, 8:9, 11:12, 16, 19, 28:ncol(data)) %>%
      rename(
        RT = `Average_Rt_min`,
        RI = `Average_RI`,
        Reference_RI = `Reference_RI`,
        Quant_mass = `Quant_mass`,
        Name = `Metabolite_name`,
        Score = `Total_score`
      ) %>%
      relocate(Reference_RI, .after = RI) %>%
      relocate(Score, .after = Name) %>%
      relocate(Formula, .after = SMILES) %>%
      mutate(
        RT = as.numeric(RT) %>% round(digits = 2),
        RI = as.numeric(RI) %>% round(),
        Score = as.numeric(Score) %>% round(),
        Delta_RI = (RI - as.numeric(Reference_RI)) %>% round(),
        Quant_mass = as.numeric(Quant_mass) %>% round(),
        across(12:ncol(.), as.numeric),
        across(12:ncol(.), round)
      ) %>%
      suppressWarnings() %>%
      relocate(Delta_RI, .after = Reference_RI)
  }

  # remove unknown
  if (keep_unknown == FALSE) {
    data <- filter(data, Name != "Unknown")
  }
  # remove EI_spectrum
  if (keep_spectrum == FALSE) {
    data <- select(data, -EI_spectrum)
  }
  # remove mean and sd
  if (keep_mean_sd == FALSE) {
    data <- select(data, -contains("_Mean"), -contains("_SD"))
  }

  # convert rowid to column
  data <- data %>%
    tibble::rowid_to_column() %>%
    rename(ID = rowid)

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
calculate_freq <- function(
    data,
    num_sample,
    sep = ",") {
  # Define a function for calculate pct for a single string.
  get_freq <- function(x) {
    # str_split returns a list and has to be unlisted into a vector
    freq <- length(unlist(strsplit(x, sep)) %>% trimws())

    return(freq)
  }

  # Apply to the whole data set.
  data <- data %>%
    mutate(
      Detection_frequency = lapply(.$Comment, get_freq) %>% unlist(),
      Detection_rate = round(Detection_frequency / num_sample * 100, 1)
    ) %>%
    relocate(Detection_frequency:Detection_rate, .after = Comment)

  return(data)
}




# section -----------------------------------------------------------
merge_batches <- function(
    df1,
    df2,
    start_mz = 50,
    end_mz = 500,
    mz_step = 1,
    ri_tol_align = 5,
    ri_tol_iden = 30,
    sim_threshold = 0.8) {
  # check if df1 and df2 have duplicated names
  df1_duplicates <- df1$Name[duplicated(df1$Name)]
  if (length(df1_duplicates) != 0) {
    message("Duplicated Name found in the first dataset, please check:")
    print(df1_duplicates)
  } else {
    message("No duplicated Name found in the first dataset.")
  }

  df2_duplicates <- df2$Name[duplicated(df2$Name)]
  if (length(df2_duplicates) != 0) {
    message("Duplicated Name found in the second dataset, please check:")
    print(df2_duplicates)
  } else {
    message("No duplicated Name found in the second dataset.")
  }

  # check if delta RI is within the tolerance for both df1 and df2
  df1_delta_ri <- df1[abs(as.numeric(df1$Delta_RI)) > ri_tol_iden, c("RI", "Reference_RI", "Delta_RI", "Name", "Score")]
  if (nrow(df1_delta_ri) != 0) {
    message("Some compounds have delta RI greater than 30 for the first dataset, please check:")
    print(df1_delta_ri)
  } else {
    message("All compounds have delta RI within the tolerance for the first dataset.")
  }

  df2_delta_ri <- df2[abs(as.numeric(df2$Delta_RI)) > ri_tol_iden, c("RI", "Reference_RI", "Delta_RI", "Name", "Score")]
  if (nrow(df2_delta_ri) != 0) {
    message("Some compounds have delta RI greater than 30 for the second dataset, please check:")
    print(df2_delta_ri)
  } else {
    message("All compounds have delta RI within the tolerance for the second dataset.")
  }

  matched_df <- data.frame() # Initialize empty dataframe for results
  matched_id <- vector() # To keep track of matched IDs in df2

  # Loop through df1 rows
  for (i in seq_len(nrow(df1))) {
    row_df1 <- df1[i, ] # Extract current row from df1
    matched_rows <- data.frame() # Initialize empty dataframe for potential matches

    # Loop through df2 rows
    for (j in seq_len(nrow(df2))) {
      row_df2 <- df2[j, ]

      # Check RI tolerances
      if (abs(row_df1$RI - row_df2$RI) <= ri_tol_align) {
        # Calculate spectrum similarity
        spec1 <- update_spectrum(
          row_df1$EI_spectrum,
          start_mz = start_mz, end_mz = end_mz, mz_step = mz_step
        )
        spec2 <- update_spectrum(
          row_df2$EI_spectrum,
          start_mz = start_mz, end_mz = end_mz, mz_step = mz_step
        )
        spectrum_similarity <- round(cosine_similarity(spec1, spec2), 2)

        # If spectrum similarity is above the threshold, store match
        if (spectrum_similarity >= sim_threshold) {
          row_df2$Spectrum_Similarity <- spectrum_similarity
          matched_rows <- bind_rows(matched_rows, row_df2) # Add match to matched_rows
        }
      }
    }

    # If multiple matches, select the one with a same Name or the highest spectrum similarity
    if (nrow(matched_rows) > 1) {
      if (row_df1$Name %in% matched_rows$Name) {
        matched_rows <- matched_rows[matched_rows$Name == row_df1$Name, ]
      } else {
        matched_rows <- matched_rows[which.max(matched_rows$Spectrum_Similarity), ]
      }
    }

    # If a match is found, merge the row from df1 with the best match from df2
    if (nrow(matched_rows) != 0) {
      if (row_df1$Name != matched_rows$Name) {
        message(paste(
          paste0("Dataset 1, ID:", row_df1$ID, ","),
          paste0("RI:", row_df1$RI, ","),
          paste0("Name:", row_df1$Name, ","),
          paste0("Quant_mass:", row_df1$Quant_mass, ","),
          "------ Dataset 2 ID:", row_df2$ID,
          paste0("RI:", matched_rows$RI, ","),
          paste0("Name:", matched_rows$Name, ","),
          paste0("Quant_mass:", matched_rows$Quant_mass, ","),
          paste0("Spectrum Similarity:", matched_rows$Spectrum_Similarity, ","),
          paste0("Dealta RI:", row_df1$RI - matched_rows$RI, ","),
          "\nPlease check if they are the same compound. Type 'y' to merge, 'n' to separate, or 's' to stop."
        )) # nolint

        # Update spectrum of the matched rows
        spec2 <- update_spectrum(
          matched_rows$EI_spectrum,
          start_mz = start_mz, end_mz = end_mz, mz_step = mz_step
        )

        # diplaly the mirrored plot for easier comparison of spectra
        plot <- plot_mirrored_spectrum(spec1, spec2)
        print(plot)

        matched_rows$Spectrum_Similarity <- NULL
        # determine if the two entries are the same compounds
        user_input <- readline()
        if (user_input == "y") {
          # determine to choose which name should be kept, 1st or 2nd
          message("Choose if the first name or second to be used after compilation. Type '1' to choose the first name, or '2' to choose the second name.") # nolint
          user_input2 <- readline()

          if (user_input2 == "1") {
            message("The first Name will be used.\n")
            # Create combined row with information from both df1 and df2
            combined_row <- data.frame(
              ID = paste(row_df1$ID, sub(".*_", "", deparse(substitute(df2))), sep = "_"),
              RT = row_df1$RT,
              RI = row_df1$RI,
              Delta_RI = row_df1$Delta_RI,
              Reference_RI = row_df1$Reference_RI,
              Quant_mass = row_df1$Quant_mass,
              Name = row_df1$Name,
              Score = row_df1$Score,
              INCHIKEY = row_df1$INCHIKEY,
              SMILES = row_df1$SMILES,
              Formula = row_df1$Formula,
              Comment = paste(row_df1$Comment, matched_rows$Comment, sep = ","),
              EI_spectrum = row_df1$EI_spectrum
            )
          } else {
            message("The second Name will be used.\n")
            # Create combined row with information from both df1 and df2
            combined_row <- data.frame(
              ID = paste(row_df1$ID, sub(".*_", "", deparse(substitute(df2))), sep = "_"),
              RT = matched_rows$RT,
              RI = matched_rows$RI,
              Delta_RI = matched_rows$Delta_RI,
              Reference_RI = matched_rows$Reference_RI,
              Quant_mass = matched_rows$Quant_mass,
              Name = matched_rows$Name,
              Score = matched_rows$Score,
              INCHIKEY = matched_rows$INCHIKEY,
              SMILES = matched_rows$SMILES,
              Formula = matched_rows$Formula,
              Comment = paste(row_df1$Comment, matched_rows$Comment, sep = ","),
              EI_spectrum = matched_rows$EI_spectrum
            )
          }

          # Append additional columns from both rows
          combined_row <- cbind(combined_row, row_df1[1, 14:ncol(row_df1)], matched_rows[1, 14:ncol(matched_rows)]) # nolint

          matched_df <- bind_rows(matched_df, combined_row) # Add the combined row to the result
          matched_id <- c(matched_id, matched_rows$ID) # Track matched df2 ID
        } else if (user_input == "n") {
          # No match found, combine df1 row with empty df2 columns (NA)
          empty_df2 <- row_df2
          empty_df2[] <- NA # Fill all columns with NA
          row_df1$ID <- paste(row_df1$ID, sub(".*_", "", deparse(substitute(df2))), sep = "_")
          # Add empty df2 columns to df1 row
          row_df1 <- cbind(row_df1, empty_df2[1, 14:ncol(empty_df2)])
          matched_df <- bind_rows(matched_df, row_df1) # Add the row to the result
        } else if (user_input == "s") {
          stop("Exit!")
        }
      } else {
        matched_rows$Spectrum_Similarity <- NULL
        combined_row <- data.frame(
          ID = paste(row_df1$ID, sub(".*_", "", deparse(substitute(df2))), sep = "_"),
          RT = row_df1$RT,
          RI = row_df1$RI,
          Delta_RI = row_df1$Delta_RI,
          Reference_RI = row_df1$Reference_RI,
          Quant_mass = row_df1$Quant_mass,
          Name = row_df1$Name,
          Score = max(c(row_df1$Score, matched_rows$Score)),
          INCHIKEY = row_df1$INCHIKEY,
          SMILES = row_df1$SMILES,
          Formula = row_df1$Formula,
          Comment = paste(row_df1$Comment, matched_rows$Comment, sep = ","),
          EI_spectrum = row_df1$EI_spectrum
        )
        # Append additional columns from both rows
        combined_row <- cbind(combined_row, row_df1[1, 14:ncol(row_df1)], matched_rows[1, 14:ncol(matched_rows)]) # nolint

        matched_df <- bind_rows(matched_df, combined_row) # Add the combined row to the result
        matched_id <- c(matched_id, matched_rows$ID) # Track matched df2 ID
      }
    } else {
      # No match found, combine df1 row with empty df2 columns (NA)
      empty_df2 <- row_df2
      empty_df2[] <- NA # Fill all columns with NA
      row_df1$ID <- as.character(row_df1$ID)
      row_df1 <- cbind(row_df1, empty_df2[1, 14:ncol(empty_df2)]) # Add empty df2 columns to df1 row
      matched_df <- bind_rows(matched_df, row_df1) # Add the row to the result
    }
  }

  # Handle unmatched rows in df2 (i.e., those that weren't paired with any df1 row)
  unmatched_df2 <- df2[!(df2$ID %in% matched_id), ]
  if (nrow(unmatched_df2) > 0) {
    unmatched_df2$ID <- sub(".*_", "", deparse(substitute(df2)))
    matched_df <- bind_rows(matched_df, unmatched_df2) # Add unmatched df2 rows to the result
  }

  matched_df <- matched_df[order(matched_df$RT), ]
  # renames rownames with numbers
  rownames(matched_df) <- seq_len(nrow(matched_df))

  duplicates <- matched_df$Name[duplicated(matched_df$Name)]
  if (length(duplicates) != 0) {
    message("Duplicated Name found in the merged dataset, please check:")
    print(duplicates)
  } else {
    message("No duplicated Name found in the merged dataset.")
  }

  message(
    "Merged dataset created!!\n",
    "There are ", nrow(matched_df), " entries in the merged dataset.\n"
  )

  return(matched_df)
}


###############################################################################
# define a function to merge duplicate rows

merge_duplicates <- function(df) {
  # Calculate duplicates
  duplicate_names <- unique(df$Name[duplicated(df$Name)])
  n_duplicates <- length(duplicate_names)

  if (n_duplicates > 0) {
    message(
      "Found ", n_duplicates, " duplicate(s):\n",
      paste(duplicate_names, collapse = "\n")
    )

    merged_data <- df %>%
      group_by(Name) %>%
      summarise(
        ID = paste(unique(ID), collapse = "_"), # Combine unique IDs
        Comment = paste(unique(Comment), collapse = ","), # Combine comments
        across(1:13, first), # Always select the first row for columns 1 to 12
        across(14:(ncol(.) - 1), \(x) sum(x, na.rm = TRUE))
      ) %>%
      ungroup() %>%
      arrange(RT) %>%
      relocate(Comment, .after = Formula) %>%
      relocate(Name, .before = Score)

    message("Duplicates merged. The merged dataset contains ", nrow(merged_data), " records.")
    return(merged_data)
  } else {
    message("No duplicates found.")
    return(df)
  }
}