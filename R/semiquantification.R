#' Select Standard
#'
#' This function selects the standard used for semi-quantification of a given
#' set of substances based on chemical structure similarity calculations between
#' the standards avaiable and the target substances
#'
#' @param std_md The data for the standards including at least the Name,
#' response variable for each standard (can be either slope of the calibration
#' curve, peak area, or any other response variable), and molecular descriptor
#' for each standard. The molecular descriptor data should be at the end of the
#' dataframe
#' @param std_res_col The column index of the response variable in std_md.
#' If you want to use all molecular descriptors, set std_res_col to FALSE
#' @param std_md1_col The column index where the molecular descriptor starts
#' in std_md
#' @param data_md The data to be determined including at least the Name and
#' molecular descriptor for each standard. The molecular descriptor data should
#' be at the end of the dataframe
#' @param data_md1_col The column index where the molecular descriptor data
#' starts in data_md
#' @param top_npct_md The percentage of top molecular descriptors to select
#' (default is 20)
#'
#' @return A metadata dataframe with the selected standard assigned to each
#' substance and their similarity value
#'
#' @importFrom randomForest randomForest
#' @importFrom randomForest importance
#' @importFrom proxy simil
#' @importFrom stats na.omit
#'
#' @export
select_std <- function(
    std_md,
    std_res_col,
    std_md1_col,
    data_md,
    data_md1_col,
    top_npct_md = 20) {
  # rename the columns to avoid errors in randomForest
  colnames(std_md) <- gsub("-", "_", colnames(std_md))
  colnames(data_md) <- gsub("-", "_", colnames(data_md))

  # construct the data for feature selection
  std4rf <- std_md[, c(std_res_col, std_md1_col:ncol(std_md))]
  colnames(std4rf)[1] <- "response" # rename the response column for later use
  # remove empty molecular descriptor
  std4rf <- std4rf[!apply(std4rf[, -1], 1, function(x) all(is.na(x))), ]
  std4rf <- std4rf[, !apply(std4rf, 2, function(x) any(is.na(x)))]

  # if you want to use all molecular descriptors, set std_res_col to FALSE
  if (std_res_col == FALSE) {
    vip_df <- data.frame(MD = colnames(std4rf))
  } else {
    # build a random forest model for calculation of importance metric
    rf_model <- randomForest::randomForest(response ~ .,
      data = std4rf,
      importance = TRUE,
      na.action = na.omit
    )

    # calculate the importance metric for each molecular descriptor
    vip_df <- randomForest::importance(rf_model) |> as.data.frame()
    vip_df <- vip_df[order(-vip_df$`%IncMSE`), ] # sort by decreasing importance
    # select top n per molecular descriptors
    vip_df <- vip_df[1:ceiling(top_npct_md / 100 * nrow(vip_df)), ]
    vip_df$MD <- rownames(vip_df)
  }

  # keep only top k molecular descriptors in the data and standard
  data_md <- data_md[, c(
    1:data_md1_col - 1,
    which(names(data_md) %in% vip_df$MD)
  )]
  std_md <- std_md[, c(1:std_md1_col - 1, which(names(std_md) %in% vip_df$MD))]


  # create a metadata dataframe for the selected standard
  data_meta <- data_md[, c(1:data_md1_col - 1)]
  data_meta$Standard <- NA
  data_meta$Similarity <- NA

  cat("Assigning standard for target compounds...\n\n")
  # calculate the similarity between the selected standard and each substance
  for (i in seq_len(nrow(data_md))) {
    cat(i, "in", nrow(data_md), "\n")
    tmp <- vector() # create a temporary vector to store tempory similarity
    data_meta$Standard[i] <- NA
    data_meta$Similarity[i] <- NA

    # take the molecular descriptor of the i substance
    a <- data_md[i, data_md1_col:ncol(data_md)]
    a_col <- names(a)[!is.na(a)] # remove empty molecular descriptor
    a <- a[, a_col]

    # calculate the similarity one by one
    for (j in seq_len(nrow(std_md))) {
      # take the molecular descriptor of the j standard
      b <- std_md[j, std_md1_col:ncol(std_md)]
      b <- b[, a_col] # keep only the common molecular descriptor
      # Calculate similarity
      tmp[j] <- proxy::simil(a, b)
    }

    # assign the standard with the highest similarity
    data_meta$Standard[i] <- std_md$Name[which.max(tmp)]
    data_meta$Similarity[i] <- max(tmp)
  }

  # remove rownames
  rownames(data_meta) <- NULL
  data_meta <- left_join(
    data_meta,
    std_md[, c("Name", "Slope", "Intercept", "R2", "Low_con", "Low_area")],
    by = c("Standard" = "Name")
  )

  cat("\nStandard assignment completed.\n")
  return(data_meta)
}

# section -----------------------------------------------------------
#' Keep specified peak areas in the dataframe
#'
#' This function processes a dataframe to keep only the specified peak areas based on sample codes
#' in the Comment column and optional flags for 'BK' and 'D8'.
#'
#' @param df A dataframe containing the data to be processed.
#' @param sam_code A dataframe containing sample codes and their corresponding names.
#' @param start_col The starting column index of peak area columns (default is 12).
#' @param keep_bk A logical indicating whether to keep 'BK' in the list of samples (default is
#' TRUE).
#' @param keep_d8 A logical indicating whether to keep 'D8' in the list of samples (default is
#' TRUE).
#'
#' @return A dataframe with the specified areas kept and others set to NA.
#'
#' @export
keep_area <- function(
    df,
    sam_code,
    start_col = 12,
    keep_bk = TRUE,
    keep_d8 = TRUE) {
  # Loop through each row of the dataframe to process comments and update areas
  for (i in seq_len(nrow(df))) {
    tryCatch(
      {
        tmp <- unlist(strsplit(df$Comment[i], ","))
      },
      warning = function(w) {
        cat("Warning at i = ", i, "\n")
        cat("Warning message:", conditionMessage(w), "\n")
        return(NULL) # Return NULL to handle potential errors gracefully
      }
    )

    # Convert sample codes to sample names using a lookup table
    tmp1 <- sam_code$Name[match(tmp, sam_code$Code)]
    # Optionally add 'BK' and 'D8' to the list of samples to keep
    if (keep_bk) {
      tmp1 <- c(tmp1, "BK")
    }
    if (keep_d8) {
      tmp1 <- c(tmp1, "D8")
    }

    # Identify columns corresponding to the compounds of interest
    comp_index <- which(grepl(
      paste(tmp1, collapse = "|"),
      colnames(df[i, start_col:ncol(df)])
    ))
    # Determine the total number of entries for the current row
    num_entries <- ncol(df[i, start_col:ncol(df)])
    tmp3 <- rep(0, num_entries)
    tmp3[comp_index] <- 1 # Mark detected compounds with 1

    # Update the dataframe with the new values
    df[i, start_col:ncol(df)] <- tmp3 * df[i, start_col:ncol(df)]
  }

  # Convert 0 values to NA for better data handling
  df[, start_col:ncol(df)][df[, start_col:ncol(df)] == 0] <- NA

  return(df)
}

# section -----------------------------------------------------------
#' Semi-quantification of Analytes
#'
#' This function performs semi-quantification of analytes in a given data frame considering signal
#' correction by Naphthalene-D8. It calculates concentrations for each batch separately. For
#' negative concentrations, it will then use the lowest concentration of the reference standard to
#' do single-point correction. If sample weights are provided, the concentrations will be expressed
#' in samples. Otherwise, the concentrations will be basedo on instrument-measured signals only.
#'
#' @param df A data frame containing peak areas and reference standards for each analyte.
#' @param sam_weight A data frame containing sample codes, weights, and final volumns.
#' @param start_col The starting column index of peak area columns (default is 12).
#' @return A data frame with updated concentrations.
#'
#' @import dplyr
#' @importFrom purrr map_dfc
#'
#' @export
calculate_con <- function(
    df,
    sam_weight,
    start_col = 12) {
  # Extract concentration columns from the data frame
  con_col <- df[, start_col:ncol(df)]

  # Identify unique batches from the column names
  batches_x <- unique(sub("_.*", "", colnames(con_col)))
  batches <- sub("^X", "", batches_x)
  num_batches <- length(batches)
  # Display the number of batches present
  cat("There are", num_batches, "batches:\n", batches, "\n")

  # Calculate concentrations for each batch separately
  for (i in seq_len(num_batches)) {
    d8_area <- df[
      grep("Naphthalene-D8", df$Name, ignore.case = TRUE),
      grep(paste(batches_x[i], "D8", sep = ".*"), colnames(df), ignore.case = TRUE)
    ]
    d8_area <- rowMeans(d8_area) # Calculate mean of d8 areas

    # Calculate concentrations using calibration curves
    for (j in seq_len(nrow(df))) {
      area <- df[j, grep(batches_x[i], colnames(df))]
      con <- (df[j, grep(batches_x[i], colnames(df))] / d8_area - df$Intercept[j]) / df$Slope[j]

      # Check for negative concentrations and recalculate using single point calibration
      for (k in seq_len(ncol(area))) {
        if (!is.na(con[1, k]) && con[1, k] < 0) {
          con[1, k] <- (area[1, k] / d8_area * df$Low_con[i]) / (df$Low_area[k] / d8_area) # nolint
        }
      }

      # Assign calculated concentrations back to the data frame
      df[j, grep(batches_x[i], colnames(df))] <- con
    }
  }

  # Remove rows and columns related to D8 or blank samples
  df <- df[
    !grepl("D8|bk", df$Name, ignore.case = TRUE),
    !grepl("D8|bk", colnames(df), ignore.case = TRUE)
  ]

  # Calculate concentrations in samples based on provided sample weights
  if (is.null(sam_weight)) {
    cat("No sample weights are provided, thus the concentrations is instrument-measured concentration but not sample concentrations.\n") # nolint
  } else {
    df_tmp <- purrr::map_dfc(
      sam_weight$Name,
      ~ {
        name_data <- dplyr::select(df, contains(.x))
        name_volumn <- sam_weight$Volumn[match(.x, sam_weight$Name)]
        name_weight <- sam_weight$Weight[match(.x, sam_weight$Name)]
        mutate(name_data, across(starts_with("X"), ~ (. * name_volumn) / name_weight))
      }
    )
    df <- dplyr::bind_cols(df[, 1:(start_col - 1)], df_tmp)
  }

  return(df)
}

# section -----------------------------------------------------------
#' Organize and summarize data after semi0quantification
#'
#' This function takes a dataframe and a sample code dataframe to organize and summarize data after
#' semi-quantification. It calculates minimum, median, maximum, and overall mean for each sample.
#'  Optionally, you can choose to replace NA values with zero or not. Moreover, you can also choose
#' wether to bind mean and standard deviation into a single column for each sample or not.
#'
#' @param df A dataframe after being processed by the 'calculate_con' function.
#' @param sam_code A dataframe containing sample codes and their corresponding names.
#' @param start_col The starting column index of peak area columns (default is 12).
#' @param digits An integer specifying the number of decimal places to round the results (default
#' is 2).
#' @param na2zero A logical value indicating whether to replace NA values with zero (default is
#' TRUE).
#' @param bind_mean_sd A logical value indicating whether to bind mean and standard deviation into
#' a single column for each sample (default is TRUE).
#'
#' @return A dataframe with summarized and organized data.
#'
#' @import dplyr
#' @importFrom purrr map_dfc
#' @importFrom stats median
#'
#' @export
organize_con <- function(
    df,
    sam_code,
    start_col = 12,
    digits = 2,
    na2zero = TRUE, # nolint
    bind_mean_sd = TRUE) {
  # Define the base dataframe excluding columns from start_col onwards
  df_base <- df[, 1:(start_col - 1)]
  # Extract sample names from the sample code dataframe
  samples <- sam_code$Name

  # Calculate mean for each sample group
  df_mean <- purrr::map_dfc(
    samples,
    ~ select(df, contains(.x)) %>%
      mutate(!!paste0(.x, "_mean") := round(rowMeans(.), digits = digits), .keep = "none")
  )

  # Calculate standard deviation for each sample group
  df_sd <- purrr::map_dfc(
    samples,
    ~ select(df, contains(.x) & !contains("_mean")) %>%
      mutate(!!paste0(.x, "_sd") := round(apply(., 1, sd), digits = digits), .keep = "none")
  )

  # Calculate minimum, median, and maximum values from the mean dataframe
  df_3m <- df_mean %>%
    mutate(Minimum = round(apply(., 1, min, na.rm = TRUE), digits = digits)) %>%
    mutate(Median = round(apply(., 1, median, na.rm = TRUE), digits = digits)) %>%
    mutate(Maximum = round(apply(., 1, max, na.rm = TRUE), digits = digits)) %>%
    select(-contains("_mean"))

  # Replace NA values with zero if specified
  if (na2zero) {
    df_mean <- as.data.frame(sapply(df_mean, function(x) replace(x, is.na(x), 0)))
  }

  # Calculate overall mean from the mean dataframe
  df_1m <- df_mean %>%
    mutate(Mean = round(apply(., 1, mean, na.rm = TRUE), digits = digits)) %>%
    select(-contains("_mean"))

  # Combine minimum, median, maximum, and overall mean dataframes
  df_4m <- dplyr::bind_cols(df_3m, df_1m)

  # Bind mean and standard deviation into a single column for each sample
  df_bind <- purrr::map_dfc(
    samples,
    ~ select(dplyr::bind_cols(df_mean, df_sd), contains(.x)) %>%
      mutate(
        !!.x := apply(., 1, function(x) {
          paste(format(na.omit(round(x, digits = digits)), nsmall = digits), collapse = " +/- ")
        }),
        .keep = "none"
      )
  )

  # Combine base dataframe with calculated statistics based on bind_mean_sd flag
  if (bind_mean_sd) {
    df <- dplyr::bind_cols(df_base, df_4m, df_bind)
  } else {
    df <- dplyr::bind_cols(df_base, df_4m, df_mean, df_sd)
  }

  # Replace zeros back to NA if any were introduced during processing
  df[df == 0] <- NA

  return(df)
}

# section -----------------------------------------------------------
#' Normalize Area by D8 Internal Standard
#'
#' This function normalizes the area of compounds in a data frame by the area of the internal
#' standard (D8). It processes each batch separately and removes rows and columns related to D8 or
#' blank samples.
#'
#' @param df A data frame containing the data to be normalized.
#' @param start_col The starting column index from which concentration columns begin. Default is 12.
#'
#' @return A data frame with normalized areas and without D8 or blank samples.
#' @export
normalize_area <- function(
    df,
    start_col = 12) {
  # Extract concentration columns from the data frame
  con_col <- df[, start_col:ncol(df)]

  # Identify unique batches from the column names
  batches_x <- unique(sub("_.*", "", colnames(con_col)))
  batches <- sub("^X", "", batches_x)
  num_batches <- length(batches)
  # Display the number of batches present
  cat("There are", num_batches, "batches:\n", batches, "\n")

  # Calculate concentrations for each batch separately
  for (i in seq_len(num_batches)) {
    d8_area <- df[
      grep("Naphthalene-D8", df$Name, ignore.case = TRUE),
      grep(paste(batches_x[i], "D8", sep = ".*"), colnames(df), ignore.case = TRUE)
    ]
    d8_area <- rowMeans(d8_area) # Calculate mean of d8 areas

    # Normalize the area to D8 for each batch
    df[, grep(batches_x[i], colnames(df))] <- df[, grep(batches_x[i], colnames(df))] / d8_area
  }

  # Remove rows and columns related to D8 or blank samples
  df <- df[
    !grepl("D8|bk", df$Name, ignore.case = TRUE),
    !grepl("D8|bk", colnames(df), ignore.case = TRUE)
  ]

  return(df)
}