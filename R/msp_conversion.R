#' Convert MSP file to dataframe
#'
#' This function reads an MSP (Mass Spectra Project) file and converts it into a dataframe format.
#'
#' @param msp_file Path to the MSP file to be converted
#' @return A dataframe containing the MSP data with columns for each field and a Spectrum column for peak data
#' @export
#' @examples
#' \dontrun{
#' df <- msp2df("example.msp")
#' }
msp2df <- function(msp_file) {
  # Read the MSP file
  lines <- readLines(msp_file)

  # Initialize variables
  compounds <- list()
  current_compound <- list()
  spectrum <- c()

  for (line in lines) {
    if (line == "") {
      # End of compound, add to list
      if (length(current_compound) > 0) {
        # Add Spectrum column if it exists
        if (length(spectrum) > 0) {
          current_compound$Spectrum <- paste(spectrum, collapse = " ")
        }
        compounds <- append(compounds, list(current_compound))
        current_compound <- list()
        spectrum <- c()
      }
    } else if (grepl("Num Peaks:", line, ignore.case = TRUE)) {
      # Handle "Num Peaks" field
      num_peaks <- sub("Num Peaks: ", "", line, ignore.case = TRUE)
      current_compound$`Num_peaks` <- num_peaks
    } else if (grepl(":", line)) {
      # Parse field
      parts <- strsplit(line, ": ")[[1]]
      field_name <- trimws(parts[1])
      field_value <- trimws(parts[2])
      current_compound[[field_name]] <- field_value
    } else if ("Num_peaks" %in% names(current_compound)) {
      # Concatenate lines after "Num Peaks" into Spectrum column
      peak_parts <- strsplit(line, "\\s+")[[1]]
      spectrum <- c(spectrum, paste(peak_parts[1], peak_parts[2], sep = ":"))
    }
  }

  # Convert list of compounds to data frame
  df <- do.call(dplyr::bind_rows, lapply(compounds, as.data.frame, stringsAsFactors = FALSE))
  names(df) <- sub("\\.+$", "", names(df))

  return(df)
}

#' Convert dataframe to MSP format
#'
#' This function converts a dataframe back to MSP (Mass Spectra Project) format and writes it to a file.
#'
#' @param df Dataframe containing the MSP data
#' @param output_file Path to the output MSP file
#' @export
#' @examples
#' \dontrun{
#' df2msp(my_data, "output.msp")
#' }
df2msp <- function(df, output_file) {
  # Open connection to output file with Windows-style line endings
  con <- file(output_file, "w", encoding = "UTF-8")

  # Process each compound
  for (i in seq_len(nrow(df))) {
    compound <- df[i, ]

    # Write all columns except Spectrum and Num_peak
    for (field in colnames(compound)) {
      if (field != "Spectrum" && field != "Num_peaks" && !is.na(compound[[field]])) {
      # if (!(field %in% c("Spectrum", "Num_peak", "Num Peaks")) && !is.na(compound[[field]])) {
        writeLines(paste0(field, ": ", compound[[field]]), con)
      }
    }

    # Handle Num Peaks
    if ("Spectrum" %in% colnames(compound) && !is.na(compound$Spectrum)) {
      # Calculate Num Peaks from Spectrum if Num_peak not available
      num_peaks <- length(strsplit(compound$Spectrum, " ")[[1]])
      writeLines(paste0("Num Peaks: ", num_peaks), con)
    }

    # Write spectrum data if present
    if ("Spectrum" %in% colnames(compound) && !is.na(compound$Spectrum)) {
      # Write each peak with tab separation
      peaks <- strsplit(compound$Spectrum, " ")[[1]]
      for (peak in peaks) {
        mz_intensity <- strsplit(peak, ":")[[1]]
        writeLines(paste(mz_intensity[1], mz_intensity[2], sep = "\t"), con)
      }
    }

    # Add two newlines between compounds
    writeLines("\r\n", con)
  }

  # Close connection
  close(con)
}