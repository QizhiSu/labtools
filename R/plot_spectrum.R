#' Convert spectral string to a named numeric intensity vector
#'
#' This function parses a space-separated string of m/z:intensity pairs (e.g., "100:30 120:90")
#' into a named numeric vector representing a binned and normalized mass spectrum.
#' It forms the foundation for downstream spectral similarity and visualization tasks.
#'
#' @param spectrum_str A character string representing the spectrum,
#' formatted as "m/z1:intensity1 m/z2:intensity2 ...".
#' @param start_mz Numeric. Starting m/z value of the output vector. Default is 50.
#' @param end_mz Numeric. Ending m/z value of the output vector. Default is 500.
#' @param mz_step Numeric. Step size (bin width) for m/z axis. Default is 1.
#' @param digits Integer. Number of decimal places to which m/z values are rounded. Default is 0.
#'
#' @return A named numeric vector of intensities, where names are m/z values (as characters).
#'
#' @examples
#' update_spectrum("98.7:10 99.9:20 101.2:15", start_mz = 98, end_mz = 102)
#'
#' @export
update_spectrum <- function(spectrum_str, start_mz = 50, end_mz = 500, mz_step = 1, digits = 0) {
  # Create empty spectrum vector
  mz_range <- seq(start_mz, end_mz, by = mz_step)
  spectrum <- numeric(length(mz_range))
  names(spectrum) <- as.character(mz_range)

  # Split spectrum string into m/z:intensity pairs
  values <- strsplit(spectrum_str, " ")[[1]]
  pairs <- strsplit(values, ":")

  # Extract m/z and intensity
  m_z <- sapply(pairs, function(x) as.numeric(x[1]))
  intensity <- sapply(pairs, function(x) as.numeric(x[2]))

  # Round m/z values
  m_z <- round(m_z, digits)

  # Filter by range
  in_range <- m_z >= start_mz & m_z <= end_mz
  m_z <- m_z[in_range]
  intensity <- intensity[in_range]

  # Remove duplicated m/z (keep first occurrence)
  unique_mz <- !duplicated(m_z)
  m_z <- m_z[unique_mz]
  intensity <- intensity[unique_mz]

  # Assign intensity values
  spectrum_match <- as.character(m_z)
  spectrum[spectrum_match] <- intensity

  return(spectrum)
}


#' Calculate cosine similarity between two spectra
#'
#' Computes the cosine similarity between two numeric spectra, typically mass spectra,
#' by matching shared m/z values and evaluating vector similarity.
#'
#' @param spec1 A named numeric vector representing spectrum 1 (m/z as names, intensity as values).
#' @param spec2 A named numeric vector representing spectrum 2.
#'
#' @return A numeric value between 0 and 1 representing the cosine similarity.
#' If there are no common m/z values, returns 0.
#'
#' @examples
#' spec1 <- update_spectrum("100:10 150:20 200:30")
#' spec2 <- update_spectrum("100:10 150:25 250:40")
#' cosine_similarity(spec1, spec2)
#'
#' @export
#'
cosine_similarity <- function(spec1, spec2) {
  common_mz <- intersect(names(spec1), names(spec2))
  if (length(common_mz) == 0) {
    return(0)
  } # No common peaks
  dot_product <- sum(spec1[common_mz] * spec2[common_mz])
  magnitude1 <- sqrt(sum(spec1[common_mz]^2))
  magnitude2 <- sqrt(sum(spec2[common_mz]^2))

  return(dot_product / (magnitude1 * magnitude2))
}


#' Plot a single mass spectrum
#'
#' Creates an interactive bar plot of a single spectrum using ggplot2 and plotly.
#' Peaks are normalized to 100% intensity. The highest peak in each m/z bin is optionally labeled.
#'
#' @param spectrum A named numeric vector representing the spectrum (m/z as names, intensity as values).
#' @param range Integer. The bin width (in m/z) used to identify local peak maxima for labeling. Default is 10.
#' @param threshold Numeric. The minimum intensity (in %) a peak must have to be labeled. Default is 1.
#' @param max_ticks Integer. Maximum number of x-axis tick marks. Default is 20.
#'
#' @return An interactive plotly plot showing the mass spectrum.
#'
#' @examples
#' spectrum <- update_spectrum("100:30 120:80 145:60 170:90 200:20")
#' plot_spectrum(spectrum)
#'
#' @importFrom ggplot2 ggplot aes geom_bar scale_y_continuous scale_x_continuous theme_classic theme element_line element_text geom_text
#' @importFrom dplyr mutate group_by ungroup filter select
#' @importFrom plotly ggplotly
#' @export
#'
plot_spectrum <- function(spectrum, range = 10, threshold = 1, max_ticks = 20) {
  # Convert spectrum to data frame
  spectrum_df <- data.frame(x = as.numeric(names(spectrum)), y = spectrum)

  # Normalize the y-values to the percentage of the highest peak
  max_y <- max(spectrum_df$y)
  spectrum_df <- spectrum_df %>%
    mutate(y = (y / max_y) * 100)

  # Identify the highest peak within each range of x ± 10
  spectrum_df <- spectrum_df %>%
    group_by(
      bin = cut(x,
        breaks = seq(floor(min(x) / range) * range, ceiling(max(x) / range) * range,
          by = range
        )
      )
    ) %>%
    mutate(is_highest_in_bin = (y == max(y))) %>%
    ungroup()

  # Define x-axis limits and breaks
  x_min <- floor(min(spectrum_df$x))
  x_max <- ceiling(max(spectrum_df$x))
  x_breaks <- seq(x_min, x_max, length.out = min(max_ticks, ceiling((x_max - x_min) / range) + 1))

  # Base ggplot with bars starting from the x-axis at 0 and axis label size 3
  plot <- ggplot(spectrum_df, aes(x, y)) +
    geom_bar(stat = "identity", fill = "#eb7771") + # Bar color set to #eb7771
    theme_classic() +
    theme(
      axis.line.y = element_line(color = "black"), # Retain y-axis line
      axis.text = element_text(size = 8), # Set axis text size
      axis.title = element_text(size = 8)
    ) + # Set axis title size
    scale_y_continuous(expand = c(0, 0), limits = c(0, 105), name = "Relative intensity (%)") +
    scale_x_continuous(expand = c(0, 0), breaks = round(x_breaks), name = "m/z")

  # Filter for the highest peaks in each bin to be labeled
  highest_peaks <- spectrum_df %>%
    filter(is_highest_in_bin, y >= threshold) %>%
    select(x, y)

  # Add annotations for the highest peaks within each x ± 10 range
  if (nrow(highest_peaks) > 0) {
    plot <- plot +
      geom_text(
        data = highest_peaks,
        aes(x = x, y = y, label = x), vjust = -0.5, color = "#080bf0", size = 3
      ) # Text color set to #080bf0, label size set to 3
  }

  # Convert to plotly for interactivity
  plotly_plot <- ggplotly(plot, tooltip = c("x", "y"))

  return(plotly_plot)
}


# section -----------------------------------------------------------
#' @importFrom ggplot2 aes geom_bar geom_hline scale_fill_manual theme_classic theme element_line element_text coord_cartesian geom_text
#' @importFrom dplyr bind_rows group_by mutate ungroup filter select
plot_mirrored_spectrum <- function(spec1, spec2, range = 20, threshold = 1, max_ticks = 20) {
  # Convert spectra to data frames
  spectrum_df1 <- data.frame(x = as.numeric(names(spec1)), y = spec1)
  spectrum_df2 <- data.frame(x = as.numeric(names(spec2)), y = -spec2)

  # Normalize the y-values to the percentage of the highest peak
  max_y1 <- max(spectrum_df1$y)
  max_y2 <- max(abs(spectrum_df2$y))
  spectrum_df1 <- spectrum_df1 %>%
    mutate(y = (y / max_y1) * 100)
  spectrum_df2 <- spectrum_df2 %>%
    mutate(y = (y / max_y2) * 100)

  # Combine the spectra
  combined_spectrum_df <- bind_rows(
    spectrum_df1 %>% mutate(spectrum = "spec1"),
    spectrum_df2 %>% mutate(spectrum = "spec2")
  )

  # Identify the highest peak within each range of x ± 10 for both spectra
  combined_spectrum_df <- combined_spectrum_df %>%
    group_by(
      spectrum,
      bin = cut(x, breaks = seq(floor(min(x) / range) * range,
        ceiling(max(x) / range) * range,
        by = range
      ))
    ) %>%
    mutate(is_highest_in_bin = (abs(y) == max(abs(y)))) %>%
    ungroup()

  # Define x-axis limits and breaks
  x_min <- floor(min(combined_spectrum_df$x))
  x_max <- ceiling(max(combined_spectrum_df$x))
  x_breaks <- seq(x_min, x_max, length.out = min(max_ticks, ceiling((x_max - x_min) / range) + 1))

  # Calculate cosine similarity
  cos_sim <- cosine_similarity(spec1, spec2)

  # Base ggplot with mirrored bars and axis label size 3
  plot <- ggplot(combined_spectrum_df, aes(x, y, fill = spectrum)) +
    geom_bar(stat = "identity", position = "identity", show.legend = FALSE) +
    geom_hline(yintercept = 0) +
    scale_fill_manual(values = c("spec1" = "blue", "spec2" = "red")) +
    scale_y_continuous(expand = c(0, 0), limits = c(-105, 105), name = "Relative intensity (%)") +
    scale_x_continuous(expand = c(0, 0), breaks = round(x_breaks), name = "m/z") +
    theme_classic() +
    theme(
      axis.line.y = element_line(color = "black"), # Retain y-axis line
      axis.text = element_text(size = 8), # Set axis text size
      axis.title = element_text(size = 8)
    ) + # Set axis title size
    coord_cartesian(clip = "off") # Allow text labels to be fully visible

  # Filter for the highest peaks in each bin to be labeled
  highest_peaks <- combined_spectrum_df %>%
    filter(is_highest_in_bin, abs(y) >= threshold) %>%
    select(x, y, spectrum)

  # Add annotations for the highest peaks in each bin
  if (nrow(highest_peaks) > 0) {
    plot <- plot +
      geom_text(
        data = highest_peaks %>% filter(spectrum == "spec1"),
        aes(x = x, y = y + 2, label = x), color = "#080bf0", size = 3
      ) +
      geom_text(
        data = highest_peaks %>% filter(spectrum == "spec2"),
        aes(x = x, y = y - 2, label = x), color = "red", size = 3
      )
  }

  # Convert to plotly for interactivity
  plotly_plot <- ggplotly(plot, tooltip = c("x", "y"))

  return(plotly_plot)
}
