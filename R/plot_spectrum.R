library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)

# section -----------------------------------------------------------
update_spectrum <- function(spectrum_str, start_mz = 50, end_mz = 500, mz_step = 1, digits = 0) {
  # create empty spectrum
  mz_range <- seq(start_mz, end_mz, by = mz_step)
  spectrum <- numeric(length(mz_range))
  names(spectrum) <- mz_range

  # update spectrum
  values <- str_split(spectrum_str, " ")[[1]]
  m_z_intensity <- do.call(rbind, str_split(values, ":"))
  m_z <- round(as.numeric(m_z_intensity[, 1]), digits = digits)
  intensity <- as.numeric(m_z_intensity[, 2])
  m_z_intensity <- data.frame(m_z, intensity) # convert to data frame
  m_z_intensity <- filter(m_z_intensity, m_z >= start_mz & m_z <= end_mz)
  m_z_intensity <- distinct(m_z_intensity, m_z, .keep_all = TRUE) # remove duplicated mz
  spectrum[names(spectrum) %in% m_z_intensity$m_z] <- m_z_intensity$intensity

  return(spectrum)
}

# section -----------------------------------------------------------
# Function to calculate cosine similarity between two spectra
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

# section -----------------------------------------------------------
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
