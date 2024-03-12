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
#' @param std_res_col The column index of the response variable in std_md
#' @param std_md1_col The column index where the molecular descriptor starts
#' in std_md
#' @param data_md The data to be determined including at least the Name and
#' molecular descriptor for each standard. The molecular descriptor data should
#' be at the end of the dataframe
#' @param data_md1_col The column index where the molecular descriptor data
#' starts in data_md
#' @param top_k_md The number of top molecular descriptors to select (default
#' is 20)
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
    top_k_md = 20) {
  # construct the data for feature selection
  std4rf <- std_md[, c(std_res_col, std_md1_col:ncol(std_md))]
  colnames(std4rf)[1] <- "response" # rename the response column for later use
  # remove empty molecular descriptor
  std4rf <- std4rf[!apply(std4rf[, -1], 1, function(x) all(is.na(x))), ]
  std4rf <- std4rf[, !apply(std4rf, 2, function(x) any(is.na(x)))]

  # build a random forest model for calculation of importance metric
  rf_model <- randomForest::randomForest(response ~ .,
    data = std4rf,
    importance = TRUE,
    na.action = na.omit
  )

  # calculate the importance metric for each molecular descriptor
  vip_df <- randomForest::importance(rf_model) |> as.data.frame()
  vip_df <- vip_df[order(-vip_df$`%IncMSE`), ] # sort by decreasing importance
  vip_df <- vip_df[1:top_k_md, ] # select top k molecular descriptors
  vip_df$MD <- rownames(vip_df)

  # keep only top k molecular descriptors in the data and standard
  data_md <- data_md[, c(1:data_md1_col, which(names(data_md) %in% vip_df$MD))]
  std_md <- std_md[, c(1:std_md1_col, which(names(std_md) %in% vip_df$MD))]

  # create a metadata dataframe for the selected standard
  data_meta <- data_md[, c(1:data_md1_col)]
  data_meta$Standard <- NA
  data_meta$Similarity <- NA

  cat("Assigning standard for target compounds...\n")
  # calculate the similarity between the selected standard and each substance
  for (i in seq_len(nrow(data_md))) {
    cat(i, "out of", nrow(data_md), "\n")
    tmp <- vector() # create a temporary vector to store tempory similarity

    # calculate the similarity one by one
    for (j in seq_len(nrow(std_md))) {
      a <- data_md[i, data_md1_col:ncol(data_md)]
      b <- std_md[j, std_md1_col:ncol(std_md)]
      tmp[j] <- proxy::simil(a, b) # calculate the similarity using the simil
      data_meta$Standard[i] <- std_md$Name[which.max(tmp)]
      data_meta$Similarity[i] <- max(tmp)
    }
  }

  cat("Standard assignment completed.\n")
  return(data_meta)
}
