#' Define a function to extract classfication information
#'
#' @param x Objects from classifireR::classification
#'
#' @return A data.frame includes classification from classifireR
extract_cla <- function(x) {
  cla <- data.frame(matrix(ncol = nrow(x), nrow = 0)) # tibble does not work
  cla <- rbind(cla, x$Classification)
  colnames(cla) <- x$Level

  return(cla)
}
