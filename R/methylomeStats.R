#' Methylome statistics
#'
#' Generate summary statistics for methylome data
#' @param data A dataframe
#' @export
#' @return a dataframe
#' @examples
#' methylomeStats(methylome)
methylomeStats <- function(data) {
  s <- broom::tidy(summary(data$depth))
  k <- moments::kurtosis(data$depth)
  return(cbind(s, data.frame(stringsAsFactors = FALSE, kurtosis = k)))
}