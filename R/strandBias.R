#' Strand sequencing bias
#'
#' This function will calculate the strand depth bias for each C base position
#' @param data input dataframe
#' @param window size of sliding window. Default is 30 bases.
#' @param fill Value used to fill in where there isn't enough data for window size. Default is 0.
#' @export
#' @return a dataframe
#' @examples
#' strandBias(methylome)
strandBias <- function(data, window = 30, fill = 0) {
  data$depth <- ifelse(data$base == "G", -(data$depth), data$depth)
  
  # split into different dataframes for each strand
  wat <- dplyr::filter(data, base == "C")
  cri <- dplyr::filter(data, base == "G")
  
  # fill in position information
  filled_wat <- fillBases(wat)
  filled_cri <- fillBases(cri)

  # calculate rolling mean for every position
  filled_wat$meanCovC <- zoo::rollapply(filled_wat$depth,
                                        width = window, mean,
                                        fill = fill, na.rm = TRUE)
  filled_wat <- dplyr::select(filled_wat, position, meanCovC)
  filled_cri$meanCovG <- zoo::rollapply(filled_cri$depth,
                                        width = window, mean,
                                        fill = fill, na.rm = TRUE)
  filled_cri <- dplyr::select(filled_cri, position, meanCovG)
  
  # Convert NaN to NA
  filled_wat$meanCovC[is.nan(filled_wat$meanCovC)] <- NA
  filled_cri$meanCovG[is.nan(filled_cri$meanCovG)] <- NA
  
  # intersect with original data
  data <- dplyr::left_join(data, filled_wat, by = "position")
  data <- dplyr::left_join(data, filled_cri, by = "position")

  # calculate bias at each position
  data <- dplyr::rowwise(data)
  data <- dplyr::mutate(data, strandBias = (meanCovC + meanCovG))
  data <- dplyr::select(data, -(meanCovC:meanCovG))

  return(data)
}

fillBases <- function(d) {
  mn <- min(d$position, na.rm = TRUE)
  mx <- max(d$position, na.rm = TRUE)
  newdf <- data.frame(position = seq(from = mn, to = mx),
                      stringsAsFactors = FALSE)
  n <- dplyr::full_join(d, newdf, by = "position")
  n <- dplyr::arrange(n, position)
  return(n)
}