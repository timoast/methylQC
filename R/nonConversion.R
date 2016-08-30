#' Non-conversion
#'
#' Calculate the bisulfite non-conversion rate
#' @param data A dataframe
#' @export
#' @return a dataframe
#' @examples
#' nonConversion(data)
nonConversion <- function(data) {
  d <- dplyr::filter(data, chr == "L")
  d <- dplyr::group_by(d, context)
  d <- dplyr::mutate(d, `Non-conversion rate` = sum(mC) / sum(depth) * 100)
  d <- dplyr::select(d, context, `Non-conversion rate`)
  d <- unique(d)
  d <- dplyr::filter(d, context %in% c("CG", "CHG", "CHH"))
  d <- dplyr::arrange(d, context)
  return(d)
}
