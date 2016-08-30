#' Load Data
#'
#' This function will load the dataset to be used in methylQC functions
#' @param path Path to the data
#' @param datasource Which bisulfite package was used? Defaults to BSseeker2
#' @keywords data
#' @export
#' @return a dataframe
#' @examples
#' loadData("~/Documents/mC_data.CGmap.gz")

loadData <- function(path, datasource="BSseeker2") {
  if(datasource == "BSseeker2") {
    header <- c("chr", "base", "position", "context", "twoBaseContext", "mC", "C", "depth")
  } else {
    return()
  }
  dat <- readr::read_tsv(path, col_names = header)
  return(dat)
}