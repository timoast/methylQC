#' Load Data
#'
#' This function will load the dataset to be used in methylQC functions
#' @param path Path to the data. Input can be gzipped or raw text.
#' @param datasource Which bisulfite package was used? Defaults to BSseeker2
#' @param forceZipped force data to be read as gzip. Default is to detect from the file extension.
#' @keywords data
#' @export
#' @return a datatable
#' @examples
#' loadData("data-raw/methylome.CGmap.gz")
loadData <- function(path, datasource="BSseeker2", forceZipped = FALSE) {
  header <- c("chr", "base", "position", "context", "twoBaseContext", "mC", "C", "depth")
  if(forceZipped == TRUE){
    dat <- data.table::fread(paste('gzip -dc ', path), header = TRUE, col.names = header)
    return(dat)
  }
  ext <- tools::file_ext(path)
  if(ext == "gz"){
    dat <- data.table::fread(paste('gzip -dc ', path), header = TRUE, col.names = header)
  } else {
    dat <- data.table::fread(path, sep = "\t", header = TRUE, col.names = header)
  }
  return(dat)
}