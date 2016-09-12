#' DNA methylation data
#'
#' A dataset containing Arabidopsis DNA methylation data for the first 1 million cytosines of chromosome 1
#' and the lambda spike-in.
#'
#' @format A data frame with 1024146 rows and 8 columns
#' \itemize{
#'   \item chr: chromosome
#'   \item base: the sequenced base
#'   \item position: chromosome position
#'   \item context: three-base cytosine sequence context
#'   \item twoBaseContext: dinucleotide sequence of the sequenced cytosine
#'   \item mC: percentage of reads that were methylated at the position
#'   \item C: number of sequenced reads at the position that were not methylated
#'   \item depth: total number of times the cytosine was sequenced
#' }
#' @source Dowen RH, Pelizzola M, Schmitz RJ, Lister R, Dowen JM, Nery JR, et al. Widespread dynamic DNA methylation in response to biotic stress. Proc Natl Acad Sci USA. National Acad Sciences; 2012;109: E2183–E2191. doi:10.1073/pnas.1209329109/-/DCSupplemental
"methylome"

#' Arabidopsis cytosine data
#'
#' A dataframe containing cytosine counts for each sequence context and chromosome
#'
#' @format A data frame with 8 rows and 3 columns
#' \itemize{
#'   \item CG: Number of cytosines in the CG context
#'   \item CHG: Number of cytosines in the CHG context
#'   \item CHH: Number of cytosines in the CHH context
#' }
#' @source \url{ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Arabidopsis_thaliana/Ensembl/TAIR10/Arabidopsis_thaliana_Ensembl_TAIR10.tar.gz}
"arabidopsis"