#' Cytosines
#'
#' Count the number of cytosines in each sequence context, for each chromosome
#' @param data Path to fasta file
#' @export
#' @return a dataframe
#' @examples
#' cytosines("./genome.fa")
cytosines <- function(data) {
  
  # regex for each context, each strand
  cg <- '(CG.)|(.CG)'
  chg <- '(C[ATC]G)|(G[ATC]C)'
  chh <- '(C[ATC]{2})|([ATC]{2}G)'
  
  # read fasta file
  fasta <- Biostrings::readDNAStringSet(data)
  
  # get frequency of each trinucleotide
  trinuc <- Biostrings::trinucleotideFrequency(fasta, step = 1)
  
  # convert to dataframe
  trinuc <- as.data.frame(trinuc)
  
  # subset to include only Cs
  cg_info <- trinuc[, grep(cg, names(trinuc))]
  chg_info <- trinuc[, grep(chg, names(trinuc))]
  chh_info <- trinuc[, grep(chh, names(trinuc))]
  
  # sum
  sum_cg <- rowSums(cg_info)
  sum_chg <- rowSums(chg_info)
  sum_chh <- rowSums(chh_info)
  all <- data.frame(cg = sum_cg, chg = sum_chg, chh = sum_chh)
  all <- rbind(all, colSums(all))
  
  return(all)
}
