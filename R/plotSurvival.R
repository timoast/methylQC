

convertSurvival <- function(d) {
  for(i in 1:nrow(d)) {
    tot <- sum(d[i:nrow(d),2])
    d[i,3] <- tot
  }
  return(d)
}

getChrom <- function(chr) {
  if(chr == "all"){
    return(1)
  } else {
    digit <- as.numeric(unlist(strsplit(chr, "chr"))[[2]])
    return(digit + 1)
  }
}

#' Plot Survival
#'
#' Plot the diminishing percentage of cytosines with increasing sequencing depth
#' @param data A dataframe
#' @param species Species name. Defaults to Arabidopsis
#' @param chromosome Chromosome number. Defaults to all (all chromosomes)
#' @export
#' @return ggplot2 plot object
#' @examples
#' plotBrowser(dat)
plotSurvival <- function(data, species="Arabidopsis", chromosome = "all") {
  library(dplyr)
  genomes <- data.frame(stringsAsFactors = F,
                        arabidopsis = c(42859753, 10856525, 7063739, 8521037, 6727440, 9691012),
                        human = 0,
                        brachypodium = 125937308,
                        brachy = 125937308,
                        mouse = 0)
  chrom_id <- getChrom(chromosome)
  species <- tolower(species)
  if(chromosome == "all") {
    d <- data
  } else {
    d <- subset(data, chr == chromosome)
  }
  id <- getChrom(chromosome)
  nC <- genomes[[species]][[id]]
  d <- data %>%
    group_by(depth) %>%
    mutate(count = n(), perc = count / nC * 100) %>%
    select(depth, perc) %>%
    unique() %>%
    arrange(depth) %>%
    convertSurvival() %>%
    select(-perc)
  colnames(d)[2] <- "Cytosines"
  p <- ggplot2::ggplot(d, aes(depth, Cytosines)) +
    geom_line() + geom_point() +
    theme_bw() + ylab("Percentage of Cytosines") +
    ylim(c(0, 100)) + xlim(c(0, 100))
  return(p)
}