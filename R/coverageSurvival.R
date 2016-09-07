#' Coverage Survival
#'
#' Find the diminishing percentage of cytosines with increasing sequencing depth
#' @param data A dataframe
#' @param cytosines Dataframe with cytosine counts. Can be generated using methylQC::cytosines()
#' @param chromosome Chromosome number. Defaults to all (all chromosomes)
#' @export
#' @return a dataframe
#' @examples
#' data <- loadData("methylome.CGmap")
#' cytosines <- cytosines("genome.fa")
#' coverage <- coverageSurvival(data, cytosines, chromosome = "chr1")
coverageSurvival <- function(data, cytosines, chromosome = "all") {
  # Get the total number of cytosines in each context for the chromosome we are looking at
  if(chromosome %in% c("lambda", "chrL", "L")) {
    nC <- data.frame(CG = 5925, CHG = 5385, CHH = 11503)
    d <- subset(data, data$chr == chromosome)
  } else if(chromosome == "all") {
    # take last row (sum of all cytosines)
    nC <- cytosines[nrow(cytosines),]
    d <- subset(data, !(data$chr %in% c("L", "chrL", "lambda")))
  } else {
    id <- as.numeric(unlist(strsplit(chromosome, "chr"))[[2]])
    nC <- cytosines[id,]
    d <- subset(data, chr == chromosome)
  }
  
  d <- subset(d, d$context != "--")

  # Count the coverage for each cytosine in each context
  d <- dplyr::group_by(d, depth, context)
  d <- dplyr::mutate(d, count = n())
  d <- dplyr::ungroup(d)
  
  # add info for all contexts combined
  temp <- dplyr::select(d, count, depth)
  temp <- unique(temp)
  temp <- dplyr::group_by(temp, depth)
  colnames(temp)[1] <- "n"
  temp <- dplyr::mutate(temp, count = sum(n), context = "all")
  temp <- dplyr::ungroup(temp)
  temp <- dplyr::select(temp, -n)
  temp <- unique(temp)
  temp$perc <- (temp$count / sum(nC)) * 100
  temp <- dplyr::select(temp, depth, context, perc)

  d$perc <- d$count / ifelse(d$context == "CG", nC$CG, ifelse(d$context == "CHG", nC$CHG, nC$CHH)) * 100
  d <- dplyr::select(d, depth, context, perc)
  d <- unique(d)
  
  # join with all context data
  d <- rbind(d, temp)
  d <- dplyr::arrange(d, depth)
  d <- convertSurvival(d)
  
  return(d)
}

# Convert counts of coverage levels (eg 1000 C sites with coverage == 12x coverage)
# to the sum of all sites wiht given coverage or higher (eg 3000 C sites with >= 12x coverage)
convertSurvival <- function(d) {
  temp <- tidyr::spread(d, context, perc)
  lim <- nrow(temp)
  for(i in 1:nrow(temp)) {
    temp[i,"totAll"] <- sum(temp[i:lim,2], na.rm = T)
    temp[i,"totCG"] <- sum(temp[i:lim,3], na.rm = T)
    temp[i,"totCHG"] <- sum(temp[i:lim,4], na.rm = T)
    temp[i,"totCHH"] <- sum(temp[i:lim,5], na.rm = T)
  }
  temp <- dplyr::select(temp, -(all:CHH))
  temp <- tidyr::gather(temp, Context, Cytosines, totAll:totCHH)
  temp$Context <- gsub("tot", "", temp$Context)
  return(temp)
}
