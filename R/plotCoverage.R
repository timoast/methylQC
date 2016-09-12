#' Plot coverage
#'
#' Plot the diminishing percentage of cytosines with increasing sequencing depth
#' @param data A dataframe. Can be generated using methylQC::coverageSurvival()
#' @export
#' @return a ggplot2 object
#' @examples
#' coverage <- coverageSurvival(data = methylome, cytosines = arabidopsis, chromosome = "L")
#' plotCoverage(coverage)
plotCoverage <- function(data) {
  cg <- '#B4B464'
  chg <- '#6665AD'
  chh <- '#B29492'
  
  p <- ggplot2::ggplot(data, aes(depth, Cytosines, color = Context)) +
    geom_line() + geom_point() + scale_color_manual(values = c("black", cg, chg, chh)) +
    theme_bw() + ylab("Percentage of Cytosines") +
    ylim(c(0, 100)) + xlim(c(0, 100))
  return(p)
}
