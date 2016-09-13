#' Plot Browser
#'
#' This function will plot a browser view of sequencing coverage
#' @param data A dataframe
#' @param start The start position. Defaults to 5000
#' @param stop The stop position. Defaults to 15000
#' @export
#' @import ggplot2
#' @return ggplot2 object
#' @examples
#' plotBrowser(methylome)
plotBrowser <- function(data, start=5000, stop=15000) {
  # filter to get the specified region (assuming single chromosome)
  d <- dplyr::filter_(data, ~position > start, ~position < stop)
  
  # convert minus strand to negative depth values
  d$depth <- ifelse(d$base == "G", -(d$depth), d$depth)
  
  # find medians
  pos_median <- median(subset(d, d$base == "C")$depth)
  neg_median <- median(subset(d, d$base == "G")$depth)
  
  color <- "#4292C6"
  # make plots
  p <- ggplot(d, aes(position, depth, base)) +
    geom_hline(yintercept = 0, color="grey") +
    geom_hline(yintercept = c(pos_median, neg_median),
               color="grey", linetype=2) +
    geom_line(aes(color = base)) + theme_classic() +
    scale_color_manual(values = c(color, color)) +
    geom_smooth(se = FALSE, span = 1/10, method = "loess",
                color="red", size=1/2) +
    theme(legend.position = "none") +
    ggtitle("Coverage")
  
  return(p)
}