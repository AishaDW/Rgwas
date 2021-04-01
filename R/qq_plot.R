#' QQ plot plotted with p. values form the GWAS results and marker map
#' Important: This function was adapted from previous class
#' @param marker_map file with markers listed in the same order as the genotype file in the first column. The chromosome number in the second column
#' #' and the marker position in base pair in the thrird column
#' @param pvals vector with the p.value result from GWAS
#' @return QQ plot.



qq_plot <- function(marker_map, pvals, trait = "Simulated trait")
{
  marker_map$pvals <- -log10(t(pvals)) # Add pvalues to the data.frame and log10 transform

  exp_pval_dist <- -log10(runif(nrow(marker_map), 0, 1)) # Sample random p-values from a uniform distribution between 0 and 1

  qq_plot <- ggplot(marker_map, aes(x = sort(exp_pval_dist),
                                    y = sort(marker_map$pvals))) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    labs(title = "Q-Q Plot",
         x = "-log10(p-values) Expected",
         y = "-log10(p-values) Observed")

  return(qq_plot)
}
