#' Manhattan plot plotted with p. values form the GWAS results and marker map
#' Important: This function was adapted from previous class
#' @param marker_map file with markers listed in the same order as the genotype file in the first column. The chromosome number in the second column
#' #' and the marker position in base pair in the thrird column
#' @param p.vals vector with the p.value result from GWAS
#' @param cutoff  If cutoff is default, uses Bonferroni;0.05/number of SNPs
#' @return Manhattan plot with a horizontal line representing the significance threshold.



manhattan_plot <- function(marker_map,p.vals,cutoff=NULL,trait = "unknown")
{
  m=nrow(marker_map)
  cutoff.final=(ifelse(
    is.null(cutoff),
    0.05/m,
    cutoff
  ))

  thresh <- -log10(cutoff.final)
  marker_map$pvals <- -log10(t(p.vals)) # Add pvalues to the data.frame and log10 transform
  marker_map$comb_pos <- marker_map$Chromosome * 1e9 + marker_map$Position

  manhattan_plot <- ggplot(marker_map, aes(x = 1:nrow(marker_map), y = pvals, color = factor(Chromosome))) +
    geom_point() +
    geom_hline(yintercept= thresh,color="magenta")+
    labs(title = paste("GWAS manhattan plot for trait:", as.character(trait)),
         y = "-log10(p-value)",
         x = "Marker Position",
         color = "Chromosome")

  return(manhattan_plot)

}
