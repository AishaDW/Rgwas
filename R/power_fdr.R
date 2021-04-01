#' This function was entirely recyled from a previous class. We use it only once not in the GWAS function of the package, but for a question
#' on the homework. We wanted to experiment different ways of calculating power than the GAPIT function we have used in our GWAS functions.
#' @param P.value vector with the p.value results from the GWAS.
#' @param QTN.position is the base pair position of QTN obtained the GWASbyCor function from zzlab.
#' @param cutoff is the threshold used to define significant marker.

power_fdr <- function(P.value, QTN.position=c(),cutoff=NULL) {
  pwr <- c()
  fdr <- c()
  t1error <- c()
  NQTN <- length(QTN.position)
  nsnp <- length(P.value)
  order.SNP=order(P.value)
  cutoff.final=(ifelse(
    is.null(cutoff),
    0.05/nsnp,
    cutoff
  ))
  sig.SNP <- order.SNP[sort(P.value)<=cutoff.final]
  TP=intersect(sig.SNP,QTN.position)
  FP=setdiff(sig.SNP, QTN.position)
  for (m in (1: nsnp)) {
    detected <- intersect(order.SNP[1:m], QTN.position)
    falsePositive <- setdiff(order.SNP[1:m], QTN.position)
    pwr <- c(pwr, length(detected)/NQTN)
    fdr <- c(fdr, length(falsePositive)/m)
    t1error <- c(t1error, length(falsePositive)/nsnp)
  }
  return(list(power=pwr, fdr=fdr, type1error=t1error,FP.fdr.power=FP,TP.fdr.power=TP))

}
