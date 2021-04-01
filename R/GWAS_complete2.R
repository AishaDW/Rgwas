#' GWAS with PCA and linear dependency between covariates and PC accounted for. Unlike the function GWAS_complete from the same package
#' this function does not print the correlation plot between PCs, nor the manhattan and QQ plots. Those plots will be saved in whatever
#' working directory had been specified by the user.
#' @param pheno file with numeric phenotypic values
#' @param geno data.frame with genotype calls coded as 0,1,2.
#' @param Cov numeric data.frame with covariates values
#' @param GM genetic map of data with chr and position of each SNP
#' @param cutoff  If cutoff is default, uses Bonferroni;0.05/number of SNPs
#' @return the GWAS results and power analysis results.


GWAS_complete2 <- function(pheno=NULL, geno=NULL, Cov=NULL, QTN.position=c(), GM=NULL,cutoff=NULL){
  G=pheno
  GD=genotypes[,-1]
  n=nrow(GD)
  m=ncol(GD)
  CV=Cov[,-1]
  y=G
  PCA=prcomp(GD)


  ### Perform GWAS looping through each marker
  P=apply(GD,2, function(x)
    # CHECK FOR GENOTYPE DISTRIBUTION
    if(max(x)==min(x)){p=1
    P=p[length(p)]
    # IF MAX NOT EQUAL TO MIN
    }else{
      # CHECK FOR DEPENDENCE
      # IF NO CV INPUT
      if (is.null(CV)){X=cbind(1, PCA$x[,1:3],x)
      # IF THERE IS CV INPUT
      }else{fD <- fixDependence(PCA$x[,1:3], as.matrix(CV), tol=.Machine$double.eps^.5,rank.def=0,strict=FALSE)
      # WITH CV INPUT, CONDITION 1: NO DEPENDENCE
      if (is.null(fD)){X=cbind(1, PCA$x[,1:3],as.matrix(CV), x)

      # WITH CV INPUT, CONDITION 2: WITH DEPENDENCE
      }else {X=cbind(1, PCA$x[,1:3],x)}
      }# END FOR DEPENDECE
      # SOLVE THE LINEAR REGRESSION MATRIX:
      #  X=as.matrix(X)
      LHS=t(X)%*%X
      C=solve(LHS)
      RHS=t(X)%*%y
      b=C%*%RHS
      yb=X%*%b
      e=y-yb
      n=length(y)
      ve=sum(e^2)/(n-1)
      vt=C*ve
      t=b/sqrt(diag(vt))
      p=2*(1-pt(abs(t),n-2))
      P=p[length(p)]
    }
  )# End of loop

  P=t(as.numeric(matrix(P)))
  P.value=P
  order.SNP=order(P.value)
  #If cutoff is default, uses Bonferroni
  cutoff=NULL
  cutoff.final=(ifelse(
    is.null(cutoff),
    0.05/m,
    cutoff
  ))
  sig.SNP <- order.SNP[sort(P.value)<=cutoff.final]
  ######Power Analysis ########################################
  TP=intersect(sig.SNP,QTN.position)
  FP=setdiff(sig.SNP, QTN.position)

  bigNum=1e9
  resolution=100000
  bin=round((marker_map[,2]*bigNum+marker_map[,3])/resolution)
  result=cbind(marker_map,P.value,bin)

  QTN.bin=result[sig.SNP,]

  index.qtn.p=order(QTN.bin[,4])

  myStat=GAPIT.FDR.TypeI(
    # WS- window size
    WS=c(1e0,1e3,1e4,1e5),
    #Input: GM - m by 3  matrix for SNP name, chromosome and BP
    GM=marker_map,
    #Input: seqQTN - s by 1 vector for index of QTN on GM (+1 for GDP column wise)
    seqQTN=QTN.position,
    #Input: GWAS - SNP,CHR,BP,P,MAF
    GWAS=result)
  ############## PLots #########################################

  mh.plot <- manhattan_plot(marker_map,P.value,cutoff=NULL,trait = "Simulated")
  ggsave("mh.plot.pdf")
  ############### Generate QQ plot #################################

  q.plot <- qq_plot(marker_map, P.value, trait = "Simulated Trait")
  ggsave("q.plot.pdf")

  return(list(P.value=P.value, cutoff.final=cutoff.final, sig.SNP=sig.SNP, sig.SNP.P=P.value[sig.SNP], order.SNP=order.SNP, TP=TP, FP=FP,pwr.tst=myStat))

}
