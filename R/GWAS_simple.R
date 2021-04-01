#' GWAS with number of PCs set to three and covariate data added.
#' @param pheno file with numeric phenotypic values
#' @param geno data.frame with genotype calls coded as 0,1,2.
#' @param Cov numeric data.frame with covariates values
#' @param GM genetic map of data with chr and position of each SNP
#' @return Only p-value results from GWAS run.


GWAS_simple <- function(pheno=NULL, geno=NULL, Cov=NULL, cutoff=NULL){
  G=pheno[,-1]
  GD=geno[,-1]
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

  P=t(matrix(P))
  P.value=P
  order.SNP=order(P.value)
  #If cutoff is default, uses Bonferroni
  cutoff.final=(ifelse(
    is.null(cutoff),
    0.05/m,
    cutoff
  ))
  sig.SNP <- order.SNP[sort(P.value)<=cutoff.final]

  GWAS.Results=list(P.value=P.value, cutoff.final=cutoff.final, sig.SNP=sig.SNP, sig.SNP.P=P.value[sig.SNP], order.SNP=order.SNP)
  return(GWAS.Results)
}
