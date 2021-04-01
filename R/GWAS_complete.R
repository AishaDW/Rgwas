#' GWAS with PCA and linear dependency between covariates and PC accounted for
#' @param pheno file with numeric phenotypic values
#' @param geno data.frame with genotype calls coded as 0,1,2.
#' @param Cov numeric data.frame with covariates values
#' @param GM genetic map of data with chr and position of each SNP
#' @param cutoff  If cutoff is default, uses Bonferroni;0.05/number of SNPs
#' @return the GWAS results, the matrix with the principal component anlysis values and the power analysis results will be returned
#' in the global environment unless saved elsewhere by the user. The function will print the manhattan plot, qq plot and the correlation
#' between the three principal components.



GWAS_complete <- function(pheno=NULL, geno=NULL, Cov=NULL, GM=NULL,cutoff=NULL){
  G=pheno[,-1]
  GD=geno[,-1]
  n=nrow(GD)
  m=ncol(GD)
  CV=Cov[,-1]
  y=G
  PCA=prcomp(GD)


  #Variance Explained & Cumulative Variance Explained
  eig.val <- get_eigenvalue(PCA)
  print(head(eig.val))

  scree_plot <- print(fviz_eig(PCA, addlabels = TRUE, ylim = c(0, 10)))

  # Investigate PCA results
  PCA_plot_data <- data.frame(PCA$x)
  pca_comp_plot_12 <-
    ggplot(data = PCA_plot_data, aes(x = PC1, y = PC2)) +
    geom_point()
  pca_comp_plot_13 <-
    ggplot(data = PCA_plot_data, aes(x = PC1, y = PC3)) +
    geom_point()
  pca_comp_plot_23 <-
    ggplot(data = PCA_plot_data, aes(x = PC2, y = PC3)) +
    geom_point()
  Cor_plot <- print(grid.arrange(pca_comp_plot_12, pca_comp_plot_13, pca_comp_plot_23, nrow = 2, ncol = 2, top = "Principal Components"))

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

  ##### Minor allele frequency of associated SNPs ############
  minor_allele_freq <- apply(GD, 2, function(x) #Minor allele frequency calculation from previous class code
  {
    allele_freq1 <- (sum(x == 0)*2 + sum(x == 1)) / (sum(!is.na(x)) * 2)
    allele_freq2 <- (sum(x == 2)*2 + sum(x == 1)) / (sum(!is.na(x)) * 2)
    return(min(allele_freq1, allele_freq2))
  })
  sig.MAF<-minor_allele_freq[sig.SNP]

  # Return GWAS numeric results
  GWAS.Results <- list(P.value=P.value, cutoff.final=cutoff.final, sig.SNP=sig.SNP, sig.SNP.P=P.value[sig.SNP],
                       order.SNP=order.SNP, maf=sig.MAF)


  ####################### Create a manhattan plot ###################

  mh.plot <- print(manhattan_plot(marker_map,P.value,cutoff=NULL,trait = "Simulated Trait"))

  ############### Generate QQ plot #################################

  q.plot <- print(qq_plot(marker_map, P.value, trait = "Simulated Trait"))

  return(list(PCA.res=eig.val, Results=GWAS.Results))
}
