\name{GWAS_complete2}
\title{GWAS by GLM}
\description{
This function performs single marker GWAS uisng general linear model model. It also looks at the linear dependency between covariates and principal components and remove factors linearly dependent with each other.
}
\usage{
GWAS_complete <- function(pheno=NULL, geno=NULL, Cov=NULL, GM=NULL,cutoff=NULL)
}
\arguments{
  \item{pheno, geno, Cov, GM, cutoff}{
Pheno=phenotype data, geno=genotype data, Cov=Covariate data, GM=marker_map, cutoff=bonferonni threshold if no other specified
}
}
\details{ The manhattan plot and QQ plot are rreturned as PDF files in the working directory. Additionally in this function, the first column of the phenotype data is not removed because the simulated phenotype is used to run the function.

}
\value{
The function returns the P.value for each marker. It also returns the values from the PCA and the power analysis using the GAPIT function. Last a manhattan and QQ plot are saved in the working directory.
}
\references{
"http://www.zzlab.net/GAPIT/gapit_functions.txt"
}
\author{
Aichatou Djibo Waziri, Richard Manasseh
}
\note{
This function uses the fixDependence function from the mgvc package.
It also uses the GAPIT.FDR.TypeI function form the GAPIT package to compute power, FDR, AUC etc..
}

\examples{
GWAS_complete2 <- function(pheno=NULL, geno=NULL, Cov=NULL, GM=NULL,cutoff=NULL)
}
