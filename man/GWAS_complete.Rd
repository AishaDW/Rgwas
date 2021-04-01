\name{GWAS_complete}
\alias{GWAS_complete}

\title{GWAS_complete
}
\description{
This function performs single marker GWAS uisng general linear model model. It also looks at the linear dependency between covariates and principal components and remove factors linearly dependent with each other. The minor allele frequency of associated SNPs are calculated and returned in the list of results.
}
\usage{
GWAS_complete <- function(pheno=NULL, geno=NULL, Cov=NULL, GM=NULL,cutoff=NULL)
}

\arguments{
  \item{pheno, geno, Cov, GM, cutoff}{
Pheno=phenotype data, geno=genotype data, Cov=Covariate data, GM=marker_map, cutoff=bonferonni threshold if no other specified
}
}
\details{
The details of data format could be retrieved from the GWAS_simple function description from the same package. For the functions in this package to work conveniently, it is important to import the datasets as data tables with the column containing the individuals active because the core functions remove that first column.
}
\value{
The function returns the P.value for each marker. Last a manhattan and QQ plot are printed. We also have the correlation plots between the different PCs.
}
\references{
"http://www.zzlab.net/GAPIT/gapit_functions.txt"
}
\author{
Aichatou Djibo Waziri, Richard Manasseh (Adapted from function from previous classes)
}
\note{
This function uses the fixDependence function from the mgvc package.
It also uses the GAPIT.FDR.TypeI function form the GAPIT package to compute power, FDR, AUC etc..
}

\examples{
GWAS_complete <- function(pheno=NULL, geno=NULL, Cov=NULL, GM=NULL,cutoff=NULL)
}

