\name{GWAS_simple}
\alias{GWAS_simple}
\title{Function to perform GWAS with P.value as only output}
\usage{
GWAS_simple <- function(pheno=NULL, geno=NULL, Cov=NULL, GM=NULL, cutoff=NULL)
}
\description{This function performs GWAS by general linear model including three principal componens. The linear dependency between the covariates and the principal components is investigated and redundant factors are eliminated. A list object called GWAS.Results containing P.values is returned in the working environment.

GWAS.Results=list(P.value=P.value, cutoff.final=cutoff.final, sig.SNP=sig.SNP, sig.SNP=P.value[sig.SNP], order.SNP=order.SNP)
}
\arguments{
\item{pheno}{phenotypic data with first column containing individuals}

\item{geno}{genotypic data in a numeric format and the first column containing the individuals}

\item{Cov}{Covariate data with first column containing individuals}
\item{cutoff}{If cutoff is default, uses Bonferroni;0.05/number of SNPs}
}
