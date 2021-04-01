\name{power_fdr}
\title{Function to caculate power, FDR and type-1 error, False Positives, and True Positives}
\usage{
power_fdr(P.value, QTN.position = NULL, cutoff = NULL)
}
\arguments{
\item{QTN.position}{position of QTN if known}

\item{cutoff}{If cutoff is default, uses Bonferroni;0.05/number of SNPs}

\item{P}{list of SNPs order by ascending  p-value}
}
\value{
list of power, FDR, type-1 error, False Positives, and True Positives
}
\description{ Similar to the GAPIT.FDR function
Function to caculate power, FDR and type-1 error, False Positives, and True Positives
}
\note{We did not write this function. It has been adapted from previous class}
