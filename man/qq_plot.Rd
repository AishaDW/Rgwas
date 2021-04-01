\name{qq_plot}
\alias{qq_plot}

\title{
qq-plot
}
\description{
Creates a quantile-quantile plot
}
\usage{
qq_plot <- function(marker_map, pvals, trait = "unknown")
}

\arguments{
  \item{marker_map}{
genetic map with marker in first column, chromosome number and base pair position in subsequent columns
}
\item{p.vals}{P.vaues ffrom the GWAS analysis}
\item{cutoff}{Significance threshold. Bonferonni, unless otherwise specified}
\item{trait}{However you want to name the trait GWAS is being performed on}
}


\value{ The function will return a QQ plot

}
\references{
Statistical genomics Spring 2020
}
\author{
Adapted from the Lab notes of Mathew McGowan
}
