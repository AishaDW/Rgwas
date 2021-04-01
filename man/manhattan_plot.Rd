\name{manhattan_plot}
\alias{manhattan_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
manhattan_plot
}
\description{
Creates a manhattan plot
}
\usage{
manhattan_plot <- function(marker_map,p.vals,cutoff=NULL,trait = "unknown")
}

\arguments{
  \item{marker_map}{
genetic map with marker in first column, chromosome number and base pair position in subsequent columns
}
\item{p.vals}{P.vaues ffrom the GWAS analysis}
\item{cutoff}{Significance threshold. Bonferonni, unless otherwise specified}
\item{trait}{However you want to name the trait GWAS is being performed on}
}

\value{
The function will return a manhattan plot
}
\references{
Statistical genomics Spring 2020
}
\author{
Adapted from the Lab notes of Mathew McGowan
}


