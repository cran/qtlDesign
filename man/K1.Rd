\name{K1}
\alias{K1}
\alias{K12}
\alias{Kstat}
\alias{eff1}
\title{Calculate scores for minimum moment abberations.}
\description{
Calculate the MMA K1, K12, and the standardized dissimilarity score (eff1).
}
\usage{
Kstat(genomat, type = 1)
K1(genomat)
K12(genomat)
eff1(n, nmark, s1)
}
\arguments{
  \item{genomat}{Genotype matrix.}
  \item{n}{Desired sample size.}
  \item{type}{Type of dissimilarity measure desired (first or second
    moment).}
  \item{nmark}{Number of markers.}  
  \item{s1}{Dissimilarity score from \code{K1} or \code{K12}.}
}
\value{
  Score or standardized score based on selected marker list.
  \code{K1} and \code{K12} call \code{Kstat} with \code{type} = 1 and 2,
  respectively. Kstat computes the minimum moment abberation
  score. \code{eff1} computes the standardized genetic dissimilarity.
}
\references{Jin C, Lan H, Attie AD, Churchill GA, Bulutuglo D, Yandell BS (2004) Selective phenotyping for increased efficiency in genetic mapping studies. Genetics 168: 2285-2293.}
\author{Brian S. Yandell (\url{mailto:byandell@wisc.edu})}
\seealso{\code{\link{mma}}, \code{\link[qtl]{read.cross}}}
%\examples{}
\keyword{design}% at least one, from doc/KEYWORDS

