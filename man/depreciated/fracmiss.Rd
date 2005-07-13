\name{fracmiss}
\alias{fracmiss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate fraction of missing information for single observation}
\description{
  Function to calculate the fraction of missing information for a single
  individual with an observed phenotype, and possibly flanking marker
  information.
}
\usage{
fracmiss(y, delta, m=NULL, theta=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{Value of phenotype}
  \item{delta}{QTL effect}
  \item{m}{Marker genotype(s)}
  \item{theta}{Recombination distance from flanking marker(s)}
}
\details{
}
\value{Fraction of missing information for a single observation given
  the observed value of the phenotype and neighbouring markers.
}
\references{Sen, Satagopan, and Churchill (2004), QTL study design from
  an information perspective,
  \url{http://repositories.cdlib.org/cbmb/QTLdesign}.
}
\author{Saunak Sen, Jaya Satagopan, and Gary Churchill}
\note{}
\seealso{\code{\link{missinfo.sim}}}
\examples{
fracmiss(1,0)
fracmiss(1,0,1,0.1)
fracmiss(1,0,c(1,1),c(0.1,0.1))
}
\keyword{design}
