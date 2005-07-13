\name{missinfo.sim}
\alias{missinfo.sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate fraction of missing information by simulation}
\description{Calculates fraction of missing information in a backcross
  by simulation.  
}
\usage{
missinfo.sim(delta, n, alpha, theta=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{QTL effect}
  \item{n}{Sample size}
  \item{alpha}{Selection fraction}
  \item{theta}{Recombination fraction(s) from flanking markers}
}
\details{Data is  simulated from a backcross with n individuals
  segregating one QTL.
}
\value{Fraction of missing information.
}
\references{Sen, Satagopan, and Churchill (2004), QTL study design from
  an information perspective,
  \url{http://repositories.cdlib.org/cbmb/QTLdesign}.
}
\author{Saunak Sen, Jaya Satagopan, and Gary Churchill}
\note{}

\seealso{\code{\link{missinfo}}}
\examples{
missinfo.sim(0,10000,0.5)
missinfo.sim(0,10000,0.5,0.1)
missinfo.sim(0,10000,0.5,c(0.1,0.1))
}
\keyword{design}
