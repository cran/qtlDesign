\name{Information-cost ratios}
\alias{info2cost.bc.null}
\alias{optalpha.bc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Functions to calculate information-cost ratios, and to find
  optimal selection fractions
}
\description{
Functions to calculate information cost-ratios under the null
hypothesis, and to find optimal selection fractions given cost.
}
\usage{
info2cost.bc.null(alpha, cost, d=0, G=NULL)
optalpha.bc(cost,d=0,G=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{Selection fraction; proportion of individuals genotyped}
  \item{cost}{Cost in units of raising individual}
  \item{d}{Marker spacing in Morgans}
  \item{G}{Genome size in Morgans}
}
\details{
}
\value{Either the information to cost ratio, or the optimal selection fraction.
}
\references{Sen, Satagopan, and Churchill (2004), QTL study design from
  an information perspective,
  \url{http://repositories.cdlib.org/cbmb/QTLdesign}.
}
\author{Saunak Sen, Jaya Satagopan, and Gary Churchill}

\note{Functions for backcross only right now.}

\seealso{\code{\link{info.bc.null}}}
\examples{
info2cost.bc.null(0.5,1)
info2cost.bc.null(0.5,1,0.1,14.50)
optalpha.bc(1)
optalpha.bc(0.001,0.1,14.50)
}
\keyword{design}
