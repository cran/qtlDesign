\name{Information-cost functions}
\alias{info2cost}
\alias{info2cost.bc}
\alias{info2cost.f2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Functions to calculate information-cost ratios}
\description{
Functions to calculate information cost-ratios.
}
\usage{
info2cost(sel.frac,cost,d,G=NULL,cross)
info2cost.bc(sel.frac,cost,d,G=NULL)
info2cost.f2(sel.frac,cost,d,G=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sel.frac}{Selection fraction; proportion of individuals genotyped}
  \item{cost}{Genotyping cost in units of raising an individual.  When
  \code{d=0} (dense genotyping), it is the cost of genotyping an
  individual.  When \code{d!=0}, it is the cost of a single marker
  genotype in an individual.}
  \item{d}{Marker spacing in centiMorgans}
  \item{G}{Genome size in Morgans}
  \item{cross}{Cross type, "bc"or "f2"} 
}
\details{The information calculations are done under the null hypothesis
  of no QTL effect.}

\value{For \code{d!=0} it calculates the ratio of information in the
  middle of a marker interval of length \code{d} cM to the cost of
  genotyping the cross.  For \code{d=0}, it calculates the ratio of
  information at any locus to the cost of genotyping the cross.}
\references{Sen S, Satagopan JM, Churchill GA (2005) 	
Quantitative trait locus study design from an information perspective.
Genetics, 170:447-64.}
\author{Saunak Sen, Jaya Satagopan, Karl Broman, and Gary Churchill}
\seealso{\code{\link{info}}}
\examples{
info2cost.bc(0.5,1)
info2cost.bc(0.5,1,10,1450)
}
\keyword{design}
