\name{Optimal selection fraction}
\alias{optalpha}
\alias{optalpha.bc}
\alias{optalpha.f2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Optimal selection fraction}
\description{
Functions to find optimal selection fractions given cost.
}
\usage{
optalpha(cost,d=0,G=NULL,cross)
optalpha.bc(cost,d=0,G=NULL)
optalpha.f2(cost,d=0,G=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cost}{Cost per genotype in units of raising individual}
  \item{d}{Marker spacing in Morgans}
  \item{G}{Genome size in Morgans}
  \item{cross}{Cross type, "bc" or "f2"}
}
\details{The function \code{optimize} is used to search for the optima.}

\value{The optimal selection  fraction.}
\references{Sen S, Satagopan JM, Churchill GA (2005) 	
Quantitative trait locus study design from an information perspective.
Genetics, 170:447-64.}
\author{Saunak Sen, Jaya Satagopan, Karl Broman, and Gary Churchill}
\seealso{\code{\link{optimize}}}
\examples{
optalpha(1,cross="bc")
optalpha(0.001,10,1450,cross="bc")
optalpha(0.001,10,1450,cross="f2")
}
\keyword{design}
