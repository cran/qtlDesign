\name{Utility}
\alias{recomb}
\alias{genetic.dist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Utility functions}
\description{Utility functions
}
\usage{
recomb(d)
genetic.dist(theta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{Genetic distance in Morgans}
  \item{theta}{Recombination fraction}
}
%\details{}
\value{\code{recomb} returns the recombination fraction
  corresponding to a genetic distance in Morgans.  \code{genetic.dist}
  returns the genetic distance in Morgans for a recombination fraction.
}

\references{Sen S, Satagopan JM, Churchill GA (2005) 	
Quantitative trait locus study design from an information perspective.
Genetics, 170:447-64.}
\author{Saunak Sen, Jaya Satagopan, Karl Broman, and Gary Churchill}
\note{We assume Haldane mapping function for the genetic distance.}

\examples{
recomb(0.1)
genetic.dist(0.1)
}
\keyword{utilities}
