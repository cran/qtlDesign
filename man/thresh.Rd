\name{Thresholds and tail probabilities}
\alias{tailprob}
\alias{thresh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculating thresholds and tail probabilities for genome scans}
\description{
Provides genome-wide thresholds and tail probabilities for the maxima of
genome scans using Poisson approximations.}
\usage{
tailprob(t,G,cross,type="1",d=0.01,cov.dim=0)
thresh(G,cross,type="1",p=c(0.10,0.05,0.01),d=0.01,cov.dim=0,
       interval=c(1,40))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{Genome size in centiMorgans.}
  \item{t}{LOD value for which tail probability is desired.}
  \item{p}{Vector giving the genome-wide Type I error for which
    thresholds are desired.}
  \item{cross}{String indicating cross type which is "bc", for 
    backcross, "f2" for intercross.}
  \item{type}{Type of LOD score for which threshold is desired.  Right
    now the only option is "1", but more options will be added in the
    future.}
  \item{d}{Marker spacing in centiMorgans.}
  \item{cov.dim}{Dimension of interacting covariate.  Set to 0 right now.}
  \item{interval}{Interval over which to search for LOD threshold.}
  }
\details{The tail probabilities are calculated using the method of
  Dupuis and Siegmund (1999).  The thresholds are calculated by solving
  the tail probability equation using \code{uniroot}.  At this time only
  one-dimensional thresholds are calculated, but this function will be
  extended in the future.
}

\value{The function \code{tailprob} returns the probability that the
  genome-wide maximum LOD score exceeds a particular value.  The
  function \code{thresh} returns genome-wide LOD thresholds
  corresponding to a particular Type I error rate.
}

\references{Dupuis J and Siegmund D (1999) Statistical methods for
  mapping quantitative trait loci froma dense set of markers. Genetics
  151:373-386.}
\author{Saunak Sen, Jaya Satagopan, Karl Broman, and Gary
  Churchill}
\seealso{\code{\link{uniroot}}.}
\examples{
tailprob(t=3,G=1440,cross="f2",d=10)
thresh(G=1440,cross="bc",d=10)
}
\keyword{design}
