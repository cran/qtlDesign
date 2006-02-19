\name{Information}
\alias{info}
\alias{info.bc}
\alias{info.f2}
\alias{deflate}
\alias{deflate.bc}
\alias{deflate.f2}
\alias{nullinfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Information under null hypothesis of equal means}
\description{
  Functions to calculate the information under the null hypothesis of
  no effect.  Functions for discount factors for incomplete genotyping.
}
\usage{
info(sel.frac,theta=0,cross)
info.bc(sel.frac,theta=0)
info.f2(sel.frac,theta=0)
deflate(theta,cross)
deflate.bc(theta)
deflate.f2(theta)
nullinfo(sel.frac)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cross}{Cross type, either "bc" for backcross, or "f2" for intercross.}
  \item{sel.frac}{Selection fraction; proportion of extremes genotyped}
  \item{theta}{Recombination fraction between flanking markers}
}
\details{The \code{nullinfo} function calculates the information
  content per observation for any contrast between genotype means when
  densely genotyping an \code{sel.frac} fraction of
  the extreme phenotypic individuals.  The information content is
  calculated under the null hypothesis of no difference between the
  genotype means.  For small differences in genotype means, the
  information content will be approximately equal to the null, but in
  general, the information estimate under the null is the lower bound.

  The \code{info} function calculates  the information per observation
  for  backcross, and F2 intercross  under the null hypothesis of equal
  gentoype means.  The   information is calculated for a point in the
  middle of an interval  spanned by markers separated by a recombination
  fraction \code{theta}.  The function \code{deflate} calculates a
  deflation factor for the information attenuation in the middle of a
  marker interval relative to a completely typed location.
}
\value{Information per individual for information functions, and the
  discount factor for the discount functions.
}
\references{Sen S, Satagopan JM, Churchill GA (2005) 	
Quantitative trait locus study design from an information perspective.
Genetics, 170:447-64.}
\author{Saunak Sen, Jaya Satagopan, Karl Broman, and Gary Churchill}
\note{Information is calculated under the equal means assumption.  This
  approximation is very good in practice, and is slightly conservative.
If the difference between the means is large, these functions will
underestimate the information.  For power calculations, that is okay.}
\examples{
nullinfo(0.5)
info(0.5,cross="bc")
info(0.5,cross="f2")
info(0.5,0.1,cross="bc")
info(0.5,0.1,cross="f2")
deflate(0.1,"bc")
deflate(0.1,"f2")
}
\keyword{design}
