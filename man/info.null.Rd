\name{Information under null}
\alias{info.null}
\alias{info.bc.null}
\alias{info.f2.null}
\alias{deflate.bc}
\alias{deflate.f2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Information under null hypothesis of equal means}
\description{
  Functions to calculate the information under the null hypothesis of
  no effect.  Functions for discount factors for incomplete genotyping.
}
\usage{
info.null(alpha)
info.bc.null(alpha,theta=0)
info.f2.null(alpha,theta=0)
deflate.bc(theta)
deflate.f2(theta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{Selection fraction; proportion of extremes genotyped}
  \item{theta}{Recombination fraction between flanking markers}
}
\details{The \code{info.null} function calculates the information
  content per observation for any contrast between genotype means when
  densely genotyping an \code{alpha} fraction of
  the extreme phenotypic individuals.  The information content is
  calculated under the null hypothesis of no difference between the
  genotype means.  For small differences in genotype means, the
  information content will be approximately equal to the null, but in
  general, the information estimate under the null is the lower bound.

  The \code{info.bc.null} and \code{info.f2.null} functions calculate
  the information per observation for backcross, and F2 intercross
  respectively under the null hypothesis of equal gentoype means.  The
  information is calculated for a point in the middle of an interval
  spanned by markers separated by a recombination fraction \code{theta}.
}
\value{Information per individual for information functions, and the
  discount factor for the discount functions.
}
\references{Sen, Satagopan, and Churchill (2004), QTL study design from
  an information perspective,
  \url{http://repositories.cdlib.org/cbmb/QTLdesign}.
}
\author{Saunak Sen, Jaya Satagopan, and Gary Churchill}
\note{Only functions for backcross available right now.}
\seealso{\code{\link{info2cost.bc.null}}}
\examples{
info.null(0.5)
info.bc.null(0.5)
info.f2.null(0.5)
info.bc.null(0.5,0.1)
info.f2.null(0.5,0.1)
deflate.bc(0.1)
deflate.f2(0.1)
}
\keyword{design}
