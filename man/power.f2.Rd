\name{power.f2}
\alias{power.f2}
\alias{detectable.f2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Power calculations for F2 intercross}
\description{
Power calculation and minimum detectable effect sizes for F2 intercross
}
\usage{
power.f2(n, model, thresh = 3, alpha = 1, theta = 0, effective.n = FALSE)
detectable.f2(n, model="add", power=0.8, thresh = 3, alpha = 1, theta =
0, delta=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Sample size}
  \item{model}{Genetic model written as a string which should be either
    "add" or "dom" or a numeric vector \code{c(a,d)}, which are the
    additive (half the difference between
     the homozygotes) and dominance (the difference
     between the heterozygote and the midpoint between the
     homozygotes) components
    respectively. The value "add" corresponds to \code{d=0}, and "dom"
    corresponds to \code{d=a}.  The genotype means will be
    \code{-a-d/2}, \code{d/2}, and \code{a-d/2}. }
  \item{power}{Desired power level}
  \item{thresh}{LOD threshold for declaring significance}
  \item{alpha}{Selection fraction}
  \item{theta}{Width of marker interval}
  \item{effective.n}{Logical flag indicating whether effective sample size
    should be returned}
  \item{delta}{Logical flag indicating whether the genetic model should
    be returned instead of the proportion of variance explained.}
}

\details{ These calculations are done assuming that the asymptotic
chi-square regimes apply.  The function will not work if the effective
sample size is less than 30.  First we calculate the effective sample
size using the width of the marker interval and the selection fraction.
The QTL is assumed to be in the middle of the marker interval.  Then we
use the fact that the non-centrality parameter of the likelihood ration
test is can be calculated given the information matrix and the QTL
effects.  The chi-squared approximation is used to calculate the power.
The minimum detectable effect size is obtained by numerical equation
solving using \code{uniroot}.}

\value{For \code{power.f2} the power is returned, unless the flag
  \code{effective.n} is set to \code{TRUE} in which case a list with
  values 
  \item{power}{Power of detecting QTL}
  \item{effective.n}{Effective sample size}
  is returned.

  For \code{detectable.f2} the minimum detectable proportion of variance
  explained is returned.
}
\references{Sen, Satagopan, and Churchill (2004), QTL study design from
  an information perspective,
  \url{http://repositories.cdlib.org/cbmb/QTLdesign}.
}
\author{Saunak Sen, Jaya Satagopan, and Gary Churchill}
\seealso{\code{\link{uniroot}}}
\examples{
power.f2(100,c(3/4,0),alpha=0.5,theta=0.1)
detectable.f2(100,alpha=0.5,theta=0.1)
}
\keyword{design}
