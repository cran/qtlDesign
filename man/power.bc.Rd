\name{power.bc}
\alias{power.bc}
\alias{detectable.bc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Power calculations for Backcross}
\description{
Power calculation and minimum detectable effect sizes for backcross
}
\usage{
power.bc(n, prop, thresh = 3, alpha = 1, theta = 0, effective.n = FALSE)
detectable.bc(n, power=0.8, thresh = 3, alpha = 1, theta = 0, delta = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Sample size}
  \item{prop}{Proportion of variance explained by the QTL}
  \item{power}{Desired power level}
  \item{thresh}{LOD threshold for declaring significance}
  \item{alpha}{Selection fraction}
  \item{theta}{Width of marker interval}
  \item{effective.n}{Logical flag indicating whether effective sample size
    should be returned}
  \item{delta}{Logical flag indicating whether the QTL effect should be
    returned instead of proportion variance explained}
}
\details{
These calculations are done assuming that the asymptotic chi-square
regimes apply.  The function will not work if the effective sample size
is less than 30.  First we calculate the effective sample size using the
width of the marker interval and the selection fraction.  The QTL is
assumed to be in the middle of the marker interval.  Then we use the fact
that the non-centrality parameter of the likelihood ration test is
\eqn{ m*\delta^2}, where \eqn{m} is the effctive sample size and
\eqn{\delta} is the QTL effect measured as the deviation of the genotype
means from the overall mean.  The chi-squared approximation is used to
calculate the power.  The minimum detectable effect size (measured in
proportion variance explained or the QTL effect) is obtained by
numerical equation solving using \code{uniroot}.
}
\value{For \code{power.bc} the power is returned, unless the flag
  \code{effective.n} is set to \code{TRUE} in which case a list with
  values 
  \item{power}{Power of detecting QTL}
  \item{effective.n}{Effective sample size}
  is returned.

  For \code{detectable.bc} the minimum detectable proportion of variance
  explained is returned, unless the flag \code{delta} is set to
  \code{TRUE} in which case the QTL effect size is returned.
}
\references{Sen, Satagopan, and Churchill (2004), QTL study design from
  an information perspective,
  \url{http://repositories.cdlib.org/cbmb/QTLdesign}.
}
\author{Saunak Sen, Jaya Satagopan, and Gary Churchill}
\seealso{\code{\link{uniroot}}}
\examples{
power.bc(100,0.2,alpha=0.5,theta=0.1)
detectable.bc(100,0.8,alpha=0.5,theta=0.1)
}
\keyword{design}
