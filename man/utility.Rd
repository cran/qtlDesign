\name{Utility}
\alias{recomb}
\alias{genetic.dist}
\alias{pmixnorm}
\alias{prop2delta.bc}
\alias{delta2prop.bc}  
\alias{prop2delta.f2}
\alias{delta2prop.f2}  
\alias{prop2gv}
\alias{gv2prop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Utility functions}
\description{Utility functions
}
\usage{
recomb(d)
genetic.dist(theta)
pmixnorm(x,mean=c(0,0),sd=c(1,1),mix.prop=0.5,level=0)
prop2delta.bc(prop)
delta2prop.bc(delta)
prop2delta.f2(prop,model)
delta2prop.f2(delta,model)
prop2gv(prop)
gv2prop(gv)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{Genetic distance in Morgans}
  \item{theta}{Recombination fraction}
  \item{x}{Real number}
  \item{mean}{Means of component mixtures}
  \item{sd}{Standard deviations of component mixtures}
  \item{mix.prop}{Mixing proportion}
  \item{level}{Subtracting constant}
  \item{prop}{Proportion of variance explained by a QTL}
  \item{delta}{QTL effect size}
  \item{model}{Type of genetic model}
  \item{gv}{Genetic variance}  
}
\details{
}
\value{\code{recomb} returns the recombination fraction
  corresponding to a genetic distance in Morgans.  \code{genetic.dist}
  returns the genetic distance in Morgans for a recombination fraction.
  \code{pmixnorm} returns the probability
  mixture of two Normal distributions with means given by \code{mean},
  standard deviations \code{sd}, and mixing proportions \code{mix.prob},
  is less than \code{x}. \code{prop2delta.bc} returns the QTL effect from
  the proportion of variance explained for a backcross.
  \code{delta2prop.bc} returns the QTL effect from the proportion of
  variance explained in a backcross.
  \code{prop2delta.f2} returns the QTL effect from
  the proportion of variance explained and genetic model for a F2 intercross.
  \code{delta2prop.f2} returns the QTL effect from the proportion of
  variance explained and the genetic model in an F2 intercross.
  \code{prop2gv} returns the genetic variance (as a fraction of the
  environmental variance) given the proportion of variance explained by
  a QTL. \code{gv2prop} performs the inverse function.
}

\references{Sen, Satagopan, and Churchill (2004), QTL study design from
  an information perspective,
  \url{http://repositories.cdlib.org/cbmb/QTLdesign}.
}
\author{Saunak Sen, Jaya Satagopan, and Gary Churchill}
\note{We assume Haldane mapping function for the genetic distance.}

\seealso{\code{\link{missinfo}}}
\examples{
recomb(0.1)
genetic.dist(0.1)
prop2delta.bc(0.1)
delta2prop.bc(0.1)
}
\keyword{utilities}
