\name{Variance and effect size}
\alias{error.var}
\alias{gmeans2effect}
\alias{prop.var}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Effect size, proportion variance explained, and error variance
  calculations}
\description{
The function \code{error.var} estimates the error variance using
estimates of the environmental variance and genetic variance.  The effect
segregating at a locus, can be calculated using \code{gmeans2effect}
These are key inputs for power calculations.  The function
\code{prop.var} calculates the proportion of variance explained by a
locus given the effect size and error variance.
}
\usage{
error.var(cross,env.var=1,gen.var=0,bio.reps=1)
gmeans2effect(cross,means)
prop.var(cross,effect,sigma2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cross}{String indicating cross type which is "bc", for
    backcross, "f2" for intercross, and "ri" for recombinant inbred lines.}
  \item{env.var}{Environmental (within genotype) variance}
  \item{gen.var}{Genetic (between genotype) variance due to all loci
    segregating between the parental lines.}
  \item{bio.reps}{Number of biological replicates per unique genotype.
    This is usually 1 for backcross and intercross, but may be larger
    for RI lines.}
  \item{means}{Vector of genotype means in the form \code{c(a,h,b)},
    where \code{a} is the mean of the "AA" homozygotes, \code{h} is
    the mean of the "AB" heterozygotes, and \code{b} is the mean of the
    "BB" homozygotes.}
  \item{effect}{The QTL effect which depends on the cross.  For
    backcross, it is the difference in means the heterozygote and
    homozygote.  For RI lines it is half the difference in means of the
    homozygotes, for intercross, it is a two component vector of the form
    \code{c(a,d)}, where \code{a} is the additive effect (half the
    difference between the homozygotes), and \code{d} is the dominance
    effect (difference between the heterozygote and the average of the
    homozygotes). The genotype means will be \code{-a-d/2}, \code{d/2},
    and \code{a-d/2.} }
  \item{sigma2}{Error variance.}
}
\details{The function \code{error.var} estimates the error variance
  segregating in a cross using estimates of the environmental (within
  genotype) variance, and the genetic (between genotype variance).  The
  environmental variance is assumed to be invariant between cross types.
  The genetic variance segregating in RI lines is assumed to be double
  that in F2 intercross, and four times that of the backcross.  This
  assumption  holds if all loci are additive.  The error variance at a
  locus of interest is aproximately \deqn{\sigma_G^2/c + \sigma_E^2/m,}{%
    sigmaG^2/c + sigmaE^2/m,}
  where \eqn{$\sigma_G^2$}{sigmaG^2} is the genetic variance
  (\code{gen.var}), \eqn{c} is a
  constant depending on the cross type (1, for RI lines, 1/2 for F2
  intercross, and 1/4 for backross), \eqn{$\sigma_E^2$}{sigmaE^2} is the
  environmental
  variance (\code{env.var}), and \eqn{m} is the number of biological
  replicates per unique genotype (\code{bio.reps}).

  The function \code{gmeans2effect} calculates the QTL effects from
  genotype means depending on the cross.

  The function \code{prop.var} calculates the proportion of variance
  attributable to a locus given the effects size(s) and the error
  variance.  The definition of effect size is in the output of
  \code{gmeans2effect} (see below).
}

\value{For \code{error.var} the value is the estimated error variance
  based on the assumptions mentioned above.  For \code{gmeans2effect}
  the value depends on the type of cross.  For RI lines it is simply the
  additive effect of the QTL which is half the difference between the
  homozygote means.  For intercross, it is a vector giving the additive and
  dominance components.  The additive component is half the difference
  between the homozygote means, and the dominance component is the
  difference between the heterozygotes and the average of the
  homozygotes.  For the backcross, it is a vector of length 2,
  \code{c(a-h,h-b)}, which is the effect of an allelic substitution of
  an "A" allele in the backcrosses to each parental strain. }

\references{Sen S, Satagopan JM, Churchill GA (2005) 	
Quantitative trait locus study design from an information perspective.
Genetics, 170:447-64.}
\author{Saunak Sen, Jaya Satagopan, Karl Broman, and Gary Churchill}
\seealso{\code{\link{powercalc}}}
\examples{
error.var(cross="bc",env.var=1,gen.var=1,bio.reps=1)
error.var(cross="f2",env.var=1,gen.var=1,bio.reps=1)
error.var(cross="ri",env.var=1,gen.var=1,bio.reps=1)
error.var(cross="ri",env.var=1,gen.var=1,bio.reps=10)
gmeans2effect(cross="f2",means=c(0,1,2))
gmeans2effect(cross="f2",means=c(0,1,1))
gmeans2effect(cross="bc",means=c(0,1,1))
gmeans2effect(cross="ri",means=c(0,1,1))
prop.var(cross="bc",effect=5,sigma2=1)
prop.var(cross="f2",effect=c(5,0),sigma2=1)
prop.var(cross="ri",effect=5,sigma2=1)
}
\keyword{design}
