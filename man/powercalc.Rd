\name{Power calculations}
\alias{powercalc}
\alias{detectable}
\alias{samplesize}
\alias{powercalc.bc}
\alias{detectable.bc}
\alias{samplesize.bc}
\alias{powercalc.f2}
\alias{detectable.f2}
\alias{samplesize.f2}
\alias{powercalc.ri}
\alias{detectable.ri}
\alias{samplesize.ri}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Power, sample size, and detectable effect size calculations }
\description{
Power, sample size, and minimum detectable effect size calculations are
performed for backcross, F2 intercross, and recombinant inbred (RI)
lines.}
\usage{
powercalc(cross,n,effect,sigma2,env.var,gen.var,thresh=3,sel.frac=1,
          theta=0,bio.reps=1)
detectable(cross,n,effect=NULL,sigma2,env.var,gen.var,power=0.8,thresh=3,
           sel.frac=1,theta=0,bio.reps=1)
samplesize(cross,effect,sigma2,env.var,gen.var,power=0.8,thresh=3,
           sel.frac=1,theta=0,bio.reps=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cross}{String indicating cross type which is "bc", for
    backcross, "f2" for intercross, and "ri" for recombinant inbred lines.}
  \item{n}{Sample size}
  \item{sigma2}{Error variance; if this argument is absent,
    \code{env.var} and \code{gen.var} must be specified.}
  \item{env.var}{Environmental (within genotype) variance}
  \item{gen.var}{Genetic (between genotype) variance due to all loci
    segregating between the parental lines.}
  \item{effect}{The QTL effect we want to detect.  For
    \code{powercalc} and \code{samplesize} this is a numeric (vector).
    For \code{detectable} it specifies the relative magnitude of the
    additive and dominance components for the intercross.
    The specification of \code{effect} depends on the cross.  For
    backcross, it is the difference in means the heterozygote and
    homozygote.  For RI lines it is half the difference in means of the
    homozygotes, for intercross, it is a two component vector of the form
    \code{c(a,d)}, where \code{a} is the additive effect (half the
    difference between the homozygotes), and \code{d} is the dominance
    effect (difference between the heterozygote and the average of the
    homozygotes). The genotype means will be \code{-a-d/2}, \code{d/2},
    and \code{a-d/2.}  For \code{detectable}, optionally for the
    intercross, one can use a string to specify the QTL effect type.
    The  strings "add" or "dom" are used to denote an additive or
    dominant model respectively for the phenotype.  It may be 
    it can be a numerical vector of the form \code{c(a,d)} indicating
    the relative magnitudes of the additive and dominance components (as
    defined above).  The default is "add".}
  \item{power}{Proportion indicating power desired}
  \item{thresh}{LOD threshold for declaring significance}
  \item{sel.frac}{Selection fraction}
  \item{theta}{Recombination fraction corresponding to a marker interval}
  \item{bio.reps}{Number of biological replicates per unique genotype.
    This is usually 1 for backcross and intercross, but may be larger
    for RI lines.}
}
\details{
These calculations are done assuming that the asymptotic chi-square
regimes apply.  A warning message is printed if the effective sample size
is less than 30 and either \code{sel.frac} is less than 1 or \code{theta}
is greater than 0.  First we calculate the effective sample size using the
width of the marker interval and the selection fraction.  The QTL is
assumed to be in the middle of the marker interval.  Then we use the fact
that the non-centrality parameter of the likelihood ration test is
\eqn{m*\delta^2}, where \eqn{m} is the effctive sample size and
\eqn{\delta} is the QTL effect measured as the deviation of the genotype
means from the overall mean.  The chi-squared approximation is used to
calculate the power.  The minimum detectable effect size is obtained by
solving the power equation numerically using \code{uniroot}.  The theory
behind the information calculations is described by Sen et. al. (2005).

A key input is the error variance, \code{sigma2} which is generally
unknown.  The user can enter the error variance directly, or estimate it
using \code{env.var} and \code{gen.var}.  The function \code{error.var}
is used to the error variance using estimates of the environmental variance
and genetic variance.  Another key input is the effect segregating in
a cross, which can be calculated using \code{gmeans2model}.  }

\value{For \code{powercalc} the power is returned, along with the
proportion of variance explained.  For \code{detectable} the effect size
detectable is returned, along with the proportion of variance explained.
For backcross and RI lines this is the effect of an allelic
substitution.  For F2 intercross the additive and dominance components
are returned.  For \code{samplesize} the sample size (rounded up to the
nearest integer) is returned along with the proportion of variance
explained.}

\references{Sen S, Satagopan JM, Churchill GA (2005) 	
Quantitative trait locus study design from an information perspective.
Genetics, 170:447-64.}
\author{Saunak Sen, Jaya Satagopan, Karl Broman, and Gary Churchill}
\seealso{\code{\link{uniroot}}. \code{\link{error.var}},
  \code{\link{gmeans2effect}}.}
\examples{
powercalc("bc",100,5,sigma2=1,sel.frac=1,theta=0)
powercalc(cross="ri",n=30,effect=5,env.var=64,gen.var=25,bio.rep=6)
detectable("bc",100,sigma2=1)
detectable(cross="ri",n=30,env.var=64,gen.var=25,bio.rep=8)
samplesize(cross="f2",effect=c(5,0),env.var=64,gen.var=25)
}
\keyword{design}
