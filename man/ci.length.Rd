\name{Confidence interval expected widths}
\alias{ci.length}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculating expected QTL confidence interval widths}
\description{
  Provides expected confidence interval widths for QTL location when we
  have  dense markers.}
\usage{
ci.length(cross,n,effect,p=0.95,sigma2=1,env.var,gen.var,bio.reps=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cross}{String indicating cross type which is "bc", for
    backcross, "f2" for intercross, and "ri" for recombinant inbred lines.}
  \item{n}{Sample size}
  \item{p}{Confidence level for desired confidence interval}
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
  \item{sigma2}{Error variance; if this argument is absent,
    \code{env.var} and \code{gen.var} must be specified.}
  \item{env.var}{Environmental (within genotype) variance}
  \item{gen.var}{Genetic (between genotype) variance due to all loci
    segregating between the parental lines.}
  \item{bio.reps}{Number of biological replicates per unique genotype.
    This is usually 1 for backcross and intercross, but may be larger
    for RI lines.}
  }
\details{With dense markers, the log likelihood follows a compound
  process.  Approximate expected confidence intervals can be calculated
  by pretending the log likelihood decays linearly with a drift rate
  that depends on the effect size and cross type.
}

\value{Returns the expected confidence interval width in cM assuming
  dense markers.}

\references{
  Dupuis J and Siegmund D (1999) Statistical methods for
  mapping quantitative trait loci from a dense set of markers. Genetics
  151:373-386.
  
  Darvasi A (1998) Experimental strategies for the genetic
  dissection of complex traits in animal models.  Nature Genetics
  18:19-24.
  
  Kong A and Wright FA (1994) Asymptotic theory for gene mapping.
  Proceedings of the National Academy of Sciences of the USA
  91:9705-9709.
}
\author{Saunak Sen}
\seealso{\code{\link{powercalc}}.}
\examples{
ci.length(cross="bc",n=400,effect=5,p=0.95,sigma2=1)
}
\keyword{design}
