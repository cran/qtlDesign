\name{qtlDesign-internal}
\alias{fracmiss0}
\alias{fracmiss1}
\alias{fracmiss2}
\alias{missinfo0}
\alias{missinfo1}
\alias{missinfo2}
\alias{missinfo0.sim}
\alias{missinfo1.sim}
\alias{missinfo2.sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{qtlDesign-internal}
\description{
Internal qtlDesign functions
}
\usage{
fracmiss0(y,delta)
fracmiss1(y,m,delta,theta)
fracmiss2(y,m1,m2,delta,theta1,theta2)
missinfo0(delta,alpha)
missinfo1(delta,alpha,theta)
missinfo2(delta,alpha,theta1,theta2)
missinfo0.sim(delta,n,alpha)
missinfo1.sim(delta,n,alpha,theta)
missinfo2.sim(delta,n,alpha,theta1,theta2)
}
\arguments{
  \item{m}{Flanking marker genotype}
  \item{m1}{Flanking marker genotype}
  \item{m2}{Flanking marker genotype}  
  \item{theta}{Recombination fraction}
  \item{theta1}{Recombination fraction}
  \item{theta2}{Recombination fraction}  
  \item{y}{Phenotype value}
  \item{n}{Sample size}
  \item{delta}{QTL effect}
  \item{alpha}{Selection fraction}
}
\details{
  These are not to be called by the user.
}

\references{Sen, Satagopan, and Churchill (2004), QTL study design from
  an information perspective,
  \url{http://repositories.cdlib.org/cbmb/QTLdesign}.
}
\author{Saunak Sen, Jaya Satagopan, and Gary Churchill}

\seealso{\code{\link{missinfo}},\code{\link{missinfo.sim}},
   \code{\link{fracmiss}} }
\examples{
}
\keyword{internal}
