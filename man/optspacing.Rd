\name{Optimal marker spacing}
\alias{optspacing}
\alias{optspacing.bc}
\alias{optspacing.f2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Optimal marker spacing}
\description{
Functions to find optimal marker spacing given cost.
}
\usage{
optspacing(cost,G=NULL,alpha,cross)
optspacing.bc(cost,G=NULL,alpha)
optspacing.f2(cost,G=NULL,alpha)
optspacing(cost,G=NULL,alpha=NULL,cross)
optspacing.bc(cost,G=NULL,alpha=NULL)
optspacing.f2(cost,G=NULL,alpha=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cost}{Cost of genotyping in units of raising an individual}
  \item{alpha}{Selection fraction; proportion of individuals genotyped}
  \item{G}{Genome size in centiMorgans}
  \item{cross}{Cross type, "bc" or "f2"}  
}
\details{The function \code{optim} is used to search for the optima.}
}
\value{In the first formm, with the selection fraction specified, the
  spacing in centiMorgans that maximizes the information to cost ratio
  in the middle of the marker interval.  In the second form, with the
  selection fraction unspecified, it returns the value of
  (\code{spacing},\code{alpha}) which maximizes the information
  to cost ratio in the middle of the marker interval.}
\references{Sen S, Satagopan JM, Churchill GA (2005) 	
Quantitative trait locus study design from an information perspective.
Genetics, 170:447-64.}
\author{Saunak Sen, Jaya Satagopan, Karl Broman, and Gary Churchill}
\seealso{\code{\link{optim}}}
\examples{
optalpha(1,cross="bc")
optalpha(0.001,0.1,14.50,cross="bc")
optalpha(0.001,0.1,14.50,cross="f2")
}
\keyword{design}
