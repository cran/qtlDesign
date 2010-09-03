\name{mma}
\alias{mma}
\title{Selective phenotyping with similarity measure 2}
\description{
Selective phenotyping with similarity measure 2 to select the most
dissimilar subset of individuals.
}
\usage{
mma(genof, p, sequent = FALSE, exact = FALSE, dismat = FALSE)
}
\arguments{
  \item{genof}{Genotype matrix.}
  \item{p}{Sample size to select.}
  \item{sequent}{Perform sequential optimization if TRUE (see below).}
  \item{exact}{Count allele differences if \code{FALSE}; binary 0 = same
    number of alleles, 1 = different if \code{TRUE}.}
  \item{dismat}{Return dissimilarity matrix if TRUE.}
}
\details{
Sequentially minimize 1st moment and then 2nd moment, swapping one
subject at a time.
\code{op} finds all the samples with same 1st moment similarity with mma
results. \code{op2} finds all the samples with the same 1st moment
similarity with every list from op result. A combination of \code{op}
and \code{op2} comes very close to exhaustive search in
practice. \code{moment2} find the best list with minimum 2nd moments
from the output of \code{op2}. Note that some warnings occurs
accompanying our return statement. The results are not affected though.

This function combines several functions in Jin's original code.
\code{mma(genof,p,sequent=TRUE} is identical to the depricated
\code{mmasequent(genof,p}.
\code{mma(genof,p,exact=TRUE} is identical to the depricated
\code{mmaM1(genof,p} (actually, \code{mma} uses dissimilarity while
\code{mmaM1} used similarity = 1 - dissimilarity).

}
\value{
  A list containing \code{cList}, \code{dismat} if that option is
  \code{TRUE} and further optimized lists (\code{op}, \code{op2},
  \code{moment2}) if \code{sequent} is \code{TRUE}. 
vector as the first item. The list of items includes:
  \item{cList}{vector of selected subjects by function mma}
  \item{op}{list containing vector of selection and update flag from
    function op}
  \item{op2}{matrix of selection by function op2}
  \item{moment2}{vector of second moment calculations}
  \item{dismat}{dissimilarity matrix}
}
\references{Jin C, Lan H, Attie AD, Churchill GA, Bulutuglo D, Yandell BS (2004) Selective phenotyping for increased efficiency in genetic mapping studies. Genetics 168: 2285-2293.}
\author{Brian S. Yandell (\url{mailto:byandell@wisc.edu})}
\seealso{\code{\link{K1}}, \code{\link[qtl]{read.cross}}}
%\examples{}
\keyword{design}% at least one, from doc/KEYWORDS
