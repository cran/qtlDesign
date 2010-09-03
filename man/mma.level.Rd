\name{mma.level}
\alias{mma.level}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MMA utility}
\description{
This routine is for internal use. It sets 3 levels to 0,1,2.
}
\usage{
mma.level(mat)
}
\arguments{
  \item{mat}{input matrix}
}
\details{Converts matrix to levels between 0 and 2.
}
\value{
Matrix of genotype levels between 0 and 2.
}
\references{Jin C, Lan H, Attie AD, Churchill GA, Bulutuglo D, Yandell BS (2004) Selective phenotyping for increased efficiency in genetic mapping studies. Genetics 168: 2285-2293.}
\author{Brian S. Yandell (\url{mailto:byandell@wisc.edu})}
\seealso{\code{\link{mma}}, \code{\link[qtl]{read.cross}}}
%\examples{}
\keyword{design}% at least one, from doc/KEYWORDS

