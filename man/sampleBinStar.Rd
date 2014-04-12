\name{sampleBinStar}
\alias{sampleBinStar}
\title{Sample a set of binary systems out of a grid}
\description{
The function extracts, without replacement, a sample of binary systems out of a grid.  
 }
\usage{
    sampleBinStar(n, grid, block, restrict=TRUE, parallel=FALSE)
     }
\arguments{
\item{n}{the number of objects to be sampled.}
\item{grid}{the grid from which the objects are sampled.}
\item{block}{the length of blocks of the same initial metallicity [Fe/H]}
\item{restrict}{logical flag indicating if only models younger than 14 Gyr
  should be considered.}
\item{parallel}{logical flag indicating if the estimations should be
  performed in parallel on a multi-core architecture.}
}
\value{

  \code{sampleBinStar} returns a matrix of \code{n} rows, sampled without
  replacement from \code{grid}.
  A maximum age difference of 10 Myr is allowed in the binary system
  selection. 

}

\examples{
require(SCEPtER)
data(stdGrid)

bl <- block(stdGrid)
sam <- sampleBinStar(10, stdGrid, bl, restrict=TRUE)
}
\keyword{manip}
