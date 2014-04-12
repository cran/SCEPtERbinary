\name{estimateBin}
\alias{estimateBin}
\title{Age estimates for the supplied binary systems}
\description{
The function computes the grid-based estimates of the age for
the supplied binary systems. 
 }
\usage{
    estimateBin(data, STAR, sigma, thr, sel, parallel=FALSE)
     }
\arguments{
\item{data}{the matrix estimation grid. The matrix should be sorted according to
  the values in the first column (the stellar effective
  temperature). }
\item{STAR}{a matrix of 18 columns, containing in each row the
  observational data of the stellar objects.} 
\item{sigma}{a vector of 14 elements, containing the error on observational
  data. The error on 
  mass and radius are expressed as relative error. If the errors on the
  observables of the two stars are identical, it is possible to provide a
  vector of length 7.}
\item{thr}{the threshold for the selection of cases over which the
  computation of likelihood is performed.}
\item{sel}{a vector of 7 elements, containing 1 or 0 for inclusion or
  exclusion of a observational data from the likelihood computation.}
\item{parallel}{logical flag indicating if the estimations should be
  performed in parallel on a multi-core architecture.}
}
\value{
  \code{STAR} should contains the observables for primary and
  secondary stars. In order, it should contains the stellar effective
  temperature, 
  the logarithm of the surface gravity, the metallicity [Fe/H], the
  average large frequency separation (in micro Hz), the maximum
  oscillation power (in micro Hz), the mass (in solar units), the radius
  (in solar units), the age (if known), and the relative age (if known).
  Identical quantities for the secondary stars should follows.
  \code{sigma} should contain, in the same order, the
  uncertainties on the observational quantities.  
  
  \code{estimateBin} returns a data frame obtained with the
  estimates
  of stellar parameters for each row of the input matrix \code{STAR}.
  In the columns of the data frame there are:
  the identifier \code{id}
  of the row of the object matrix adopted for the estimate; the
  independent star estimate of the primary mass
  \code{M1} (in solar units), radius \code{R1} (in
  solar units), age \code{age1} (in Gyr); the corresponding independent
  estimates for the secondary (\code{M2}, \code{R2}, \code{age2}); the
  joint-likelihood estimates of these quantities obtained explicitly assuming
  coevality (\code{M1b}, \code{R1b}, \code{age1b}, \code{M2b},
  \code{R2b}, \code{age2b}); the relative age \code{r} of the primary
  star (0 = ZAMS; 1 = TAMS); the best estimate of the age of the system
  \code{ageBin} under coevality assumption; the mean of the age of the two
  independent estimates of the stellar ages.  

  Age estimates are obtained by a maximum likelihood
  technique. Details on the technique can
  be found in the references reported below.
}
\references{
  
G. Valle, M. Dell'Omodarme, P.G. Prada Moroni, and S. Degl'Innocenti
(2014). \emph{Uncertainties in grid-based estimates of stellar mass and
  radius. SCEPtER: Stellar CharactEristics Pisa
  Estimation gRid}. A&A, 561, A125.

G. Valle, M. Dell'Omodarme, P.G. Prada Moroni, and S. Degl'Innocenti
(2014). \emph{Uncertainties in asteroseismic grid-based estimates of stellar ages. SCEPtER: Stellar CharactEristics Pisa
  Estimation gRid}. A&A, accepted.

G. Valle, M. Dell'Omodarme, P.G. Prada Moroni, and S. Degl'Innocenti
(2014). \emph{Grid-based estimates of stellar ages in binary
  systems. SCEPtER: Stellar CharactEristics Pisa Estimation gRid}. A&A, submitted.
}

\examples{
require(SCEPtER)
data(stdGrid)
ogrid <- stdGrid[ do.call(order, as.data.frame(stdGrid[,1])), ]

## observational constraint: Teff, [Fe/H], M, R
vsel <- c(1,0,1,0,0,1,1)
sigma <- c(100, 0.25, 0.1, 0.025, 0.05, 0.01, 0.005)

## syntetic sample of 10 bynary systems
bl <- block(stdGrid)
sam <- sampleBinStar(10, stdGrid, bl, restrict=TRUE)

## add Gaussian noise to the observations
starbp <- errorObsBin(sigma, sam)

resB <- estimateBin(ogrid, starbp, sigma, 3, vsel, parallel=FALSE)


## An observed system: PK Peg
## observation from Clausen et al. (2010). A&A, 516, A42.

obs <- c(6265,NA,-0.12,NA,NA,1.414,1.988,NA,NA, 6320,NA,-0.12,NA,NA,1.257,1.474,NA,NA)
sigma <- c(85,0,0.07,0,0,0.007/1.414,0.008/1.988, 90,0,0.07,0,0,0.005/1.257,0.017/1.474)

## for sensible estimate adopt at least nrow=10000
star <- matrix(obs, nrow=100, ncol=18, byrow=TRUE)

## add Gaussian noise
starp <- errorObsBin(sigma, star)

resB <- estimateBin(ogrid, starp, sigma, 3, vsel, parallel=FALSE)

quantile(resB$ageBin, c(0.16, 0.5, 0.84))
}
\keyword{manip}
