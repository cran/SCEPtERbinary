\name{errorObsBin}
\alias{errorObsBin}
\title{Perturbe a set of observational data}
\description{
The function returns a perturbed object, starting from observational data
and a vector of uncertainties. 
 }
\usage{
    errorObsBin(sigma, STAR, parallel=FALSE, corr=c(0,0,0,0,0,0,0))
     }
\arguments{
\item{sigma}{a vector of 14 elements containing the uncertainties on
  observational data. The error on 
  mass and radius must be expressed as relative error. If the errors on
  the observations of the two stars are identical, it is possible to
  provide a vector of length 7.
} 
\item{STAR}{a vector of 18 elements containing  the observational data.}
\item{parallel}{logical flag specifying if the computation should be
  performed in parallel on a multi-core architecture.}
\item{corr}{a vector with the assumed the correlation between corresponding
  observables of the two stars. Default to zero (no correlations).}
}
\value{
  \code{errorObsBin} returns an object obtained sampling from a
  multivariate normal distribution with vector of mean \code{STAR} and
  covariance matrix computed
  according to the value of \code{corr}. 
}
\details{
  The vector of observations \code{STAR} must contain, in the order:
  \enumerate{
  \item{The effective temperature (in K) of the star.}
  \item{The logarithm of the surface gravity of the star (in cm s^-2).}
  \item{The metallicity [Fe/H] of the star.}
  \item{The value of large frequency separation of the star, divided by
    the solar value of this quantity. }
  \item{The value of the frequency of maximum oscillation power of the
  star, divided by 
  the solar value of this quantity. }
\item{The mass (in solar units) of the star. }
\item{The radius (in solar unit) of the star. }
}

  The vector of the uncertainties on the observation \code{sigma} must
  contain the uncertainty on the above quantity. The
  uncertainty on the asteroseismic parameters, mass, and radius 
  must be expressed as relative error.
  }
\keyword{manip}
