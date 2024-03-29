% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{blockBAGENI}
\alias{blockBAGENI}
\title{Type I naive Bayesian adaptive graphical elastic net block Gibbs sampler for Gaussian graphical models.}
\usage{
blockBAGENI(
  X,
  burnin,
  iterations,
  verbose = TRUE,
  r = 0.001,
  s = 0.01,
  a = 0.001,
  b = 0.1
)
}
\arguments{
\item{X}{A numeric matrix, assumed to be generated from a multivariate Gaussian distribution.}

\item{burnin}{An integer specifying the number of burn-in iterations.}

\item{iterations}{An integer specifying the length of the Markov chain after the burn-in iterations.}

\item{verbose}{A logical determining whether the progress of the MCMC sampler should be displayed.}

\item{r}{A double specifying the value of the shape parameter for the gamma prior associated with the Bayesian graphical lasso penalty term.}

\item{s}{A double specifying the value of the scale parameter for the gamma prior associated with the Bayesian graphical lasso penalty term.}

\item{a}{A double specifying the value of the shape parameter for the inverse gamma prior associated with the Bayesian graphical ridge penalty term.}

\item{b}{A double specifying the value of the scale parameter for the inverse gamma prior associated with the Bayesian graphical ridge penalty term.}
}
\value{
A list containing precision `Omega` and covariance `Sigma` matrices
from the Markov chains.
}
\description{
Implements the Type I naive Bayesian adaptive graphical elastic net block Gibbs sampler to simulate the
posterior distribution of the precision matrix for Gaussian graphical models.
}
\examples{
# Generate true precision matrix:
p             <- 10
n             <- 500
OmegaTrue    <- pracma::Toeplitz(c(0.7^rep(1:p-1)))
SigTrue      <- pracma::inv(OmegaTrue)
# Generate expected value vector:
mu            <- rep(0,p)
# Generate multivariate normal distribution:
set.seed(123)
X             <- MASS::mvrnorm(n, mu = mu, Sigma = SigTrue)
# Generate posterior distribution:
posterior     <- blockBAGENI(X, iterations = 1000, burnin = 500)
# Estimated precision matrix using the mean of the posterior:
OmegaEst      <- apply(simplify2array(posterior$Omega), 1:2, mean)
}
