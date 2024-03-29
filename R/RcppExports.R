# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Type I naive Bayesian adaptive graphical elastic net block Gibbs sampler for Gaussian graphical models.
#'
#' Implements the Type I naive Bayesian adaptive graphical elastic net block Gibbs sampler to simulate the
#' posterior distribution of the precision matrix for Gaussian graphical models.
#'
#' @param X A numeric matrix, assumed to be generated from a multivariate Gaussian distribution.
#' @param burnin An integer specifying the number of burn-in iterations.
#' @param iterations An integer specifying the length of the Markov chain after the burn-in iterations.
#' @param a A double specifying the value of the shape parameter for the inverse gamma prior associated with the Bayesian graphical ridge penalty term.
#' @param b A double specifying the value of the scale parameter for the inverse gamma prior associated with the Bayesian graphical ridge penalty term.
#' @param r A double specifying the value of the shape parameter for the gamma prior associated with the Bayesian graphical lasso penalty term.
#' @param s A double specifying the value of the scale parameter for the gamma prior associated with the Bayesian graphical lasso penalty term.
#' @param verbose A logical determining whether the progress of the MCMC sampler should be displayed.
#' @return A list containing precision `Omega` and covariance `Sigma` matrices
#' from the Markov chains.
#' @examples
#'# Generate true precision matrix:
#'p             <- 10
#'n             <- 500
#' OmegaTrue    <- pracma::Toeplitz(c(0.7^rep(1:p-1)))
#' SigTrue      <- pracma::inv(OmegaTrue)
#'# Generate expected value vector:
#'mu            <- rep(0,p)
#'# Generate multivariate normal distribution:
#'set.seed(123)
#'X             <- MASS::mvrnorm(n, mu = mu, Sigma = SigTrue)
#'# Generate posterior distribution:
#'posterior     <- blockBAGENI(X, iterations = 1000, burnin = 500)
#'# Estimated precision matrix using the mean of the posterior:
#'OmegaEst      <- apply(simplify2array(posterior$Omega), 1:2, mean)
#' @export
blockBAGENI <- function(X, burnin, iterations, verbose = TRUE, r = 1e-3, s = 1e-2, a = 1e-3, b = 1e-1) {
    .Call(`_baygel_blockBAGENI`, X, burnin, iterations, verbose, r, s, a, b)
}

#' Type II naive Bayesian adaptive graphical elastic net block Gibbs sampler for Gaussian graphical models.
#'
#' Implements the Type II naive Bayesian adaptive graphical elastic net block Gibbs sampler to simulate the
#' posterior distribution of the precision matrix for Gaussian graphical models.
#'
#' @param X A numeric matrix, assumed to be generated from a multivariate Gaussian distribution.
#' @param burnin An integer specifying the number of burn-in iterations.
#' @param iterations An integer specifying the length of the Markov chain after the burn-in iterations.
#' @param b A double specifying the value of the rate parameter for the exponential prior associated with the Bayesian graphical ridge penalty term.
#' @param s A double specifying the value of the rate parameter for the exponential prior associated with the Bayesian graphical lasso penalty term.
#' @param verbose A logical determining whether the progress of the MCMC sampler should be displayed.
#' @return A list containing precision `Omega` and covariance `Sigma` matrices
#' from the Markov chains.
#' @examples
#'# Generate true precision matrix:
#'p             <- 10
#'n             <- 500
#' OmegaTrue    <- pracma::Toeplitz(c(0.7^rep(1:p-1)))
#' SigTrue      <- pracma::inv(OmegaTrue)
#'# Generate expected value vector:
#'mu            <- rep(0,p)
#'# Generate multivariate normal distribution:
#'set.seed(123)
#'X             <- MASS::mvrnorm(n, mu = mu, Sigma = SigTrue)
#'# Generate posterior distribution:
#'posterior     <- blockBAGENII(X, iterations = 1000, burnin = 500)
#'# Estimated precision matrix using the mean of the posterior:
#'OmegaEst      <- apply(simplify2array(posterior$Omega), 1:2, mean)
#' @export
blockBAGENII <- function(X, burnin, iterations, verbose = TRUE, s = 1e-1, b = 1e-3) {
    .Call(`_baygel_blockBAGENII`, X, burnin, iterations, verbose, s, b)
}

#' Bayesian adaptive graphical lasso block Gibbs sampler for Gaussian graphical models.
#'
#' Implements a Bayesian adaptive graphical lasso block Gibbs sampler to simulate the
#' posterior distribution of the precision matrix for Gaussian graphical models.
#'
#' @param X A numeric matrix, assumed to be generated from a multivariate Gaussian distribution.
#' @param burnin An integer specifying the number of burn-in iterations.
#' @param iterations An integer specifying the length of the Markov chain after the burn-in iterations.
#' @param r A double specifying the value of the shape parameter for the gamma prior.
#' @param s A double specifying the value of the scale parameter for the gamma prior.
#' @param verbose A logical determining whether the progress of the MCMC sampler should be displayed.
#' @return A list containing precision `Omega` and covariance `Sigma` matrices
#' from the Markov chains.
#' @examples
#'# Generate true precision matrix:
#'p             <- 10
#'n             <- 500
#' OmegaTrue    <- pracma::Toeplitz(c(0.7^rep(1:p-1)))
#' SigTrue      <- pracma::inv(OmegaTrue)
#'# Generate expected value vector:
#'mu            <- rep(0,p)
#'# Generate multivariate normal distribution:
#'set.seed(123)
#'X             <- MASS::mvrnorm(n, mu = mu, Sigma = SigTrue)
#'# Generate posterior distribution:
#'posterior     <- blockBAGL(X, iterations = 1000, burnin = 500)
#'# Estimated precision matrix using the mean of the posterior:
#'OmegaEst      <- apply(simplify2array(posterior$Omega), 1:2, mean)
#' @export
blockBAGL <- function(X, burnin, iterations, verbose = TRUE, r = 1e-2, s = 1e-6) {
    .Call(`_baygel_blockBAGL`, X, burnin, iterations, verbose, r, s)
}

#' Bayesian adaptive graphical ridge block Gibbs sampler for Gaussian graphical models.
#'
#' Implements a Bayesian adaptive graphical ridge block Gibbs sampler to simulate the
#' posterior distribution of the precision matrix for Gaussian graphical models.
#'
#' @param X A numeric matrix, assumed to be generated from a multivariate Gaussian distribution.
#' @param burnin An integer specifying the number of burn-in iterations.
#' @param iterations An integer specifying the length of the Markov chain after the burn-in iterations.
#' @param a A double specifying the value of the shape parameter for the inverse gamma prior.
#' @param b A double specifying the value of the scale parameter for the inverse gamma prior.
#' @param verbose A logical determining whether the progress of the MCMC sampler should be displayed.
#' @return A list containing precision `Omega` and covariance `Sigma` matrices
#' from the Markov chains.
#' @examples
#'# Generate true precision matrix:
#'p             <- 10
#'n             <- 500
#' OmegaTrue    <- pracma::Toeplitz(c(0.7^rep(1:p-1)))
#' SigTrue      <- pracma::inv(OmegaTrue)
#'# Generate expected value vector:
#'mu            <- rep(0,p)
#'# Generate multivariate normal distribution:
#'set.seed(123)
#'X             <- MASS::mvrnorm(n, mu = mu, Sigma = SigTrue)
#'# Generate posterior distribution:
#'posterior     <- blockBAGR(X, iterations = 1000, burnin = 500)
#'# Estimated precision matrix using the mean of the posterior:
#'OmegaEst      <- apply(simplify2array(posterior$Omega), 1:2, mean)
#' @export
blockBAGR <- function(X, burnin, iterations, verbose = TRUE, a = 1, b = 1e-2) {
    .Call(`_baygel_blockBAGR`, X, burnin, iterations, verbose, a, b)
}

#' Naive Bayesian graphical elastic net block Gibbs sampler for Gaussian graphical models.
#'
#' Implements the Bayesian graphical elastic net block Gibbs sampler to simulate the
#' posterior distribution of the precision matrix for Gaussian graphical models.
#'
#' @param X A numeric matrix, assumed to be generated from a multivariate Gaussian distribution.
#' @param burnin An integer specifying the number of burn-in iterations.
#' @param iterations An integer specifying the length of the Markov chain after the burn-in iterations.
#' @param lambda A numeric value representing the rate parameter for the double
#' exponential and exponential prior associated with the Bayesian graphical lasso penalty term.
#' @param sig A numeric value representing the standard deviation parameter for the double
#' Gaussian and truncated Gaussian prior associated with the Bayesian graphical ridge penalty term.
#' @param verbose A logical determining whether the progress of the MCMC sampler should be displayed.
#' @return A list containing precision `Omega` and covariance `Sigma` matrices
#' from the Markov chains.
#' @examples
#'# Generate true precision matrix:
#'p             <- 10
#'n             <- 500
#' OmegaTrue    <- pracma::Toeplitz(c(0.7^rep(1:p-1)))
#' SigTrue      <- pracma::inv(OmegaTrue)
#'# Generate expected value vector:
#'mu            <- rep(0,p)
#'# Generate multivariate normal distribution:
#'set.seed(123)
#'X             <- MASS::mvrnorm(n, mu = mu, Sigma = SigTrue)
#'# Generate posterior distribution:
#'posterior     <- blockBGEN(X, iterations = 1000, burnin = 500, lambda = 1, sig = 1)
#'# Estimated precision matrix using the mean of the posterior:
#'OmegaEst      <- apply(simplify2array(posterior$Omega), 1:2, mean)
#' @export
blockBGEN <- function(X, burnin, iterations, lambda = 1, sig = 1, verbose = TRUE) {
    .Call(`_baygel_blockBGEN`, X, burnin, iterations, lambda, sig, verbose)
}

#' Bayesian graphical lasso block Gibbs sampler for Gaussian graphical models.
#'
#' Implements a Bayesian graphical lasso block Gibbs sampler to simulate the
#' posterior distribution of the precision matrix for Gaussian graphical models.
#'
#' @param X A numeric matrix, assumed to be generated from a multivariate Gaussian distribution.
#' @param burnin An integer representing the number of burn-in iterations.
#' @param iterations An integer representing the length of the Markov chain post burn-in.
#' @param lambda A numeric value representing the rate parameter for the double
#' exponential and exponential prior.
#' @param verbose A logical indicating if the MCMC sampler progress should be printed.
#' @return A list containing precision `Omega` and covariance `Sigma` matrices
#' from the Markov chains.
#' @examples
#'# Generate true precision matrix:
#'p             <- 10
#'n             <- 500
#' OmegaTrue    <- pracma::Toeplitz(c(0.7^rep(1:p-1)))
#' SigTrue      <- pracma::inv(OmegaTrue)
#'# Generate expected value vector:
#'mu            <- rep(0,p)
#'# Generate multivariate normal distribution:
#'set.seed(123)
#'X             <- MASS::mvrnorm(n, mu = mu, Sigma = SigTrue)
#'# Generate posterior distribution:
#'posterior     <- blockBGL(X, iterations = 1000, burnin = 500, lambda = 0.5)
#'# Estimated precision matrix using the mean of the posterior:
#'OmegaEst      <- apply(simplify2array(posterior$Omega), 1:2, mean)
#' @export
blockBGL <- function(X, burnin, iterations, lambda = 1, verbose = TRUE) {
    .Call(`_baygel_blockBGL`, X, burnin, iterations, lambda, verbose)
}

#' Bayesian graphical ridge block Gibbs sampler for Gaussian graphical models.
#'
#' Implements a Bayesian graphical ridge block Gibbs sampler to simulate the
#' posterior distribution of the precision matrix for Gaussian graphical models.
#' 
#' @name blockBGR
#' @param X A numeric matrix, assumed to be generated from a multivariate Gaussian distribution.
#' @param burnin An integer representing the number of burn-in iterations.
#' @param iterations An integer representing the length of the Markov chain post burn-in.
#' @param sig A numeric value representing the standard deviation parameter for the double
#' Gaussian and truncated Gaussian prior.
#' @param verbose A logical indicating if the MCMC sampler progress should be printed.
#' @return A list containing precision `Omega` and covariance `Sigma` matrices
#' from the Markov chains.
#' @examples
#'# Generate true precision matrix:
#'p             <- 10
#'n             <- 500
#' OmegaTrue    <- pracma::Toeplitz(c(0.7^rep(1:p-1)))
#' SigTrue      <- pracma::inv(OmegaTrue)
#'# Generate expected value vector:
#'mu            <- rep(0,p)
#'# Generate multivariate normal distribution:
#'set.seed(123)
#'X             <- MASS::mvrnorm(n, mu = mu, Sigma = SigTrue)
#'# Generate posterior distribution:
#'posterior     <- blockBGR(X, iterations = 1000, burnin = 500, sig = 0.5)
#'# Estimated precision matrix using the mean of the posterior:
#'OmegaEst      <- apply(simplify2array(posterior$Omega), 1:2, mean)
#' @export
blockBGR <- function(X, burnin, iterations, sig = 1, verbose = TRUE) {
    .Call(`_baygel_blockBGR`, X, burnin, iterations, sig, verbose)
}
