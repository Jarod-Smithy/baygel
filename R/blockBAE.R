#' Naïve Bayesian Adaptive Graphical Elastic Net
#'
#' A naïve Bayesian adaptive graphical elastic net data-augmented block Gibbs sampler.
#'
#' @param X Numeric matrix.
#' @param burnin An integer specifying the number of burn-in iterations.
#' @param nmc An integer specifying the number of MCMC samples.
#' @param tauPrior A numeric specifying the shrinkage hyperparameter for the off-diagonal taus.
#' @param lambdaPrior A numeric specifying the shrinkage hyperparameter for the off-diagonal lambdas.
#' @name blockBAE
#'
#' @return list containing:
#' \describe{
#' \item{Omega}{A \code{p} by \code{p} by nmc array of saved posterior samples of precision matrices.}
#' }
#'
#' @importFrom pracma triu Reshape
#' @importFrom stats rgamma rnorm
#' @importFrom statmod rinvgauss
#'
#'
#' @examples
#' \donttest{
#' # Generate true covariance matrix:
#' p             <- 10
#' n             <- 50
#' SigTrue       <- pracma::Toeplitz(c(0.7^rep(1:p-1)))
#' CTrue         <- pracma::inv(SigTrue)
#' # Generate expected value vector:
#' mu            <- rep(0,p)
#' # Generate multivariate normal distribution:
#' set.seed(123)
#' X             <- MASS::mvrnorm(n,mu=mu,Sigma=SigTrue)
#' omega_post    <- blockBAE(X,burnin = 1000,nmc = 500,tauPrior = 0.5,lambdaPrior = 0.05)
#'}
#' @export

blockBAE <- function(X,burnin = 1000,nmc = 2000,tauPrior,lambdaPrior){
  total_it <- burnin + nmc
  n <- nrow(X)
  p <- ncol(X)
  S <- t(X) %*% X
  Sig <- S/n
  C <- solve(Sig)
  indmx <- Reshape(1:p^2, p, p)
  upperind <- indmx[pracma::triu(indmx, 1) > 0]
  indmx_t <- t(indmx)
  lowerind <- indmx_t[pracma::triu(indmx_t, 1) > 0]
  C_save <- array(rep(0, p * p * nmc), dim = c(p, p, nmc))
  tau <- matrix(0, p, p)
  ind_noi_all <- matrix(0, p - 1, p)
  for (i in 1:p) {
    if (i == 1)
      ind_noi <- t(2:p)
    else if (i == p)
      ind_noi <- t(1:(p - 1))
    else ind_noi <- t(c(1:(i - 1), (i + 1):(p)))
    ind_noi_all[, i] <- ind_noi
  }
  lambda_ii <- 1

  for (iter in 1:(total_it)) {
    if (iter%%100 == 0) {
      cat("Total iterations= ", iter, "Iterations since burn in= ",
          ifelse(iter - burnin > 0, iter - burnin, 0),
          "\n")
    }
    Cadjust <- pmax(abs(as.vector(C)[upperind]), 10^-12)
    s_post <- 1
    t_post <- Cadjust + lambdaPrior
    lambda <- rgamma(length(t_post), shape = s_post, rate = t_post)
    a_post = 3/2
    b_post = (Cadjust^2)/2 + tauPrior
    sigma = 1 / rgamma(length(b_post), shape = a_post, rate = b_post)

    lambda_prime <- lambda^2
    mu_prime <- pmin(lambda/Cadjust, 1e+12)
    tau_temp <- 1/rinvgauss(length(mu_prime), mu_prime, lambda_prime)

    new <- (sigma+tau_temp)/(tau_temp*sigma)
    tau[upperind] <- new
    tau[lowerind] <- new

    for (i in 1:p) {
      ind_noi <- ind_noi_all[, i]
      tau_temp <- tau[ind_noi, i]
      Sig11 <- Sig[ind_noi, ind_noi]
      Sig12 = Sig[ind_noi, i]
      Sig12_scaled <- Sig12 / sqrt(Sig[i, i])
      invC11 <- Sig11 - tcrossprod(Sig12_scaled)
      Ci <- (S[i, i] + lambda_ii + 1) * invC11
      diag(Ci) <- diag(Ci) + tau_temp
      Ci <- (Ci + t(Ci))/2
      Ci_chol <- chol(Ci)
      mu_i <- -solve(Ci) %*% S[ind_noi, i]
      beta <- mu_i + solve(Ci_chol) %*% rnorm(p - 1)
      C[ind_noi, i] <- beta
      C[i, ind_noi] <- beta
      gam <- rgamma(1, (n/2) + 1, 1/(2/(S[i, i] + lambda_ii + 1)))
      C[i, i] <- gam + t(beta) %*% invC11 %*% beta
      Sig <- solve(C)

    }
    if (iter > burnin) {
      C_save[, , iter - burnin] <- C
    }
  }
  return(list(Omega = C_save))

}
