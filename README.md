
# **baygel** <a href='https://CRAN.R-project.org/package=baygel'><img src="man/figures/logo.png" align="right" height="150"/></a>

<!-- badges: start -->

![](https://www.r-pkg.org/badges/version/baygel)
![](https://www.r-pkg.org/badges/last-release/baygel)
![](https://cranlogs.r-pkg.org/badges/baygel)
![](https://cranlogs.r-pkg.org/badges/grand-total/baygel)
[![CodeFactor](https://www.codefactor.io/repository/github/jarod-smithy/baygel/badge)](https://www.codefactor.io/repository/github/jarod-smithy/baygel)
[![R-CMD-check](https://github.com/Jarod-Smithy/baygel/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Jarod-Smithy/baygel/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

The **baygel** `R` package provides data-augmented block Gibbs samplers,
for Bayesian shrinkage methods, to return the posterior distribution of
precision matrices for *Gaussian* distributed data with *positive
definite* covariance matrix. The package is implemented within the
following literature, including [Smith et
al. (2022)](https://doi.org/10.48550/arXiv.2210.16290) and [Smith et
al. (2023)](https://doi.org/10.48550/arXiv.2306.14199). The Bayesian
(adaptive) graphical lasso block Gibbs samplers of [H. Wang
(2012)](https://doi.org/10.1214/12-BA729) are also included for
convenience.

## Installation

You can install the latest version from CRAN using:

``` r
install.packages("baygel")
```

## Loading

``` r
library(baygel)
```

## Simple example

``` r
library(baygel)

# Generate true covariance matrix:
p             <- 10
n             <- 500
OmegaTrue     <- pracma::Toeplitz(c(0.7^rep(1:p-1)))
SigTrue       <- pracma::inv(OmegaTrue)
# Generate expected value vector:
mu            <- rep(0,p)
# Generate multivariate normal distribution:
set.seed(123)
X             <- MASS::mvrnorm(n, mu = mu, Sigma = SigTrue)
# Generate posterior distribution:
posterior     <- blockBAGR(X,iterations = 1000, burnin = 500)
# Estimated precision matrix
OmegaEst      <- apply(simplify2array(posterior$Omega), 1:2, mean)
```
