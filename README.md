
# **baygel** <a href='https://CRAN.R-project.org/package=baygel'><img src="man/figures/logo.png" align="right" height="150"/></a>

![](https://www.r-pkg.org/badges/version/baygel)
![](https://www.r-pkg.org/badges/last-release/baygel)
![](https://cranlogs.r-pkg.org/badges/baygel)
![](https://cranlogs.r-pkg.org/badges/grand-total/baygel)

## Overview

The **baygel** `R` package provides data-augmented block Gibbs samplers
to return the posterior distribution of precision matrices for
*Gaussian* distributed data with *positive definite* covariance matrix.
The package is implemented within the following literature, including
[Smith et al. (2022)](https://doi.org/10.48550/arXiv.2210.16290) and
[Smith et al. (2023)](https://doi.org/10.48550/arXiv.2306.14199).

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
n             <- 50
SigTrue       <- pracma::Toeplitz(c(0.7^rep(1:p-1)))
CTrue         <- pracma::inv(SigTrue)
# Generate expected value vector:
mu            <- rep(0,p)
# Generate multivariate normal distribution:
set.seed(123)
X             <- MASS::mvrnorm(n,mu=mu,Sigma=SigTrue)
posterior     <- blockBSGR(X,iterations = 1000, burnIn = 500)
```
