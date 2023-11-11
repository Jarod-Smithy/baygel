#ifndef BayesGgmHelper_H
#define BayesGgmHelper_H

#include <RcppArmadillo.h>
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma);
arma::vec rrinvgauss(int n, double mu, double lambda);
arma::vec extract_upperind(arma::mat A);
Rcpp::NumericVector rrgamma(int n, double shape, double scale);
#endif
