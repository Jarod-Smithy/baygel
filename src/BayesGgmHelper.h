#ifndef BayesGgmHelper_H
#define BayesGgmHelper_H

arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);

  return arma::repmat(mu, 1, n).t() + Y * chol(sigma);
}

#endif
