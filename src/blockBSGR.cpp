#include <RcppArmadillo.h>
#include "BayesGgmHelper.h"
#include "progress.hpp"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]

using namespace Rcpp;

//' Block Gibbs sampler function.
//'
//' A Bayesian standard graphical ridge-type data-augmented block Gibbs sampler for simulating the posterior distribution of the concentration matrix specifying a Gaussian graphical model.
//'
//' @param X Numeric data matrix, data is assumed to be Gaussian distributed.
//' @param burnIn An integer specifying the number of burn-in iterations.
//' @param iterations An integer specifying the length of the Markov chain after the burn-in iterations.
//' @param tau A double specifying the value of the variance parameter of both the Gaussian and truncated Gaussian distribution in Bayesian graphical ridge prior.
//' @param mu A double specifying the value of the mean parameter of both the Gaussian and truncated Gaussian distribution in Bayesian graphical ridge prior.
//' @param verbose A logical determining whether the progress of the MCMC sampler should be displayed.
//' @return blockBSGR: List of precision matrices from the Markov chains.
//' @examples
//'# Generate true covariance matrix:
//'p             <- 10
//'n             <- 50
//'SigTrue       <- pracma::Toeplitz(c(0.7^rep(1:p-1)))
//'CTrue         <- pracma::inv(SigTrue)
//'# Generate expected value vector:
//'mu            <- rep(0,p)
//'# Generate multivariate normal distribution:
//'set.seed(123)
//'X             <- MASS::mvrnorm(n,mu=mu,Sigma=SigTrue)
//'posterior     <- blockBSGR(X,iterations = 1000, burnIn = 500)
//' @export
// [[Rcpp::export]]
List blockBSGR(arma::mat X, int burnIn, int iterations,double tau = 1, double mu = 0, bool verbose = true) {

  // variable declarations and initialisations
  int totIter, n, p;
  totIter = burnIn + iterations;
  n = X.n_rows;
  p = X.n_cols;

  // setup progress bar
  Progress prog(totIter, verbose);

  arma::mat S(p, p), Sigma(p,p), Omega(p, p), perms(p-1, p);
  S = X.t()* X;
  Sigma = S/n;
  Omega = inv(Sigma);

  perms.fill(NA_REAL);
  arma::vec permInt(p);
  for(int i = 0; i < p; i++){
    permInt(i) = i;
  }

  for (int i = 0; i < p; i++) {
    if (i == 0){
      perms.col(i) = permInt.subvec(i+1,p-1);
    }
    if (i == p-1){
      perms.col(i) = permInt.subvec(0,p-2);
    }
    if ((i != p-1) & (i != 0)){
      arma::vec A, B;
      A = permInt.subvec(0,i-1);
      B = permInt.subvec(i+1,p-1);
      perms.col(i) = join_vert(A, B);
    }
  }

  List SigmaMatList(iterations);
  List OmegaMatList(iterations);
  arma::mat Sigma11(p-1,p-1);
  arma::mat Omega11(p-1,p-1);
  arma::mat S11(p-1,p-1);
  arma::vec Sigma12(p-1);
  arma::vec S_sub(p-1);
  arma::mat Omega11inv;
  arma::mat var(p-1,p-1,arma::fill::value(tau));
  arma::mat Ci;
  arma::vec mui;
  arma::rowvec beta(p-1);
  NumericVector gamm;
  int idx = 0;

  for (int iter=0; iter<totIter; iter++){

    if (Progress::check_abort())
      return -1.0;

    if (verbose){
      prog.increment(); // update progress
    }

    for (int i=0; i<p; i++){
      for (int j=0; j<(p-1); j++){
        Sigma12(j) = Sigma(perms.col(i)(j), i);
        S_sub(j) = S(perms.col(i)(j), i);
      }

      for (int j=0; j<(p-1); j++){
        for (int k=0; k<(p-1); k++){
          Sigma11(j,k) = Sigma(perms.col(i)(j), perms.col(i)(k));
          Omega11(j,k) = Omega(perms.col(i)(j), perms.col(i)(k));
          S11(j,k) = S(perms.col(i)(j), perms.col(i)(k));
        }
      }

      Omega11inv =Sigma11 - (Sigma12*Sigma12.t())/Sigma(i,i);
      Ci = (S(i,i) - 2*mu + 1)*Omega11inv + diagmat(var);
      mui = -inv(Ci) * S_sub;

      beta = mvrnormArma(1, mui, inv(Ci));
      for (int j=0; j<(p-1); j++){
        Omega(perms.col(i)(j),i)=Omega(i,perms.col(i)(j))=beta(j);
      }
      NumericVector gamm = Rcpp::rgamma(1, n/2 + 1, 1/((S(i,i) - 2*mu + 1)/2));
      Omega.submat(i,i,i,i) = gamm(0) + beta * Omega11inv * beta.t();
      Sigma = inv_sympd(Omega);

    }
    if (iter>=burnIn){
      OmegaMatList[idx] = Omega;
      idx += 1;
    }
  }
  return OmegaMatList;
}


