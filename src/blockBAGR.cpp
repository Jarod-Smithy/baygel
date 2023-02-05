#include <RcppArmadillo.h>
#include "BayesGgmHelper.h"
#include "progress.hpp"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]

using namespace Rcpp;

//' Block Gibbs sampler function.
//'
//' A Bayesian adaptive graphical ridge-type data-augmented block Gibbs sampler for simulating the posterior distribution of the concentration matrix specifying a Gaussian graphical model.
//'
//' @param X Numeric data matrix, data is assumed to be Gaussian distributed.
//' @param burnIn An integer specifying the number of burn-in iterations.
//' @param iterations An integer specifying the length of the Markov chain after the burn-in iterations.
//' @param s A double specifying the value of the prior inverse gamma's shape parameter.
//' @param t A double specifying the value of the prior inverse gamma's scale parameter.
//' @param verbose A logical determining whether the progress of the MCMC sampler should be displayed.
//' @return blockBAGR: List of precision matrices from the Markov chains.
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
//'posterior     <- blockBAGR(X,iterations = 1000, burnIn = 500)
//' @export
// [[Rcpp::export]]
List blockBAGR(arma::mat X, int burnIn, int iterations,double s = 1,double t = 1, bool verbose = true) {

  // variable declarations and initialisations
  int totIter, n, p;
  totIter = burnIn + iterations;
  n = X.n_rows;
  p = X.n_cols;

  // setup progress bar
  Progress prog(totIter, verbose);

  arma::mat S(p, p), Sigma(p,p), Omega(p, p), perms(p-1, p), tau(p, p, arma::fill::zeros);
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
      perms.col(i) = permInt.subvec(i+1,p-1);//;
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
  arma::mat OmegaTemp;
  arma::uvec lw_idx;
  arma::rowvec tauI(p-1);
  arma::mat Sigma11(p-1,p-1);
  arma::mat Omega11(p-1,p-1);
  arma::mat S11(p-1,p-1);
  arma::vec Sigma12(p-1);
  arma::vec S_sub(p-1);
  arma::mat Omega11inv;
  arma::mat Ci;
  arma::vec mui;
  arma::rowvec beta(p-1);
  NumericVector gamm;
  double s_post;
  double t_post;
  int idx = 0;

  for (int iter=0; iter<totIter; iter++){

    if (Progress::check_abort())
      return -1.0;

    if (verbose){
      prog.increment(); // update progress
    }

    arma::uvec lw_idx = arma::trimatl_ind(arma::size(Omega), -1);
    OmegaTemp = Omega(lw_idx);
    arma::vec low_tri_tau_vals(OmegaTemp.size());
    for (int i=0; i<low_tri_tau_vals.size(); i++){
      s_post = 0.5 + s;
      t_post = std::max<double>(1e-6,OmegaTemp(i)*OmegaTemp(i))/2.0 + t;
      low_tri_tau_vals(i)= 1.0/Rcpp::rgamma(1, s_post, 1/t_post)(0);
    }
    arma::uvec tau_lw_idx = arma::trimatl_ind(arma::size(tau), -1);
    tau.elem(tau_lw_idx) = low_tri_tau_vals;
    tau = tau.t() + tau;
    for (int i=0; i<p; i++){
      for (int j=0; j<(p-1); j++){
        tauI(j) = tau(perms.col(i)(j), i);
        Sigma12(j) = Sigma(perms.col(i)(j), i);
        S_sub(j) = S(perms.col(i)(j), i);
      }

      for (int j=0; j<(p-1); j++){
        for (int k=0; k<(p-1); k++){
          Sigma11(j,k) = Sigma(perms.col(i)(j), perms.col(i)(k));
          Omega11(j,k) = Omega(perms.col(i)(j), perms.col(i)(k));
        }
      }

      Omega11inv =Sigma11 - (Sigma12*Sigma12.t())/Sigma(i,i);
      Ci = (S(i,i) + 1.0)*Omega11inv + diagmat(1/tauI);
      mui = -inv(Ci) * S_sub;
      beta = mvrnormArma(1, mui, inv(Ci));
      for (int j=0; j<(p-1); j++){
        Omega(perms.col(i)(j),i)=Omega(i,perms.col(i)(j))=beta(j);
      }
      NumericVector gamm = Rcpp::rgamma(1, n/2 + 1, 1/((S(i,i) + 1.0)/2));
      Omega.submat(i,i,i,i) = gamm(0) + beta * Omega11inv * beta.t();
      Sigma = inv(Omega);
    }
    if (iter>=burnIn){
      OmegaMatList[idx] = Omega;
      idx += 1;
    }
  }
  return OmegaMatList;
}


