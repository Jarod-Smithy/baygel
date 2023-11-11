#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);

  return arma::repmat(mu, 1, n).t() + Y * chol(sigma);
}

arma::vec rrinvgauss(int n, double mu, double lambda){
  arma::vec random_vector(n);
  double z,y,x,u;
  for(int i=0; i<n; ++i){
    z=R::rnorm(0,1);
    y=z*z;
    x=mu+0.5*mu*mu*y/lambda - 0.5*(mu/lambda)*sqrt(4*mu*lambda*y+mu*mu*y*y);
    u=R::runif(0,1);
    if(u <= mu/(mu+x)){
      random_vector(i)=x;
    }else{
      random_vector(i)=mu*mu/x;
    };
  }
  return(random_vector);
}

arma::vec extract_upperind(arma::mat A){
  // Extract upper triangular part and set diagonal to zero
  arma::mat upper = trimatu(A);
  upper.diag().zeros();

  // Get indices where values > 0
  arma::uvec indices = find(vectorise(upper) > 0);

  // Extract values using indices
  arma::vec values = vectorise(upper.elem(indices));

  return values;
}

NumericVector rrgamma(int n, double shape, double scale) {
  return(rgamma(n, shape, scale));
}
