//trunc_mean.cpp
#include <RcppArmadillo.h>
#include "huber_mean.h"

// [[Rcpp::export]]
double rcpp_huber_mean(const arma::vec& x, double tau) {
  double eps = 1e-8;
  int n = x.size();
  double mu_new = mean(x);
  double mu_old = 0;
  arma::vec r;
  arma::vec a = arma::zeros<arma::vec>(n);
  arma::vec b = arma::ones<arma::vec>(n);
  arma::vec c = arma::zeros<arma::vec>(n);
  while (fabs(mu_new-mu_old) > eps){
    mu_old = mu_new;
    r = x - mu_new;
    for(int i = 0; i < n; ++i){
      if(fabs(r[i]) > tau){
        a[i] = 1;
        b[i] = tau/fabs(r[i]);
      }
      c[i] = x[i]*b[i];
    }
    mu_new = sum(c)/sum(b);
  }
  return mu_new;
}
