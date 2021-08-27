//DP_VAR1.cpp
#include <RcppArmadillo.h>
#include <limits>
#include "error_pred_seg_VAR1.h"
#include "DP_VAR1.h"

// [[Rcpp::export]]
Rcpp::List rcpp_DP_VAR1(const arma::mat& X_futu, const arma::mat& X_curr, double alpha, double gamma, double lambda, int delta){
  arma::uword n = X_curr.n_cols;
  double b = 0;
  double dist = 0;
  arma::vec bestvalue = arma::zeros<arma::vec>(n+1);
  arma::vec partition = arma::zeros<arma::vec>(n+1);
  bestvalue(0) = -gamma;
  for(int i = 1; i < n+1; ++i){
    bestvalue(i) = std::numeric_limits<int>::max();
    for(int l = 1; l < i+1; ++l){
      if(i - l > delta){
        dist = rcpp_error_pred_seg_VAR1(X_futu, X_curr, l, i, alpha, lambda, delta);
        b = bestvalue(l-1) + gamma + dist;
      }else{
        b = std::numeric_limits<int>::max();
      }
      if (b < bestvalue(i)){
        bestvalue(i) = b;
        partition(i) = l-1;
      }
    }
  }
  int R = n;
  int L = partition(R);
  while(R > 0){
    R = L;
    L = partition(R);
  }
  return Rcpp::List::create(Rcpp::Named("partition")=partition.subvec(1,n));
}
