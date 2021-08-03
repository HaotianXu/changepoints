//D_P.cpp
#include <RcppArmadillo.h>
#include <limits>
#include "D_P.h"

// [[Rcpp::export]]
Rcpp::List rcpp_D_P(const arma::vec& y, double gamma) {
  int N = y.size();
  double b = 0;
  arma::vec bestvalue = arma::zeros<arma::vec>(N+1);
  arma::vec partition = arma::zeros<arma::vec>(N+1);
  arma::vec yhat = arma::zeros<arma::vec>(N+1);
  bestvalue(0) = -gamma;
  for(int r = 1; r < N+1; ++r){
    bestvalue(r) = std::numeric_limits<int>::max();
    for(int l = 1; l < r+1; ++l){
      b = bestvalue(l-1) + gamma + arma::norm(y.subvec(l-1,r-1) - mean(y.subvec(l-1,r-1)), 2);
      if (b <= bestvalue(r)){
        bestvalue(r) = b;
        partition(r) = l-1;
      }
    }
  }
  int R = N;
  int L = partition(R);
  while(R > 0){
    for(int t = L+1; t < R+1; ++t){
      yhat(t) = mean(y.subvec(L,R-1));
    }
    R = L;
    L = partition(R);
  }
  return Rcpp::List::create(Rcpp::Named("partition")=partition.subvec(1,N),
                            Rcpp::Named("yhat")=yhat.subvec(1,N));
}
