//DP_poly.cpp
#include <RcppArmadillo.h>
#include "basis_poly.h"
#include "DP_poly.h"

// [[Rcpp::export]]
Rcpp::List rcpp_DP_poly(const arma::vec& y, int r, double gamma, int delta){
  int n = y.size();
  double b = 0;
  double dist = 0;
  arma::vec bestvalue = arma::zeros<arma::vec>(n+1);
  arma::vec partition = arma::zeros<arma::vec>(n+1);
  arma::vec yhat = arma::zeros<arma::vec>(n+1);
  arma::mat u_mat;
  arma::mat proj_mat;
  bestvalue(0) = -gamma;
  for(int i = 1; i < n+1; ++i){
    bestvalue(i) = R_PosInf;
    for(int l = 1; l < i+1; ++l){
      if(i - l > 2*delta){
        u_mat = rcpp_basis_poly(n, l, i, r);
        proj_mat = u_mat * (u_mat.t() * u_mat).i() * u_mat.t();
        dist = arma::norm(y.subvec(l-1,i-1) - proj_mat * y.subvec(l-1,i-1), 2);
        b = bestvalue(l-1) + gamma + dist*dist;
      }else{
        b = R_PosInf;
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
    u_mat = rcpp_basis_poly(n, L, R-1, r);
    proj_mat = u_mat * (u_mat.t() * u_mat).i() * u_mat.t();
    yhat.subvec(L,R-1) = proj_mat * y.subvec(L,R-1);
    R = L;
    L = partition(R);
  }
  return Rcpp::List::create(Rcpp::Named("partition")=partition.subvec(1,n),
                            Rcpp::Named("yhat")=yhat.subvec(0,n-1));
}
