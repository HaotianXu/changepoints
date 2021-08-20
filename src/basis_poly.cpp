//D_P.cpp
#include <RcppArmadillo.h>
#include <limits>
#include "basis_poly.h"

// [[Rcpp::export]]
arma::mat rcpp_basis_poly(int n, int s, int e, int r){
  arma::mat basis_mat(e-s+1, r+1, arma::fill::zeros);
  for(int i = 0; i < r+1; ++i){
    basis_mat.col(i) = arma::pow(arma::linspace(s, e, e-s+1)/n, i);
  }
  return basis_mat;
}
