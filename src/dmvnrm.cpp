// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include "dmvnrm.h"

static double const log2pi = std::log(2.0 * M_PI);

/* C++ version of the dtrmv BLAS function */
void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;
  
  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

// [[Rcpp::export]]
arma::vec dmvnrm_arma_fast(arma::mat const &x,  
                           arma::rowvec const &mean,  
                           arma::mat const &sigma) { 
  using arma::uword;
  uword const n = x.n_rows, 
  xdim = x.n_cols;
  arma::vec out(n);
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag())), 
  constants = -(double)xdim/2.0 * log2pi, 
  other_terms = rootisum + constants;
  
  arma::rowvec z;
  for (uword i = 0; i < n; i++) {
    z = (x.row(i) - mean);
    inplace_tri_mat_mult(z, rooti);
    out(i) = other_terms - 0.5 * arma::dot(z, z);     
  }  
  return exp(out);
}


// [[Rcpp::export]]
arma::vec rcpp_dmvnorm_mixt(arma::mat const &evalpoints, arma::mat const &mus, arma::mat const &sigmas, arma::rowvec const &props){
  int p = evalpoints.n_cols;
  int n = evalpoints.n_rows;
  int k = props.size();
  arma::vec dens(n,arma::fill::zeros);
  arma::vec y(n,arma::fill::zeros);
  for(int i = 0; i < k; ++i){
    y = props(i) * dmvnrm_arma_fast(evalpoints, mus.row(i), sigmas.rows(i*p, (i+1)*p-1));
    dens += y;
  }
  return dens;
}