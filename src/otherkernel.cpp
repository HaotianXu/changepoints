// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include "otherkernel.h"

// [[Rcpp::export]]
arma::vec epanker(arma::mat const &x,  
                   arma::rowvec const &mean,  
                   double const &bw) { 
  using arma::uword;
  uword const n = x.n_rows, 
    p = x.n_cols;
  arma::vec out(n);
  arma::rowvec z;
  double c_p = tgamma((p+1)/2) / (pow(M_PI, p/2) * tgamma(p/2 + 1));
  for (uword i = 0; i < n; i++) {
    z = (x.row(i) - mean)/bw;
    double d = arma::dot(z, z);
    if (d > 1.0){
      out(i) = 0;
    }else{
      out(i) = c_p * (1.0 - d);
    }
  }  
  return out;
}


// [[Rcpp::export]]
arma::vec rcpp_epanker_mixt(arma::mat const &evalpoints, arma::mat const &mus, double const &bw, arma::rowvec const &props){
  int n = evalpoints.n_rows;
  int k = props.size();
  arma::vec dens(n,arma::fill::zeros);
  arma::vec y(n,arma::fill::zeros);
  for(int i = 0; i < k; ++i){
    y = props(i) * epanker(evalpoints, mus.row(i), bw);
    dens += y;
  }
  return dens;
}


// [[Rcpp::export]]
arma::vec biweiker(arma::mat const &x,  
                           arma::rowvec const &mean,  
                           double const &bw) { 
  using arma::uword;
  uword const n = x.n_rows, 
    p = x.n_cols;
  arma::vec out(n);
  arma::rowvec z;
  double c_p = pow(15, -p/2) / pow(M_PI, p/2);
  for (uword i = 0; i < n; i++) {
    z = (x.row(i) - mean)/bw;
    double d = arma::dot(z, z);
    if (d > 1.0){
      out(i) = 0;
    }else{
      out(i) = c_p * (1.0 - d)*(1.0 - d);
    }
  }  
  return out;
}


// [[Rcpp::export]]
arma::vec rcpp_biweiker_mixt(arma::mat const &evalpoints, arma::mat const &mus, double const &bw, arma::rowvec const &props){
  int n = evalpoints.n_rows;
  int k = props.size();
  arma::vec dens(n,arma::fill::zeros);
  arma::vec y(n,arma::fill::zeros);
  for(int i = 0; i < k; ++i){
    y = props(i) * biweiker(evalpoints, mus.row(i), bw);
    dens += y;
  }
  return dens;
}

// [[Rcpp::export]]
arma::vec triweiker(arma::mat const &x,  
                   arma::rowvec const &mean,  
                   double const &bw) { 
  using arma::uword;
  uword const n = x.n_rows, 
    p = x.n_cols;
  arma::vec out(n);
  arma::rowvec z;
  for (uword i = 0; i < n; i++) {
    z = (x.row(i) - mean)/bw;
    double d = arma::dot(z, z);
    if (d > 1.0){
      out(i) = 0;
    }else{
      double c_p = (p+2) / (2 * pow(bw, p+1) * tgamma(p/2 + 1));
      out(i) = c_p * (1.0 - d)*(1.0 - d)*(1.0 - d);
    }
  }  
  return out;
}


// [[Rcpp::export]]
arma::vec rcpp_triweiker_mixt(arma::mat const &evalpoints, arma::mat const &mus, double const &bw, arma::rowvec const &props){
  int n = evalpoints.n_rows;
  int k = props.size();
  arma::vec dens(n,arma::fill::zeros);
  arma::vec y(n,arma::fill::zeros);
  for(int i = 0; i < k; ++i){
    y = props(i) * triweiker(evalpoints, mus.row(i), bw);
    dens += y;
  }
  return dens;
}