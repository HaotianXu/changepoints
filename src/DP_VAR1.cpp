//DP_VAR1.cpp
#include <RcppArmadillo.h>
#include "DP_VAR1.h"
#include "lasso.h"

// [[Rcpp::export]]
Rcpp::List rcpp_error_pred_seg_VAR1(const arma::mat& X_futu, const arma::mat& X_curr, int s, int e, const arma::vec& lambda, int delta, double eps){
  int p = X_curr.n_rows;
  arma::mat tran_hat;
  Rcpp::List lassofit;
  arma::mat X;
  arma::mat X_futu_temp;
  arma::vec y;
  double error = 0;
  if(e - s > 2*delta){
    X = X_curr.cols(s-1, e-1).t();
    X_futu_temp = X_futu.cols(s-1, e-1);
    y = X_futu_temp.row(0);
    lassofit = rcpp_lasso_seq(X, y, lambda/sqrt(e-s), eps);
    tran_hat = Rcpp::as<arma::mat>(lassofit["beta_mat"]);
    for(int m = 1; m < p; ++m){
      y = X_futu_temp.row(m);
      lassofit = rcpp_lasso_seq(X, y, lambda/sqrt(e-s), eps);
      tran_hat = arma::join_rows(tran_hat, Rcpp::as<arma::mat>(lassofit["beta_mat"]));
    }
    error = norm(tran_hat.t() * X_curr.cols(s-1, e-1) - X_futu.cols(s-1, e-1), "fro");
  }else{
    error = R_PosInf;
  }
  return Rcpp::List::create(Rcpp::Named("MSE")= error * error,
                            Rcpp::Named("tran_hat")= tran_hat.t());
}



// [[Rcpp::export]]
Rcpp::List rcpp_DP_VAR1(const arma::mat& X_futu, const arma::mat& X_curr, double gamma, const arma::vec& lambda, int delta, double eps){
  int n = X_futu.n_cols;
  double b = 0;
  arma::vec bestvalue = arma::zeros<arma::vec>(n+1);
  arma::vec partition = arma::zeros<arma::vec>(n+1);
  bestvalue(0) = -gamma;
  for(int i = 1; i < n+1; ++i){
    //Rcpp::Rcout << "r is" << std::endl << i << std::endl;
    bestvalue(i) = R_PosInf;
    for(int l = 1; l < i+1; ++l){
      //Rcpp::Rcout << "l is" << std::endl << l << std::endl;
      b = bestvalue(l-1) + gamma + Rcpp::as<double>(rcpp_error_pred_seg_VAR1(X_futu, X_curr, l, i, lambda, delta, eps)["MSE"]);
      if (b < bestvalue(i)){
        bestvalue(i) = b;
        partition(i) = l-1;
        //Rcpp::Rcout << partition(i) << std::endl;
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
