// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <limits>
#include "error_pred_seg_VAR1.h"

// [[Rcpp::export]]
double rcpp_error_pred_seg_VAR1(const arma::mat& X_futu, const arma::mat& X_curr, int s, int e, double alpha, double lambda, int delta){
  int p = X_curr.n_rows;
  arma::sp_mat tran_hat;
  arma::mat X;
  arma::vec y;
  double error = 0;
  Rcpp::Environment glmnet("package:glmnet");
  Rcpp::Function lasso = glmnet["glmnet"];
  Rcpp::Environment stats("package:stats");
  Rcpp::Function coef = stats["coef"];
  if(e - s > 2*delta){
    X = X_curr.cols(s-1, e-1).t();
    y = X.col(0);
    tran_hat = Rcpp::as<arma::sp_mat>(coef(lasso(X, y, Rcpp::Named("alpha", alpha), Rcpp::Named("lambda", lambda))));
    for(int m = 1; m < p; ++m){
      y = X.col(m);
      tran_hat = arma::join_rows(tran_hat, Rcpp::as<arma::sp_mat>(coef(lasso(X, y, Rcpp::Named("alpha", alpha), Rcpp::Named("lambda", lambda)))));
    }
    error = norm(tran_hat.rows(1, p).t() * X_curr.cols(s-1, e-1) - X_futu.cols(s-1, e-1), "fro");
  }else{
    error = std::numeric_limits<int>::max();
  }
  return error * error;
}
