//DP_VAR1.cpp
#include <RcppArmadillo.h>
#include "DP_regression.h"
#include "lasso.h"

// [[Rcpp::export]]
Rcpp::List rcpp_error_pred_seg_regression(const arma::vec& y, const arma::mat& X, int s, int e, const arma::vec& lambda, int delta, double eps){
  int p = X.n_cols;
  int n = X.n_rows;
  arma::vec beta_hat;
  Rcpp::List lassofit;
  arma::mat X_temp;
  arma::vec y_temp;
  double error = 0.0;
  if(e - s > 2*delta){
    X_temp = X.rows(s-1, e-1);
    y_temp = y.subvec(s-1, e-1);
    lassofit = rcpp_lasso_seq(X_temp, y_temp, lambda*sqrt(std::max(log(std::max(n,p)), (e-s+0.0)))*sqrt(log(std::max(n,p)))/(e-s), eps);
    beta_hat = Rcpp::as<arma::vec>(lassofit["beta_mat"]);
    error = accu(square(y_temp - X_temp * beta_hat));
  }else{
    error = R_PosInf;
  }
  return Rcpp::List::create(Rcpp::Named("MSE")= error,
                            Rcpp::Named("beta_hat")= beta_hat);
}



// [[Rcpp::export]]
Rcpp::List rcpp_DP_regression(const arma::vec& y, const arma::mat& X, double gamma, const arma::vec& lambda, int delta, double eps){
  int n = y.n_elem;
  int p = X.n_cols;
  double b = 0;
  arma::vec bestvalue = arma::zeros<arma::vec>(n+1);
  arma::vec partition = arma::zeros<arma::vec>(n+1);
  bestvalue(0) = -gamma*log(std::max(n,p));
  for(int i = 1; i < n+1; ++i){
    //Rcpp::Rcout << "r is" << std::endl << i << std::endl;
    bestvalue(i) = R_PosInf;
    for(int l = 1; l < i+1; ++l){
      //Rcpp::Rcout << "l is" << std::endl << l << std::endl;
      b = bestvalue(l-1) + gamma*log(std::max(n,p)) + Rcpp::as<double>(rcpp_error_pred_seg_regression(y, X, l, i, lambda, delta, eps)["MSE"]);
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
