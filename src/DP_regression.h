#ifndef DP_REGRESSION_H
#define DP_REGRESSION_H

Rcpp::List rcpp_error_pred_seg_regression(const arma::vec& y, const arma::mat& X, int s, int e, const arma::vec& lambda, int delta, double eps);
Rcpp::List rcpp_DP_regression(const arma::vec& y, const arma::mat& X, double gamma, const arma::vec& lambda, int delta, double eps);

#endif
