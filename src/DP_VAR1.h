#ifndef DP_VAR1_H
#define DP_VAR1_H

Rcpp::List rcpp_DP_VAR1(const arma::mat& X_futu, const arma::mat& X_curr, double gamma, const arma::vec& lambda, int delta, double eps);
double rcpp_error_pred_seg_VAR1(const arma::mat& X_futu, const arma::mat& X_curr, int s, int e, const arma::vec& lambda, int delta, double eps);

#endif
