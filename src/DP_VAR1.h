#ifndef DP_VAR1_H
#define DP_VAR1_H

Rcpp::List rcpp_DP_VAR1(const arma::mat& X_futu, const arma::mat& X_curr, double alpha, double gamma, double lambda, int delta);

#endif
