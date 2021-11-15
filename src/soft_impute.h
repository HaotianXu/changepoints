#ifndef SOFT_IMPUTE_H
#define SOFT_IMPUTE_H

Rcpp::List rcpp_soft_threshold(arma::mat& x_mat, double lambda);
Rcpp::List rcpp_soft_impute(const Rcpp::List& data_incomplete_list, const Rcpp::List& eta_list, double lambda, double rho, int it_max);
double rcpp_lambda(int s, int e, int t, double alpha, double rho, double m, double d, double C_lambda);
double rcpp_CUSUM(const Rcpp::List& data_incomplete_list, const Rcpp::List& eta_list, int s, int e, int t, double alpha, double rho, double m, double C_lambda, int delta);

#endif
