#ifndef lasso_H
#define lasso_H

double rcpp_soft_threshold_scalar(double x, double lambda);
double rcpp_lasso_standardized_obj(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& beta, double lambda);
arma::colvec rcpp_lasso_standardized(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& beta_start, double lambda, double eps);
arma::mat rcpp_lasso_standardized_seq(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& lambda_seq, double eps);
Rcpp::List rcpp_standardizeXY(const arma::mat& X, const arma::colvec& Y);
Rcpp::List rcpp_lasso_seq(const arma::mat& X, const arma::colvec& Y, const arma::colvec& lambda_seq, double eps);

#endif