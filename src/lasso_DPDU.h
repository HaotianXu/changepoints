#ifndef lasso_DPDU_H
#define lasso_DPDU_H

double rcpp_lassoDPDU_standardized_obj(const arma::mat& Mtilde, const arma::colvec& Vtilde, const arma::colvec& beta, int n, double lambda);
arma::colvec rcpp_lassoDPDU_standardized(const arma::mat& Mtilde, const arma::colvec& Vtilde, const arma::colvec& beta_start, int n, double lambda, double eps);
arma::mat rcpp_lassoDPDU_standardized_seq(const arma::mat& Mtilde, const arma::colvec& Vtilde, int n, const arma::colvec& lambda_seq, double eps);
Rcpp::List rcpp_lassoDPDU(const arma::mat& Mtilde, const arma::colvec& Vtilde, const arma::vec& Xmeans, const double& Ymean, const arma::vec& weights, int n, double lambda, double eps);
Rcpp::List rcpp_DPDU_regression(const arma::vec& y, const arma::mat& X, double lambda, int zeta, double eps);

#endif