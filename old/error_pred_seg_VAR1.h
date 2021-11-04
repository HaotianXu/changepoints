#ifndef error_pred_seg_VAR1_H
#define error_pred_seg_VAR1_H

double rcpp_error_pred_seg_VAR1(const arma::mat& X_futu, const arma::mat& X_curr, int s, int e, double alpha, double lambda, int delta);

#endif