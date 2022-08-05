// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rcpp_error_pred_seg_VAR1
Rcpp::List rcpp_error_pred_seg_VAR1(const arma::mat& X_futu, const arma::mat& X_curr, int s, int e, const arma::vec& lambda, int delta, double eps);
RcppExport SEXP _changepoints_rcpp_error_pred_seg_VAR1(SEXP X_futuSEXP, SEXP X_currSEXP, SEXP sSEXP, SEXP eSEXP, SEXP lambdaSEXP, SEXP deltaSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X_futu(X_futuSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X_curr(X_currSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type e(eSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_error_pred_seg_VAR1(X_futu, X_curr, s, e, lambda, delta, eps));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_DP_VAR1
Rcpp::List rcpp_DP_VAR1(const arma::mat& X_futu, const arma::mat& X_curr, double gamma, const arma::vec& lambda, int delta, double eps);
RcppExport SEXP _changepoints_rcpp_DP_VAR1(SEXP X_futuSEXP, SEXP X_currSEXP, SEXP gammaSEXP, SEXP lambdaSEXP, SEXP deltaSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X_futu(X_futuSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X_curr(X_currSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_DP_VAR1(X_futu, X_curr, gamma, lambda, delta, eps));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_DP_poly
Rcpp::List rcpp_DP_poly(const arma::vec& y, int r, double gamma, int delta);
RcppExport SEXP _changepoints_rcpp_DP_poly(SEXP ySEXP, SEXP rSEXP, SEXP gammaSEXP, SEXP deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< int >::type delta(deltaSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_DP_poly(y, r, gamma, delta));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_error_pred_seg_regression
Rcpp::List rcpp_error_pred_seg_regression(const arma::vec& y, const arma::mat& X, int s, int e, const arma::vec& lambda, int delta, double eps);
RcppExport SEXP _changepoints_rcpp_error_pred_seg_regression(SEXP ySEXP, SEXP XSEXP, SEXP sSEXP, SEXP eSEXP, SEXP lambdaSEXP, SEXP deltaSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type e(eSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_error_pred_seg_regression(y, X, s, e, lambda, delta, eps));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_DP_regression
Rcpp::List rcpp_DP_regression(const arma::vec& y, const arma::mat& X, double gamma, const arma::vec& lambda, int delta, double eps);
RcppExport SEXP _changepoints_rcpp_DP_regression(SEXP ySEXP, SEXP XSEXP, SEXP gammaSEXP, SEXP lambdaSEXP, SEXP deltaSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_DP_regression(y, X, gamma, lambda, delta, eps));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_DP_univar
Rcpp::List rcpp_DP_univar(const arma::vec& y, double gamma, int delta);
RcppExport SEXP _changepoints_rcpp_DP_univar(SEXP ySEXP, SEXP gammaSEXP, SEXP deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< int >::type delta(deltaSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_DP_univar(y, gamma, delta));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_basis_poly
arma::mat rcpp_basis_poly(int n, int s, int e, int r);
RcppExport SEXP _changepoints_rcpp_basis_poly(SEXP nSEXP, SEXP sSEXP, SEXP eSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type e(eSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_basis_poly(n, s, e, r));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_huber_mean
double rcpp_huber_mean(const arma::vec& x, double tau);
RcppExport SEXP _changepoints_rcpp_huber_mean(SEXP xSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_huber_mean(x, tau));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_soft_threshold_scalar
double rcpp_soft_threshold_scalar(double x, double lambda);
RcppExport SEXP _changepoints_rcpp_soft_threshold_scalar(SEXP xSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_soft_threshold_scalar(x, lambda));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_lasso_standardized_obj
double rcpp_lasso_standardized_obj(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& beta, double lambda);
RcppExport SEXP _changepoints_rcpp_lasso_standardized_obj(SEXP XtildeSEXP, SEXP YtildeSEXP, SEXP betaSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Xtilde(XtildeSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type Ytilde(YtildeSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_lasso_standardized_obj(Xtilde, Ytilde, beta, lambda));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_lasso_standardized
arma::colvec rcpp_lasso_standardized(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& beta_start, double lambda, double eps);
RcppExport SEXP _changepoints_rcpp_lasso_standardized(SEXP XtildeSEXP, SEXP YtildeSEXP, SEXP beta_startSEXP, SEXP lambdaSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Xtilde(XtildeSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type Ytilde(YtildeSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type beta_start(beta_startSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_lasso_standardized(Xtilde, Ytilde, beta_start, lambda, eps));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_lasso_standardized_seq
arma::mat rcpp_lasso_standardized_seq(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& lambda_seq, double eps);
RcppExport SEXP _changepoints_rcpp_lasso_standardized_seq(SEXP XtildeSEXP, SEXP YtildeSEXP, SEXP lambda_seqSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Xtilde(XtildeSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type Ytilde(YtildeSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type lambda_seq(lambda_seqSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_lasso_standardized_seq(Xtilde, Ytilde, lambda_seq, eps));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_standardizeXY
Rcpp::List rcpp_standardizeXY(const arma::mat& X, const arma::colvec& Y);
RcppExport SEXP _changepoints_rcpp_standardizeXY(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_standardizeXY(X, Y));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_lasso_seq
Rcpp::List rcpp_lasso_seq(const arma::mat& X, const arma::colvec& Y, const arma::colvec& lambda_seq, double eps);
RcppExport SEXP _changepoints_rcpp_lasso_seq(SEXP XSEXP, SEXP YSEXP, SEXP lambda_seqSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type lambda_seq(lambda_seqSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_lasso_seq(X, Y, lambda_seq, eps));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_lassoDPDU_standardized_obj
double rcpp_lassoDPDU_standardized_obj(const arma::mat& Mtilde, const arma::colvec& Vtilde, const arma::colvec& beta, int n, double lambda);
RcppExport SEXP _changepoints_rcpp_lassoDPDU_standardized_obj(SEXP MtildeSEXP, SEXP VtildeSEXP, SEXP betaSEXP, SEXP nSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Mtilde(MtildeSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type Vtilde(VtildeSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_lassoDPDU_standardized_obj(Mtilde, Vtilde, beta, n, lambda));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_lassoDPDU_standardized
arma::colvec rcpp_lassoDPDU_standardized(const arma::mat& Mtilde, const arma::colvec& Vtilde, const arma::colvec& beta_start, int n, double lambda, double eps);
RcppExport SEXP _changepoints_rcpp_lassoDPDU_standardized(SEXP MtildeSEXP, SEXP VtildeSEXP, SEXP beta_startSEXP, SEXP nSEXP, SEXP lambdaSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Mtilde(MtildeSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type Vtilde(VtildeSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type beta_start(beta_startSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_lassoDPDU_standardized(Mtilde, Vtilde, beta_start, n, lambda, eps));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_lassoDPDU_standardized_seq
arma::mat rcpp_lassoDPDU_standardized_seq(const arma::mat& Mtilde, const arma::colvec& Vtilde, int n, const arma::colvec& lambda_seq, double eps);
RcppExport SEXP _changepoints_rcpp_lassoDPDU_standardized_seq(SEXP MtildeSEXP, SEXP VtildeSEXP, SEXP nSEXP, SEXP lambda_seqSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Mtilde(MtildeSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type Vtilde(VtildeSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type lambda_seq(lambda_seqSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_lassoDPDU_standardized_seq(Mtilde, Vtilde, n, lambda_seq, eps));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_lassoDPDU
Rcpp::List rcpp_lassoDPDU(const arma::mat& Mtilde, const arma::colvec& Vtilde, const arma::vec& Xmeans, const double& Ymean, const arma::vec& weights, const arma::colvec& beta_start, int n, double lambda, double eps);
RcppExport SEXP _changepoints_rcpp_lassoDPDU(SEXP MtildeSEXP, SEXP VtildeSEXP, SEXP XmeansSEXP, SEXP YmeanSEXP, SEXP weightsSEXP, SEXP beta_startSEXP, SEXP nSEXP, SEXP lambdaSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Mtilde(MtildeSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type Vtilde(VtildeSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Xmeans(XmeansSEXP);
    Rcpp::traits::input_parameter< const double& >::type Ymean(YmeanSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type beta_start(beta_startSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_lassoDPDU(Mtilde, Vtilde, Xmeans, Ymean, weights, beta_start, n, lambda, eps));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_DPDU_regression
Rcpp::List rcpp_DPDU_regression(const arma::vec& y, const arma::mat& X, double lambda, int zeta, double eps);
RcppExport SEXP _changepoints_rcpp_DPDU_regression(SEXP ySEXP, SEXP XSEXP, SEXP lambdaSEXP, SEXP zetaSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type zeta(zetaSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_DPDU_regression(y, X, lambda, zeta, eps));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_lassoDPDU_error
double rcpp_lassoDPDU_error(const arma::vec& y, const arma::mat& X, const arma::colvec& beta_hat);
RcppExport SEXP _changepoints_rcpp_lassoDPDU_error(SEXP ySEXP, SEXP XSEXP, SEXP beta_hatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type beta_hat(beta_hatSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_lassoDPDU_error(y, X, beta_hat));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_soft_threshold
Rcpp::List rcpp_soft_threshold(arma::mat& x_mat, double lambda);
RcppExport SEXP _changepoints_rcpp_soft_threshold(SEXP x_matSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x_mat(x_matSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_soft_threshold(x_mat, lambda));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_soft_impute
Rcpp::List rcpp_soft_impute(const Rcpp::List& data_incomplete_list, const Rcpp::List& eta_list, double lambda, double rho, int it_max);
RcppExport SEXP _changepoints_rcpp_soft_impute(SEXP data_incomplete_listSEXP, SEXP eta_listSEXP, SEXP lambdaSEXP, SEXP rhoSEXP, SEXP it_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type data_incomplete_list(data_incomplete_listSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type eta_list(eta_listSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< int >::type it_max(it_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_soft_impute(data_incomplete_list, eta_list, lambda, rho, it_max));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_lambda
double rcpp_lambda(int s, int e, int t, double alpha, double rho, double m, double d, double C_lambda);
RcppExport SEXP _changepoints_rcpp_lambda(SEXP sSEXP, SEXP eSEXP, SEXP tSEXP, SEXP alphaSEXP, SEXP rhoSEXP, SEXP mSEXP, SEXP dSEXP, SEXP C_lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type e(eSEXP);
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    Rcpp::traits::input_parameter< double >::type C_lambda(C_lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_lambda(s, e, t, alpha, rho, m, d, C_lambda));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_CUSUM
double rcpp_CUSUM(const Rcpp::List& data_incomplete_list, const Rcpp::List& eta_list, int s, int e, int t, double alpha, double rho, double m, double C_lambda, int delta);
RcppExport SEXP _changepoints_rcpp_CUSUM(SEXP data_incomplete_listSEXP, SEXP eta_listSEXP, SEXP sSEXP, SEXP eSEXP, SEXP tSEXP, SEXP alphaSEXP, SEXP rhoSEXP, SEXP mSEXP, SEXP C_lambdaSEXP, SEXP deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type data_incomplete_list(data_incomplete_listSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type eta_list(eta_listSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type e(eSEXP);
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type C_lambda(C_lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type delta(deltaSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_CUSUM(data_incomplete_list, eta_list, s, e, t, alpha, rho, m, C_lambda, delta));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_changepoints_rcpp_error_pred_seg_VAR1", (DL_FUNC) &_changepoints_rcpp_error_pred_seg_VAR1, 7},
    {"_changepoints_rcpp_DP_VAR1", (DL_FUNC) &_changepoints_rcpp_DP_VAR1, 6},
    {"_changepoints_rcpp_DP_poly", (DL_FUNC) &_changepoints_rcpp_DP_poly, 4},
    {"_changepoints_rcpp_error_pred_seg_regression", (DL_FUNC) &_changepoints_rcpp_error_pred_seg_regression, 7},
    {"_changepoints_rcpp_DP_regression", (DL_FUNC) &_changepoints_rcpp_DP_regression, 6},
    {"_changepoints_rcpp_DP_univar", (DL_FUNC) &_changepoints_rcpp_DP_univar, 3},
    {"_changepoints_rcpp_basis_poly", (DL_FUNC) &_changepoints_rcpp_basis_poly, 4},
    {"_changepoints_rcpp_huber_mean", (DL_FUNC) &_changepoints_rcpp_huber_mean, 2},
    {"_changepoints_rcpp_soft_threshold_scalar", (DL_FUNC) &_changepoints_rcpp_soft_threshold_scalar, 2},
    {"_changepoints_rcpp_lasso_standardized_obj", (DL_FUNC) &_changepoints_rcpp_lasso_standardized_obj, 4},
    {"_changepoints_rcpp_lasso_standardized", (DL_FUNC) &_changepoints_rcpp_lasso_standardized, 5},
    {"_changepoints_rcpp_lasso_standardized_seq", (DL_FUNC) &_changepoints_rcpp_lasso_standardized_seq, 4},
    {"_changepoints_rcpp_standardizeXY", (DL_FUNC) &_changepoints_rcpp_standardizeXY, 2},
    {"_changepoints_rcpp_lasso_seq", (DL_FUNC) &_changepoints_rcpp_lasso_seq, 4},
    {"_changepoints_rcpp_lassoDPDU_standardized_obj", (DL_FUNC) &_changepoints_rcpp_lassoDPDU_standardized_obj, 5},
    {"_changepoints_rcpp_lassoDPDU_standardized", (DL_FUNC) &_changepoints_rcpp_lassoDPDU_standardized, 6},
    {"_changepoints_rcpp_lassoDPDU_standardized_seq", (DL_FUNC) &_changepoints_rcpp_lassoDPDU_standardized_seq, 5},
    {"_changepoints_rcpp_lassoDPDU", (DL_FUNC) &_changepoints_rcpp_lassoDPDU, 9},
    {"_changepoints_rcpp_DPDU_regression", (DL_FUNC) &_changepoints_rcpp_DPDU_regression, 5},
    {"_changepoints_rcpp_lassoDPDU_error", (DL_FUNC) &_changepoints_rcpp_lassoDPDU_error, 3},
    {"_changepoints_rcpp_soft_threshold", (DL_FUNC) &_changepoints_rcpp_soft_threshold, 2},
    {"_changepoints_rcpp_soft_impute", (DL_FUNC) &_changepoints_rcpp_soft_impute, 5},
    {"_changepoints_rcpp_lambda", (DL_FUNC) &_changepoints_rcpp_lambda, 8},
    {"_changepoints_rcpp_CUSUM", (DL_FUNC) &_changepoints_rcpp_CUSUM, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_changepoints(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
