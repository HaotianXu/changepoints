#include <RcppArmadillo.h>
#include "lasso_DPDU.h"
#include "lasso.h"


// [[Rcpp::export]]
double rcpp_lassoDPDU_standardized_obj(const arma::mat& Mtilde, const arma::colvec& Vtilde, const arma::colvec& beta, int n, double lambda){
  // Calculate and return objective value
  double obj = (arma::as_scalar(beta.t()*Mtilde*beta)-2*accu(Vtilde%beta))/n + lambda * accu(abs(beta));
  return obj;
}

// [[Rcpp::export]]
arma::colvec rcpp_lassoDPDU_standardized(const arma::mat& Mtilde, const arma::colvec& Vtilde, const arma::colvec& beta_start, int n, double lambda, double eps){
  // Create beta_last and beta_new vectors
  int p = Mtilde.n_cols;
  arma::colvec beta_last = beta_start;
  arma::colvec beta_new = beta_start;
  // Loss difference variable
  double loss_diff = 100;
  //double err_update = 0.0;
  // Update beta
  while (loss_diff >= eps) {
    // Rcpp::Rcout << "loss diff" << std::endl << loss_diff << std::endl;
    // Rcpp::Rcout << "beta is" << std::endl << beta_new << std::endl;
    beta_last = beta_new;
    // Objective value old
    double loss_old = rcpp_lassoDPDU_standardized_obj(Mtilde, Vtilde, beta_last, n, lambda);
    for(int i = 0; i < p; i++){
      // Update beta
      beta_new(i) = rcpp_soft_threshold_scalar(arma::as_scalar(beta_last(i) + (Vtilde(i) - Mtilde.row(i)*beta_new) / n), lambda);
      // Rcpp::Rcout << "r is" << std::endl << beta_new(i) << std::endl;
    }
    // Objective value new
    double loss_new = rcpp_lassoDPDU_standardized_obj(Mtilde, Vtilde, beta_new, n, lambda);
    
    // Objective value difference
    loss_diff = loss_old - loss_new;
  }
  // Return updated beta
  return beta_new;
}



// [[Rcpp::export]]
arma::mat rcpp_lassoDPDU_standardized_seq(const arma::mat& Mtilde, const arma::colvec& Vtilde, int n, const arma::colvec& lambda_seq, double eps){
  // Creates beta matrix 
  int row = Mtilde.n_cols;
  int col = lambda_seq.size();
  arma::mat beta_mat(row, col);
  // Initialize beta start
  arma::colvec beta_start(row,arma::fill::zeros);
  // Calculate beta for each lambda
  for(int i = 0; i < col; i++){
    arma::colvec beta_new = rcpp_lassoDPDU_standardized(Mtilde, Vtilde, beta_start, n, lambda_seq(i), eps);
    beta_mat.col(i) = beta_new;
    beta_start = beta_new;
  }
  // Return beta matrix
  return beta_mat;
}



// [[Rcpp::export]]
Rcpp::List rcpp_lassoDPDU(const arma::mat& Mtilde, const arma::colvec& Vtilde, const arma::vec& Xmeans, const double& Ymean, const arma::vec& weights, const arma::colvec& beta_start, int n, double lambda, double eps){
  arma::vec lasso_fit = rcpp_lassoDPDU_standardized(Mtilde, Vtilde, beta_start, n, lambda, eps);
  double loss = arma::as_scalar(lasso_fit.t()*Mtilde*lasso_fit-2*accu(Vtilde%lasso_fit));
  arma::vec beta_hat = diagmat(1/weights) * lasso_fit;
  double beta0 = Ymean - arma::as_scalar(Xmeans.t() * beta_hat);
  // Return beta matrix
  return Rcpp::List::create(Rcpp::Named("lasso_fit")=lasso_fit,
                            Rcpp::Named("loss")=loss,
                            Rcpp::Named("beta_hat")=beta_hat,
                            Rcpp::Named("beta0")=beta0);
}


// [[Rcpp::export]]
Rcpp::List rcpp_DPDU_regression(const arma::vec& y, const arma::mat& X, double lambda, int zeta, double eps){
  int N = y.n_elem;
  int p = X.n_cols;
  double b = 0;
  arma::vec bestvalue = arma::zeros<arma::vec>(N+1);
  arma::vec partition = arma::zeros<arma::vec>(N+1);
  bestvalue(0) = -zeta;
  arma::cube M_new(p,p,N,arma::fill::zeros);
  arma::cube M_old(p,p,N,arma::fill::zeros);
  arma::mat V_new(N,p,arma::fill::zeros);
  arma::mat V_old(N,p,arma::fill::zeros);
  arma::vec Ymean_new(N,arma::fill::zeros);
  arma::vec Ymean_old(N,arma::fill::zeros);
  arma::mat Xmeans_new(N,p,arma::fill::zeros);
  arma::mat Xmeans_old(N,p,arma::fill::zeros);
  arma::mat Mcentered(p,p,arma::fill::zeros);
  arma::vec weights(p,arma::fill::zeros);
  arma::mat Mtilde(p,p,arma::fill::zeros);
  arma::vec Vtilde(p,arma::fill::zeros);
  Rcpp::List lassofit;
  arma::mat beta_mat_temp(p+1,N,arma::fill::zeros);
  arma::mat beta_mat_temp2(p+1,N,arma::fill::zeros);
  arma::mat beta_mat(p+1,N,arma::fill::zeros);
  arma::colvec beta_start(p,arma::fill::zeros);
  int n = 0;
  for(int i = 1; i < N+1; ++i){
    bestvalue(i) = R_PosInf;
    for(int l = 1; l < i+1; ++l){
      n = i-l+1;
      Ymean_new(l-1) = (Ymean_old(l-1)*(n-1) + y(i-1))/n;
      Xmeans_new.row(l-1) = (Xmeans_old.row(l-1)*(n-1) + X.row(i-1))/n;
      M_new.slice(l-1) = M_old.slice(l-1) + X.row(i-1).t()*X.row(i-1);
      V_new.row(l-1) = V_old.row(l-1) + y(i-1)*X.row(i-1);
      Mcentered = M_new.slice(l-1) - n*(Xmeans_new.row(l-1).t()*Xmeans_new.row(l-1));
      weights = sqrt(Mcentered.diag()/n);
      Mtilde = diagmat(1/weights)*Mcentered*diagmat(1/weights);
      Vtilde = (V_new.row(l-1) - n*Ymean_new(l-1)*Xmeans_new.row(l-1)).t()/weights;
      if(n >= zeta){
        lassofit = rcpp_lassoDPDU(Mtilde, Vtilde, Xmeans_new.row(l-1).t(), Ymean_new(l-1), weights, beta_start, n, lambda*sqrt(std::max(log(std::max(N,p)), (n-1.0)))*sqrt(log(std::max(N,p)))/(n), eps);
        beta_start = Rcpp::as<arma::vec>(lassofit["lasso_fit"]);
        beta_mat_temp(0,l-1) = Rcpp::as<double>(lassofit["beta0"]);
        beta_mat_temp(arma::span(1,p),l-1) = Rcpp::as<arma::vec>(lassofit["beta_hat"]);
        b = bestvalue(l-1) + zeta + Rcpp::as<double>(lassofit["loss"]);
      }else{
        b = 0;
      }
      if (b < bestvalue(i)){
        bestvalue(i) = b;
        partition(i) = l-1;
        beta_mat_temp2.col(i-1) = beta_mat_temp.col(l-1);
      }
    }
    M_old = M_new;
    V_old = V_new;
    Ymean_old = Ymean_new;
    Xmeans_old = Xmeans_new;
  }
  int R = N;
  int L = partition(R);
  while(R > 0){
    for(int t = L; t < R; ++t){
      beta_mat.col(t) = beta_mat_temp2.col(R-1);
    }
    R = L;
    L = partition(R);
  }
  return Rcpp::List::create(Rcpp::Named("beta_mat")=beta_mat,
                            Rcpp::Named("partition")=partition.subvec(1,N));
}




// [[Rcpp::export]]
double rcpp_lassoDPDU_error(const arma::vec& y, const arma::mat& X, const arma::colvec& beta_hat){
  double error = accu(square(y - X * beta_hat));
  return error;
}





