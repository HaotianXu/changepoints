#include <RcppArmadillo.h>
#include "lasso.h"


// [[Rcpp::export]]
double rcpp_soft_threshold_scalar(double x, double lambda){
  if(x > lambda){
    return(x-lambda);
  }else if(x < -lambda){
    return(x + lambda);
  }else{
    return(0);
  }
}


// [[Rcpp::export]]
double rcpp_lasso_standardized_obj(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& beta, double lambda){
  int n = Xtilde.n_rows;
  // Initialize parameters
  arma::mat X = Xtilde;
  arma::colvec Y = Ytilde;
  // Calculate and return objective value
  double obj = accu(square(Y - X * beta)/(2 * n)) + lambda * accu(abs(beta));
  return obj;
}

// [[Rcpp::export]]
arma::colvec rcpp_lasso_standardized(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& beta_start, double lambda, double eps){
  // Create beta_last and beta_new vectors
  int n = Xtilde.n_rows;
  int p = Xtilde.n_cols;
  arma::colvec beta_last = beta_start;
  arma::colvec beta_new = beta_start;
  // Loss difference variable
  double loss_diff = 100;
  // Full initial residual
  arma::colvec res_vec = Ytilde - Xtilde * beta_start;
  // Update beta
  while (loss_diff >= eps) {
    beta_last = beta_new;
    // Objective value old
    double loss_old = rcpp_lasso_standardized_obj(Xtilde, Ytilde, beta_last, lambda);
    for(int i = 0; i < p; i++){
      // Update beta
      beta_new(i) = rcpp_soft_threshold_scalar(arma::as_scalar(beta_last(i) + (Xtilde.col(i).t() * res_vec / n)), lambda);
      // Partial residual
      res_vec = res_vec + Xtilde.col(i) * (beta_last(i) - beta_new(i));
    }
    // Objective value new
    double loss_new = rcpp_lasso_standardized_obj(Xtilde, Ytilde, beta_new, lambda);
    
    // Objective value difference
    loss_diff = loss_old - loss_new;
  }
  // Return updated beta
  return beta_new;
}  

// [[Rcpp::export]]
arma::mat rcpp_lasso_standardized_seq(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& lambda_seq, double eps){
  // Creates beta matrix 
  int row = Xtilde.n_cols;
  int col = lambda_seq.size();
  arma::mat beta_mat(row, col);
  // Initialize beta start
  arma::colvec beta_start(row);
  beta_start.fill(0.0);
  // Calculate beta for each lambda
  for(int i = 0; i < col; i++){
    arma::colvec beta_new = rcpp_lasso_standardized(Xtilde, Ytilde, beta_start, lambda_seq(i), eps);
    beta_mat.col(i) = beta_new;
    beta_start = beta_new;
  }
  // Return beta matrix
  return beta_mat;
}


// [[Rcpp::export]]
Rcpp::List rcpp_standardizeXY(const arma::mat& X, const arma::colvec& Y){
  double Ymean = mean(Y);
  arma::vec Ytilde = Y - Ymean;
  arma::vec Xmeans = arma::conv_to<arma::vec>::from(mean(X, 0));
  arma::mat X_centred = X.each_row() - Xmeans.t();
  arma::vec weights = sqrt(arma::conv_to<arma::vec>::from(sum(X_centred % X_centred, 0))/X_centred.n_rows);
  arma::mat Xtilde = X_centred * diagmat(1/weights);
  // Xtilde - centered and appropriately scaled X
  // Ytilde - centered Y
  // Ymean - the mean of original Y
  // Xmeans - means of columns of X (vector)
  // weights - defined as sqrt(X_j^{\top}X_j/n) after centering of X but before scaling
  return Rcpp::List::create(Rcpp::Named("Ymean")=Ymean,
                            Rcpp::Named("Ytilde")=Ytilde,
                            Rcpp::Named("Xmeans")=Xmeans,
                            Rcpp::Named("Xtilde")=Xtilde,
                            Rcpp::Named("weights")=weights);
}



// [[Rcpp::export]]
Rcpp::List rcpp_lasso_seq(const arma::mat& X, const arma::colvec& Y, const arma::colvec& lambda_seq, double eps){
  Rcpp::List standardized = rcpp_standardizeXY(X, Y);
  arma::mat Xtilde = standardized["Xtilde"];
  arma::vec Ytilde = standardized["Ytilde"];
  arma::vec weights = standardized["weights"];
  arma::vec Xmeans = standardized["Xmeans"];
  double Ymean = standardized["Ymean"];
  arma::mat lasso_fit = rcpp_lasso_standardized_seq(Xtilde, Ytilde, lambda_seq, eps);
  
  arma::mat beta_mat = diagmat(1/weights) * lasso_fit;
  arma::vec beta0_vec = Ymean - arma::conv_to<arma::vec>::from(Xmeans.t() * beta_mat);
  // Return beta matrix
  return Rcpp::List::create(Rcpp::Named("lambda_seq")=lambda_seq,
                            Rcpp::Named("beta_mat")=beta_mat,
                            Rcpp::Named("beta0_vec")=beta0_vec);
}


