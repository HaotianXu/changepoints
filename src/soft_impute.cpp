#include <RcppArmadillo.h>
#include "soft_impute.h"

// [[Rcpp::export]]
Rcpp::List rcpp_soft_threshold(arma::mat& x_mat, double lambda){
  arma::mat U;
  arma::vec d;
  arma::mat V;
  svd(U, d, V, x_mat);
  int dim = d.size();
  arma::vec zero_vec = arma::zeros<arma::vec>(dim);
  arma::vec d_thre = max(d-lambda, zero_vec);
  return Rcpp::List::create(Rcpp::Named("u")=U,
                     Rcpp::Named("d")=d_thre,
                     Rcpp::Named("v")=V);
}



// [[Rcpp::export]]
Rcpp::List rcpp_soft_impute(const Rcpp::List& data_incomplete_list, const Rcpp::List& eta_list, double lambda, double rho, int it_max = 10000){
  int dim = Rcpp::as<arma::mat>(data_incomplete_list[0]).n_rows;
  int obs_num = data_incomplete_list.size();
  int it = 0;
  arma::mat MOld_mat = arma::mat(dim, dim, arma::fill::zeros);
  arma::mat MNew_mat = arma::mat(dim, dim, arma::fill::zeros);
  arma::mat diff_mat = arma::mat(dim, dim, arma::fill::zeros);
  arma::mat diff_miss_mat = arma::mat(dim, dim, arma::fill::zeros);
  Rcpp::List temp;
  arma::mat x_mat = arma::mat(dim, dim, arma::fill::zeros);
  arma::mat zero_mat = arma::mat(dim, dim, arma::fill::zeros);
  arma::mat one_mat = arma::mat(dim, dim, arma::fill::ones);
  int FLAG = 0;
  while((FLAG == 0) & (it < it_max)){
    it += 1;
    x_mat = zero_mat;
    diff_miss_mat = zero_mat;
    for(int i = 0; i < obs_num; ++i){
      x_mat += Rcpp::as<arma::mat>(data_incomplete_list[i]) + (one_mat - Rcpp::as<arma::mat>(eta_list[i])) % MOld_mat;
    }
    x_mat /= obs_num;
    temp = rcpp_soft_threshold(x_mat, lambda);
    MNew_mat = Rcpp::as<arma::mat>(temp["u"]) * arma::diagmat(Rcpp::as<arma::vec>(temp["d"])) * Rcpp::as<arma::mat>(temp["v"]).t();
    diff_mat = MNew_mat - MOld_mat;
    for(int i = 0; i < obs_num; ++i){
      diff_miss_mat += (one_mat - Rcpp::as<arma::mat>(eta_list[i])) % diff_mat;
    }
    diff_miss_mat /= obs_num;
    if((arma::norm(diff_miss_mat, 2) < lambda/3) && (abs(diff_mat).max() < rho)){
      FLAG = 1;
    }
    MOld_mat = MNew_mat;
    MOld_mat.elem(arma::find(MOld_mat > rho)).fill(rho);
    MOld_mat.elem(arma::find(MOld_mat < -rho)).fill(-rho);
    // Rcout << "Armadillo matrix is" << std::endl << MOld_mat << std::endl;
  }
  return Rcpp::List::create(Rcpp::Named("u")=temp["u"],
                            Rcpp::Named("d")=temp["d"],
                            Rcpp::Named("v")=temp["v"]);
}


// [[Rcpp::export]]
double rcpp_lambda(int s, int e, int t, double alpha, double rho, double m, double d, double C_lambda){
  double res;
  res = C_lambda*(m*sqrt(d*rho) + sqrt(log(32*t*t/alpha)))/sqrt(e-s+1);
  return(res);
}


// [[Rcpp::export]]
double rcpp_CUSUM(const Rcpp::List& data_incomplete_list, const Rcpp::List& eta_list, int s, int e, int t, double alpha, double rho, double m, double C_lambda, int delta){
  double result;
  int d = Rcpp::as<arma::mat>(data_incomplete_list[0]).n_rows;
  arma::mat result_mat;
  Rcpp::List temp_st;
  Rcpp::List temp_te;
  if((t-s < 2*delta-1) | (e-t < 2*delta)){
    result = 0;
  }else{
    temp_st = rcpp_soft_impute(data_incomplete_list[Rcpp::Range(s-1, t-1)], eta_list[Rcpp::Range(s-1, t-1)], rcpp_lambda(s-1, t-1, e-s+1, alpha, rho, m, d, C_lambda), 1);
    temp_te = rcpp_soft_impute(data_incomplete_list[Rcpp::Range(t, e-1)], eta_list[Rcpp::Range(t, e-1)], rcpp_lambda(t, e-1, e-s+1, alpha, rho, m, d, C_lambda), 1);
    result_mat = Rcpp::as<arma::mat>(temp_st["u"]) * arma::diagmat(Rcpp::as<arma::vec>(temp_st["d"])) * Rcpp::as<arma::mat>(temp_st["v"]).t() - Rcpp::as<arma::mat>(temp_te["u"]) * arma::diagmat(Rcpp::as<arma::vec>(temp_te["d"])) * Rcpp::as<arma::mat>(temp_te["v"]).t();
    result = arma::norm(result_mat, "fro");
  }
  return(result);
}
