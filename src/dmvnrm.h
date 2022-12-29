#ifndef dmvnrm_H
#define dmvnrm_H

arma::vec dmvnrm_arma_fast(arma::mat const &x, arma::rowvec const &mean, arma::mat const &sigma);
arma::vec rcpp_dmvnorm_mixt(arma::mat const &evalpoints, arma::mat const &mus, arma::mat const &sigmas, arma::rowvec const &props);

#endif