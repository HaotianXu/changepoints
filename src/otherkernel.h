#ifndef otherkernel_H
#define otherkernel_H

arma::vec biweiker(arma::mat const &x, arma::rowvec const &mean, double const &bw);
arma::vec rcpp_biweiker_mixt(arma::mat const &evalpoints, arma::mat const &mus, double const &bw, arma::rowvec const &props);
arma::vec triweiker(arma::mat const &x, arma::rowvec const &mean, double const &bw);
arma::vec rcpp_triweiker_mixt(arma::mat const &evalpoints, arma::mat const &mus, double const &bw, arma::rowvec const &props);
arma::vec epanker(arma::mat const &x, arma::rowvec const &mean, double const &bw);
arma::vec rcpp_epanker_mixt(arma::mat const &evalpoints, arma::mat const &mus, double const &bw, arma::rowvec const &props);


#endif