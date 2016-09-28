#ifndef alR_H
#define alR_H

double bw(const arma::vec& x, const int& type);

Rcpp::NumericVector kdeGaussInt2(Rcpp::NumericVector mu, double h, Rcpp::NumericVector q1, Rcpp::NumericVector q2, bool quantile);

#endif