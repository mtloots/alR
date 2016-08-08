#ifndef alR_H
#define alR_H

Rcpp::NumericVector kdeGaussInt2(Rcpp::NumericVector mu, double h, Rcpp::NumericVector q1, Rcpp::NumericVector q2, bool quantile);

double Silverman(const arma::vec& x);

#endif