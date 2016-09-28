#include <RcppArmadillo.h>
using namespace Rcpp;
#include "alR.h"
// [[Rcpp::depends(RcppArmadillo)]] 

//' Objective function for KDE moment matching.
//'
//' The moments of kernel density estimates are matched.
//'
//' @param beta Vector of regression parameters.
//' @param gamma Design matrix.
//' @param momy A numerical vector of moments to be matched to (same length as \code{beta}).
//' @param kdeGaussMom A R function object to be passed that calculates exact moments.
//' @param type An integer specifying the bandwidth selection method used, see \code{\link{bw}}.
//'
//' @return Sum of squared differences between \code{gamma}*\code{beta} and \code{momy}.
//'
//' @export
// [[Rcpp::export]]
double momKDE(const arma::vec& beta, const arma::mat& gamma, const arma::vec& momy, Function kdeGaussMom, const int& type)
{
arma::vec X = gamma*beta;
double h_gamma = bw(X, type);
arma::vec mom_gamma = as<NumericVector>(wrap(kdeGaussMom(beta.n_rows, X, h_gamma)));
return arma::sum(pow(mom_gamma-momy, 2));
}