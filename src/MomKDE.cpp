#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]] 
#include "alR.h"

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
//' @return Square root of sum of squared differences between \code{gamma}*\code{beta} and \code{momy} (Eucledian distance).
//'
//' @export
// [[Rcpp::export]]
double momKDE(const arma::vec& beta, const arma::mat& gamma, const arma::vec& momy, Function kdeGaussMom, const int& type)
{
arma::vec X = gamma*beta;
double h_gamma = bw(as<NumericVector>(wrap(X)), type);
arma::vec mom_gamma = as<NumericVector>(wrap(kdeGaussMom(beta.n_rows, X, h_gamma)));
return pow(arma::sum(pow(mom_gamma-momy, 2)), 0.5);
}