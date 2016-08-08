#include <RcppArmadillo.h>
using namespace Rcpp;
#include "alR.h"
// [[Rcpp::depends(RcppArmadillo)]]

//' Objective function for KDE arc length matching.
//'
//' The arc lengths over specified intervals, in the domain of kernel density estimates, are matched.
//'
//' @param beta Vector of regression parameters.
//' @param gamma Design matrix.
//' @param aly A numerical vector of arc length segments to be matched to (same length as \code{beta}).
//' @param q1 A vector of points (not quantiles) specifying the lower limit of the arc length segments.
//' @param q2 A vector of points (not quantiles) specifying the upper limit of the arc length segments.
//'
//' @return Sum of squared differences between \code{gamma}*\code{beta} and \code{aly}.
//'
//' @export
// [[Rcpp::export]]
double alrKDE(const arma::vec& beta, const arma::mat& gamma, const arma::vec& aly, const arma::vec& q1, const arma::vec& q2)
{
arma::vec X = gamma*beta;
double h_gamma = Silverman(X);
arma::vec al_gamma = as<NumericVector>(wrap(kdeGaussInt2(as<NumericVector>(wrap(X)), h_gamma, as<NumericVector>(wrap(q1)), as<NumericVector>(wrap(q2)), false)));
return sum(pow(al_gamma-aly, 2));
}