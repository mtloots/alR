#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]] 

//' Silverman's rule of thumb.
//'
//' Calculate the bandwidth estimator using Silverman's rule of thumb.
//'
//' @param x A vector of data points.
//'
//' @return Bandwidth estimator based on Silverman's rule of thumb.
//'
//' @examples
//' set.seed(1)
//' x <- rnorm(100)
//' Silverman(x)
//'
//' @export
// [[Rcpp::export]]
double Silverman(const arma::vec& x)
{
return arma::stddev(x)*pow(4.0/(3*x.n_rows), 0.2);
}