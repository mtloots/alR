#include <Rcpp.h>
using namespace Rcpp;
#include "alR.h"


//' The Generalised Pareto Distribution (GPD).
//'
//' Density, distribution function, quantile function and random generation for the generalised Pareto distribution (GPD) with location parameter \code{mu}, scale parameter \code{sigma}, and shape parameter \code{alpha}. .
//'
//' @param x,q A data point, or quantile, at which the GPD should be evaluated (may be vectors).
//' @param p A probability, at which the GPD should be evaluated.
//' @param mu The location parameter of the GPD.
//' @param sigma The scale parameter of the GPD (should be strictly positive).
//' @param alpha The shape parameter of the GPD (valid for all real numbers).
//' @param n The number of random numbers which should be generated for the GPD.
//' @rdname GPD
//' @return dGPD: The estimated value of the density function of the GPD at the point \code{x}.
//' @examples
//' library(alR)
//' dGPD(0.5, 0, 1, 2)
//' @export
// [[Rcpp::export]]
NumericVector dGPD(NumericVector x, double mu, double sigma, double alpha)
{
if (sigma<=0.0) stop("The scale parameter of the GPD should be strictly positive.");
if (alpha<=0.0 && is_true(any(x<mu))) warning("For the chosen parameters, x>=%f.", mu);
if (alpha>0.0 && (is_true(any(x<mu)) || is_true(any(x>mu+(sigma/alpha))))) warning("For the chosen parameters, %f<=x<=%f.", mu, mu+(sigma/alpha));

if (alpha==0) return (1/sigma)*exp(-(x-mu)/sigma);
else return (1/sigma)*pow(1-alpha*(x-mu)/sigma, (1/alpha)-1);
}


//' @rdname GPD
//' @examples
//' pGPD(0.5, 0, 1, 2)
//' @return pGPD: The value of the cumulative distribution function of the GPD at the point \code{x}.
//' @export
// [[Rcpp::export]]
NumericVector pGPD(NumericVector q, double mu, double sigma, double alpha)
{
if (sigma<=0.0) stop("The scale parameter of the GPD should be strictly positive.");
if (alpha<=0.0 && is_true(any(q<mu))) warning("For the chosen parameters, x>=%f.", mu);
if (alpha>0.0 && (is_true(any(q<mu)) || is_true(any(q>mu+(sigma/alpha))))) warning("For the chosen parameters, %f<=x<=%f.", mu, mu+(sigma/alpha));

if (alpha==0) return 1-exp(-(q-mu)/sigma);
else return 1-pow(1-alpha*(q-mu)/sigma, 1/alpha);
}


//' @rdname GPD
//' @examples
//' qGPD(0.5, 0, 1, 2)
//' @return qGPD: The \code{x}th quantile of the GPD.
//' @export
// [[Rcpp::export]]
double qGPD(double p, double mu, double sigma, double alpha)
{
if (sigma<=0.0) stop("The scale parameter of the GPD should be strictly positive.");
if (p<0 || p>1) stop("Probabilities are only defined in the range [0, 1].");

if (alpha==0) return mu-sigma*std::log(1-p);
else return mu+(sigma/alpha)*(1-pow(1-p, alpha));
}


//' @rdname GPD
//' @examples
//' rGPD(10, 0, 1, 2)
//' @return rGPD: \code{n} random numbers from the GPD.
//' @export
// [[Rcpp::export]]
NumericVector rGPD(const int n, double mu, double sigma, double alpha)
{
if (sigma<=0.0) stop("The scale parameter of the GPD should be strictly positive.");

NumericVector x(n);
NumericVector u = runif(n);

for (int i=0; i<n; i++)
{
x[i] = qGPD(u[i], mu, sigma, alpha);
}
return x;
}