#include <Rcpp.h>
#include <R_ext/Applic.h>
using namespace Rcpp;
#include "alR.h"


typedef struct { NumericVector mu; double h; } Params;
typedef struct { NumericVector mu; double h; double q; } Qarams;


//' Gaussian kernel density estimator.
//'
//' Estimate a density function using a kernel density estimator with a Gaussian kernel.
//'
//' The cumulative distribution function is calculated using the numerical integration C code implimented for R's integrate functions, i.e. using Rdqagi.  For this approximation, subdiv = 100 (100 subdivisions), and eps_abs = eps_rel = 1e-10, i.e. the absolute and relative errors respectively.
//'
//' The quantiles of the Gaussian kernel density estimator are calculated using Brent's method.  This method requires an interval in which a solution is saught.  The objective funcion for which a zero is saught is \code{\link{dkdeGauss}}-\code{x}, where \code{x} is the quantile saught.  The first interval in which a solution is searched for, corresponds to the range of \code{mu}, and is expanded in multiples thereof in consequtive steps.  The maximum number of iterations is set at 1000, and the accuracy saught between iterations, is set at 1e-10.
//'
//' @param x A data point, or quantile, at which the kernel density estimator should be evaluated.
//' @param mu A vector of data points on which the kernel density estimator is based.
//' @param h The kernel density estimator bandwidth.
//' @rdname kdeGauss
//' @return dkdeGauss: The estimated value of the density function at the point x.
//' @examples
//' library(alR)
//' x <- rnorm(100)
//' h_x <- bw(x, type=1)
//' dkdeGauss(0, x, h_x)
//' @export
// [[Rcpp::export]]
double dkdeGauss(double x, NumericVector mu, double h)
{
return (1.0/(mu.size()*pow(2*PI, 0.5)*h))*sum(exp(-0.5*pow((x-mu)/h, 2)));
}

void pdfkdeGauss(double *t, int n, void *ex)
{
Params *param = (Params *) ex;
NumericVector mu = param->mu;
double h = param->h;

for(int i = 0; i < n; i++)
{
t[i] = dkdeGauss(t[i], mu, h);
}
} 


//' @rdname kdeGauss
//' @examples
//' pkdeGauss(0, x, h_x)
//' @return pkdeGauss: The estimated value of the cumulative distribution function at the point \code{x}.
//' @export
// [[Rcpp::export]]
double pkdeGauss(double x, NumericVector mu, double h)
{
int n = mu.size();
double kde=0;

for (int i=0; i<n; i++)
{
kde += R::pnorm(x, mu[i], h, true, false);
}

return kde/n;
}

double qkde(double x, void *ex)
{
Qarams *param = (Qarams *) ex;
NumericVector mu = param->mu;
double h = param->h;
double q = param->q;

return double(pkdeGauss(x, mu, h))-q;
}


//' @rdname kdeGauss
//' @examples
//' qkdeGauss(0.5, x, h_x)
//' @return qkdeGauss: A list with the following components:
//' \itemize{
//' \item result: The \code{x}th quantile of the Gaussian kernel density estimator.
//' \item value: The value of the cumulative distribution function of the Gaussian kernel density estimator at the \code{x}th quantile.
//' \item obj.fun: The value of the objective function resulting from Brent's method; should be less than 1e-10.
//' \item iterations: Number of iterations for Brent's method in order to achieve the desired accuracy.
//' \item steps: Number of range expansions of the search boundaries for Brent's method.
//' }
//' @export
// [[Rcpp::export]]
List qkdeGauss(double x, NumericVector mu, double h)
{
Qarams param = {mu, h, x};
double minmu = min(mu);
double maxmu = max(mu);
double rangemu = maxmu-minmu;
int steps = 0;
double lower, upper;
List bf;

for (int i=1; i<=1000; i++)
{
steps++;
lower = steps*minmu-(steps-1)*rangemu/2;
upper = steps*maxmu+(steps-1)*rangemu/2;
bf = brents_fun(qkde, &param, lower, upper, 1e-10, 1000);
if (as<int>(bf["msg"]) == 0)
{
i = 1001;
}
}

return List::create(_["result"] = as<double>(bf["result"]),
_["value"] = pkdeGauss(as<double>(bf["result"]), mu, h),
_["obj.fun"] = as<double>(bf["value"]),
_["iterations"] = as<int>(bf["iterations"]),
_["steps"] = steps);
}