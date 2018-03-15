#include <Rcpp.h>
#include <R_ext/Applic.h>
using namespace Rcpp;
#include "alR.h"

typedef struct { double mu; double sigma; double alpha; } Params;


void apdfGPD(double *t, int n, void *ex)
{
Params *param = (Params *) ex;
double mu = param->mu;
double sigma = param->sigma;
double alpha = param->alpha;

if (alpha==0)
{
for(int i = 0; i < n; i++)
{
t[i] = pow(1+pow((1/sigma)*exp(-(t[i]-mu)/sigma), 2), 0.5);
}
}
else
{
for(int i = 0; i < n; i++)
{
t[i] = pow(1+pow(((1-alpha)/pow(sigma, 2))*pow(1-alpha*(t[i]-mu)/sigma, (1/alpha)-2), 2), 0.5);
}
}
}


//' Arc length of GPD PDF.
//'
//' Calculate the arc length for a univariate generalised Pareto probability density function over a specified interval.
//'
//' The arc length of a univariate generalised Pareto probability density function is approximated using the numerical integration C code implimented for R's integrate functions, i.e. using Rdqags.  For this approximation, subdiv = 100 (100 subdivisions), and eps_abs = eps_rel = 1e-10, i.e. the absolute and relative errors respectively.
//'
//' @param mu A real number specifying the location parameter.
//' @param sigma A strictly positive real number specifying the scale parameter.
//' @param alpha A real number specifying the shape parameter.
//' @param q1 The point (or vector for \code{GPDInt2}) specifying the lower limit of the arc length integral.
//' @param q2 The point (or vector for \code{GPDInt2}) specifying the upper limit of the arc length integral.
//' @param quantile Logical, TRUE/FALSE, whether \code{q1} and \code{q2} are quantiles, or actual points in the domain.
//'
//' @return GPDInt: A list with the following components:
//' \itemize{
//' \item value: The resultant arc length.
//' \item abs.err: The absolute error between iterations.
//' subdivisions: Number of subdivisions used in the numerical approximation.
//' \item neval: Number of function evaluations used by the numerical approximation.
//' }
//'
//' @examples
//' library(alR)
//' mu <- 0
//' sigma <- 1
//' alpha <- 2
//' GPDInt(mu, sigma, alpha, 0.025, 0.975, TRUE)
//' GPDInt(mu, sigma, alpha, 0.001, 0.5, FALSE)
//'
//' @export
// [[Rcpp::export]]
List GPDInt(double mu, double sigma, double alpha, double q1, double q2, bool quantile)
{
double p1, p2;
Params param = {mu, sigma, alpha};

if (quantile)
{
p1 = qGPD(q1, mu, sigma, alpha);
p2 = qGPD(q2, mu, sigma, alpha);
}
else
{
if (sigma<=0.0) stop("The scale parameter of the GPD should be strictly positive.");
if (alpha<=0.0 && (q1<mu || q2<mu)) warning("For the chosen parameters, x>=%f.", mu);
if (alpha>0.0 && (q1<mu || q2<mu || q1>mu+(sigma/alpha) || q2>mu+(sigma/alpha))) warning("For the chosen parameters, %f<=x<=%f.", mu, mu+(sigma/alpha));

p1 = q1;
p2 = q2;
}

return Rcpp_integrate(apdfGPD, &param, p1, p2);
}


//' @rdname GPDInt
//' @return GPDInt2: A vector having length equal to that of the vector of lower quantile bounds, containing the arc lengths requested for a generalised Pareto probability density function.
//'
//' @examples
//' GPDInt2(mu, sigma, alpha, c(0.025, 0.5), c(0.5, 0.975), TRUE)
//' GPDInt2(mu, sigma, alpha, c(-1.96, 0), c(0, 1.96), FALSE)
//'
//' @export
// [[Rcpp::export]]
NumericVector GPDInt2(double mu, double sigma, double alpha, NumericVector q1, NumericVector q2, bool quantile)
{
int i, n = q1.size();
NumericVector result(n);

for (i=0; i<n; i++)
{
result[i] = as<double>(GPDInt(mu, sigma, alpha, q1[i], q2[i], quantile)["value"]);
}

return result;
}