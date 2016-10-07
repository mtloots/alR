#include <Rcpp.h>
#include <R_ext/Applic.h>
using namespace Rcpp;
#include "alR.h"

typedef struct { double mu; double sigma; } Params;


void apdfGauss(double *t, int n, void *ex)
{
Params *param = (Params *) ex;
double mu = param->mu;
double sigma = param->sigma;

for(int i = 0; i < n; i++)
{
t[i] = sqrt(1+(1/(2*PI*pow(sigma, 6)))*pow(t[i]-mu, 2)*exp(-pow((t[i]-mu)/sigma, 2)));
}
}


//' Arc length of Gaussian PDF.
//'
//' Calculate the arc length for a univariate Gaussian probability density function over a specified interval.
//'
//' The arc length of a univariate Gaussian probability density function is approximated using the numerical integration C code implimented for R's integrate functions, i.e. using Rdqags.  For this approximation, subdiv = 100 (100 subdivisions), and eps_abs = eps_rel = 1e-10, i.e. the absolute and relative errors respectively.
//'
//' @param mu A real number specifying the location parameter.
//' @param sigma A positive real number specifying the scale parameter.
//' @param q1 The point (or vector for \code{GaussInt2}) specifying the lower limit of the arc length integral.
//' @param q2 The point (or vector for \code{GaussInt2}) specifying the upper limit of the arc length integral.
//' @param quantile Logical, TRUE/FALSE, whether \code{q1} and \code{q2} are quantiles, or actual points in the domain.
//'
//' @return GaussInt: A list with the following components:
//' \itemize{
//' \item value: The resultant arc length.
//' \item abs.err: The absolute error between iterations.
//' subdivisions: Number of subdivisions used in the numerical approximation.
//' \item neval: Number of function evaluations used by the numerical approximation.
//' }
//'
//' @examples
//' library(alR)
//' mu <- 2
//' sigma <- 3.5
//' GaussInt(mu, sigma, 0.025, 0.975, TRUE)
//' GaussInt(mu, sigma, -1.96, 1.96, FALSE)
//'
//' @export
// [[Rcpp::export]]
List GaussInt(double mu, double sigma, double q1, double q2, bool quantile)
{
double p1, p2;
Params param = {mu, sigma};

if (quantile)
{
p1 = R::qnorm(q1, mu, sigma, 1, 0);
p2 = R::qnorm(q2, mu, sigma, 1, 0);
}
else
{
p1 = q1;
p2 = q2;
}

return Rcpp_integrate(apdfGauss, &param, p1, p2);
}


//' @rdname GaussInt
//' @return GaussInt2: A vector having length equal to that of the vector of lower quantile bounds, containing the arc lengths requested for a Gaussian probability density function.
//'
//' @examples
//' GaussInt2(mu, sigma, c(0.025, 0.5), c(0.5, 0.975), TRUE)
//' GaussInt2(mu, sigma, c(-1.96, 0), c(0, 1.96), FALSE)
//'
//' @export
// [[Rcpp::export]]
NumericVector GaussInt2(double mu, double sigma, NumericVector q1, NumericVector q2, bool quantile)
{
int i, n = q1.size();
NumericVector result(n);

for (i=0; i<n; i++)
{
result[i] = GaussInt(mu, sigma, q1[i], q2[i], quantile)[0];
}

return result;
}