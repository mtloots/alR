#include <Rcpp.h>
using namespace Rcpp;
#include <R_ext/Applic.h>
#include "alR.h"


typedef struct { NumericVector mu; double h; } Params;


void apdfkdeGauss(double *t, int n, void *ex)
{
Params *param = (Params *) ex;
NumericVector mu = param->mu;
double h = param->h;

for(int i = 0; i < n; i++)
{
t[i] = sqrt(1+(1/(2*PI*pow(mu.size(), 2)*pow(h, 6)))*pow(sum((t[i]-mu)*exp(-0.5*pow((t[i]-mu)/h, 2))), 2));
}
}


//' Arc length of Gaussian KDE.
//'
//' Calculate the arc length for a univariate Gaussian kernel density estimator over a specified interval.
//'
//' For \code{kdeGaussInt} and \code{kdeGaussInt2}, the arc length of a univariate Gaussian kernel density estimator is approximated using the numerical integration C code implimented for R's integrate functions, i.e. using Rdqags.  For this approximation, subdiv = 100 (100 subdivisions), and eps_abs = eps_rel = 1e-10, i.e. the absolute and relative errors respectively.
//'
//' For \code{kdeGaussIntApprox}, the arc length is approximated by constructing the KDE, and then calculated as the sum of a finite collection of straight lines, based on the Pythagorean theorem.
//'
//' @param mu A vector of data points on which the kernel density estimator is based.
//' @param h The kernel density estimator bandwidth.
//' @param q1 The point (or vector for \code{kdeGaussInt2}) specifying the lower limit of the arc length integral.
//' @param q2 The point (or vector for \code{kdeGaussInt2}) specifying the upper limit of the arc length integral.
//' @param quantile Logical, TRUE/FALSE, whether \code{q1} and \code{q2} are quantiles, or actual points in the domain.
//'
//' @return kdeGaussInt: A list with the following components:
//' \itemize{
//' \item value: The resultant arc length.
//' \item abs.err: The absolute error between iterations.
//' subdivisions: Number of subdivisions used in the numerical approximation.
//' \item neval: Number of function evaluations used by the numerical approximation.
//' }
//'
//' @examples
//' library(alR)
//' mu <- rnorm(100)
//' h <- bw(mu, type=1)
//' kdeGaussInt(mu, h, 0.025, 0.975, TRUE)
//' kdeGaussInt(mu, h, -1.96, 1.96, FALSE)
//'
//' @export
// [[Rcpp::export]]
List kdeGaussInt(NumericVector mu, double h, double q1, double q2, bool quantile)
{
double p1, p2;
Params param = {mu, h};

if (quantile)
{
p1 = as<double>(qkdeGauss(q1, mu, h)["result"]);
p2 = as<double>(qkdeGauss(q2, mu, h)["result"]);
}
else
{
p1 = q1;
p2 = q2;
}

return Rcpp_integrate(apdfkdeGauss, &param, p1, p2);
}


//' @rdname kdeGaussInt
//' @return kdeGaussInt2: A vector having length equal to that of the vector of lower quantile bounds, containing the arc lengths requested for a Gaussian kernel density estimator.
//'
//' @examples
//' kdeGaussInt2(mu, h, c(0.025, 0.5), c(0.5, 0.975), TRUE)
//' kdeGaussInt2(mu, h, c(-1.96, 0), c(0, 1.96), FALSE)
//'
//' @export
// [[Rcpp::export]]
NumericVector kdeGaussInt2(NumericVector mu, double h, NumericVector q1, NumericVector q2, bool quantile)
{
int i, n = q1.size();
NumericVector result(n);

for (i=0; i<n; i++)
{
result[i] = as<double>(kdeGaussInt(mu, h, q1[i], q2[i], quantile)["value"]);
}

return result;
}


//' @rdname kdeGaussInt
//' @return kdeGaussIntApprox: The resultant arc length.
//' @examples
//' kdeGaussIntApprox(mu, h, 0.025, 0.975, TRUE)
//' kdeGaussIntApprox(mu, h, -1.96, 1.96, FALSE)
//'
//' @export
// [[Rcpp::export]]
double kdeGaussIntApprox(NumericVector mu, double h, double q1, double q2, bool quantile)
{
double Q1, Q2;

if (quantile)
{
Q1 = qsamp(mu, q1);
Q2 = qsamp(mu, q2);
}
else
{
Q1 = q1;
Q2 = q2;
}

NumericVector x = mu[(mu >= Q1) & (mu <= Q2)];
int n = x.size();
std::sort(x.begin(), x.end());
NumericVector xHat(n);

for (int i=0; i<n; i++)
{
xHat[i] = dkdeGauss(x[i], x, h);
}

return sum(sqrt(pow(diff(x), 2)+pow(diff(xHat), 2)));
}

//' @rdname kdeGaussInt
//' @return kdeGaussIntApprox2: A vector having length equal to that of the vector of lower quantile bounds, containing the discrete arc lengths requested for a Gaussian kernel density estimator.
//'
//' @examples
//' kdeGaussIntApprox2(mu, h, c(0.025, 0.5), c(0.5, 0.975), TRUE)
//' kdeGaussIntApprox2(mu, h, c(-1.96, 0), c(0, 1.96), FALSE)
//'
//' @export
// [[Rcpp::export]]
NumericVector kdeGaussIntApprox2(NumericVector mu, double h, NumericVector q1, NumericVector q2, bool quantile)
{
int i, n = q1.size();
NumericVector result(n);

for (i=0; i<n; i++)
{
result[i] = kdeGaussIntApprox(mu, h, q1[i], q2[i], quantile);
}

return result;
}