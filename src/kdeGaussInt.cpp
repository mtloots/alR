#include <Rcpp.h>
#include <R_ext/Applic.h>
using namespace Rcpp;
#include "alR_ext.h"


typedef struct { NumericVector mu; double h; } Params;


List Rcpp_integrate(integr_fn f, void *ex, double lower, double upper, int subdiv = 100, double eps_abs = 1e-10, double eps_rel = 1e-10)
{
int lenw = 4 * subdiv;
int *iwork = new int[subdiv];
double *work = new double[lenw];
double value;
double abs_err;
int subdiv_used;
int neval;
int error;

Rdqags(f, ex, &lower, &upper, &eps_abs, &eps_rel, &value, &abs_err, &neval, &error, &subdiv, &lenw, &subdiv_used, iwork, work);

delete [] iwork;
delete [] work;

return List::create(_["value"] = value,
_["abs.err"] = abs_err,
_["subdivisions"] = subdiv_used,
_["neval"] = neval);
}


void apdfkdeGauss(double *t, int n, void *ex)
{
Params *param = (Params *) ex;
NumericVector mu = param->mu;
double h = param->h;

for(int i = 0; i < n; i++)
{
t[i] = sqrt(1+(1/(2*PI*pow(mu.size(), 2))*pow(h, 6))*sum(pow((t[i]-mu)*exp(-0.5*pow((t[i]-mu)/h, 2)), 2)));
}
}


//' Arc length of Gaussian KDE.
//'
//' Calculate the arc length for a univariate Gaussian kernel density estimator over a specified interval.
//'
//' The arc length of a univariate Gaussian kernel density estimator is approximated using the numerical integration C code implimented for R's integrate functions, i.e. using Rdqags.  For this approximation, subdiv = 100 (100 subdivisions), and eps_abs = eps_rel = 1e-10, i.e. the absolute and relative errors respectively.
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
p1 = qkdeGauss(q1, mu, h)[0];
p2 = qkdeGauss(q2, mu, h)[0];
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
result[i] = kdeGaussInt(mu, h, q1[i], q2[i], quantile)[0];
}

return result;
}