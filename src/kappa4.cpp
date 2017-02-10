#include <Rcpp.h>
using namespace Rcpp;
#include "alR.h"


typedef struct { double mu; double sigma; double h; } params;

typedef struct { NumericVector x; NumericVector F; double x_min; double x_max; int N;} Sparams;

double gkappa4(double x, double mu, double sigma, double k)
{
return 1-k*(x-mu)/sigma;
}


//' Four-parameter kappa distribution.
//'
//' Functions for the four-parameter kappa distribution.
//'
//' @param x A data point, or quantile, at which the four-parameter kappa distribution should be evaluated.
//' @param mu A real value representing the location of the distribution.
//' @param sigma A positive real number representing the scale parameter of the distribution.
//' @param h,k Real numbers representing shape parameters of the distribution.
//' @param n Number of random variates to generate.
//' @rdname kappa4
//' @examples
//' pkappa4(1, 1, 2, 0.5, 2)
//' @return pkappa4: The cumulative distribution function at the point \code{x}.
//' @export
// [[Rcpp::export]]
double pkappa4(double x, double mu, double sigma, double h, double k)
{
return pow(1-h*pow(gkappa4(x, mu, sigma, k), 1.0/k), 1.0/h);
}


//' @rdname kappa4
//' @return dkappa4: The density function at the point x.
//' @examples
//' dkappa4(1, 1, 2, 0.5, 2)
//' @export
// [[Rcpp::export]]
double dkappa4(double x, double mu, double sigma, double h, double k)
{
return (1.0/sigma)*pow(pkappa4(x, mu, sigma, h, k), 1-h)*pow(gkappa4(x, mu, sigma, k), (1.0/k)-1);
}


//' @rdname kappa4
//' @examples
//' qkappa4(0.25, 1, 2, 0.5, 2)
//' @return qkappa4: The \code{x}th quantile of the distribution.
//' @export
// [[Rcpp::export]]
double qkappa4(double x, double mu, double sigma, double h, double k)
{
if (h == 0)
{
return mu+(sigma/k)*(1-pow(log(1/x), k));
}
else
{
return mu+(sigma/k)*(1-pow((1-pow(x, h))/h, k));
}
}


//' @rdname kappa4
//' @examples
//' rkappa4(10, 1, 2, 0.5, 2)
//' @return rkappa4: Randomly generated numbers from the distribution.
//' @export
// [[Rcpp::export]]
NumericVector rkappa4(int n, double mu, double sigma, double h, double k)
{
NumericVector u = runif(n);
NumericVector k4(n);

for (int i=0; i<n; i++)
{
k4[i] = qkappa4(u[i], mu, sigma, h, k);
}

return k4;
}


double ddkappa4(double x, double mu, double sigma, double h, double k)
{
return (1/sigma)*dkappa4(x, mu, sigma, h, k)*pow(gkappa4(x, mu, sigma, k), -1)*((1-h)*pow(pow(gkappa4(x, mu, sigma, k), -1.0/k)-h, -1.0)-(1-k));
}


//' @rdname kappa4
//' @return dddkappa4: The second derivative of dkappa4.
//' @examples
//' dddkappa4(1, 1, 2, 0.5, 2)
//' @export
// [[Rcpp::export]]
double dddkappa4(double x, double mu, double sigma, double h, double k)
{
return pow(gkappa4(x, mu, sigma, k), -1.0)*(((1-h)/sigma)*pow(pkappa4(x, mu, sigma, h, k), -h)*pow(gkappa4(x, mu, sigma, k), 1.0/k)*(2*ddkappa4(x, mu, sigma, h, k)-(1.0/pow(sigma, 2))*pow(pkappa4(x, mu, sigma, h, k), 1-2*h)*pow(gkappa4(x, mu, sigma, k), 2.0/k-2))-((1-k)/sigma)*(ddkappa4(x, mu, sigma, h, k)+(k/pow(sigma, 2))*pow(pkappa4(x, mu, sigma, h, k), 1-h)*pow(gkappa4(x, mu, sigma, k), 1.0/k-2)));
}


//' @rdname kappa4
//' @return kappa4cond: The resultant induction period (IP).
//' @examples
//' kappa4cond(1, 2, 0.5, 2)
//' @export
// [[Rcpp::export]]
double kappa4cond(double mu, double sigma, double h, double k)
{
return mu+(sigma/k)*(1-pow((1-k)/(1-h*k), k-1));
}


double findK(int n, double *k, void *ex)
{
params *param = (params *) ex;
double mu = param->mu;
double sigma = param->sigma;
double h = param->h;
double cond = 0;

// Condition 4: Avoid -infty at lower bound
if ((k[0] >= 0) & (h <= 0))
{
cond += 1;
}

double x = kappa4cond(mu, sigma, h, k[0]);

return std::abs(dddkappa4(x, mu, sigma, h, k[0]))+cond;
}


//' @rdname kappa4
//' @return kappa4tc: A list with the following components:
//' \itemize{
//' \item $par: The k shape parameter corresponding to a given h parameter for the time-conductivity problem.
//' \item $abstol: The absolute tolerance for the numerical optimisation.
//' \item $fail: A code relating to the optimisation routine.
//' \item $fncount: Number of function calls.
//' }
//' @examples
//' kappa4tc(-4, 0, 1)
//' @export
// [[Rcpp::export]]
List kappa4tc(double h, double mu, double sigma)
{
params param = {mu, sigma, h};
double xin[1];

if (h < -1)
{
xin[0] = 1/(2*h);
}
else if (h < -0.89)
{
xin[0] = -0.34;
}
else
{
xin[0] = 1;
}

return Rcpp_nmmin(1, findK, xin, &param, 0.0, 1e-25, 1e-25);
}