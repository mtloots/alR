#include <Rcpp.h>
#include <R_ext/Applic.h>
using namespace Rcpp;
#include "alR.h"

typedef struct { double mu; double sigma; double h; double k; double tau; } Params;


void acdfKappa4(double *t, int n, void *ex)
{
Params *param = (Params *) ex;
double mu = param->mu;
double sigma = param->sigma;
double h = param->h;
double k = param->k;
double tau = param->tau;

for(int i = 0; i < n; i++)
{
t[i] = pow(1+pow(dkappa4(t[i], mu, sigma, h, k)/tau, 2), 0.5);
}
}


//' Arc length of four-parameter kappa CDF.
//'
//' Calculate the arc length for a univariate four-parameter kappa cumulative distribution function over a specified interval.
//'
//' The arc length of a univariate four-parameter kappa cumulative distribution function is approximated using the numerical integration C code implimented for R's integrate functions, i.e. using Rdqags.  For this approximation, subdiv = 100 (100 subdivisions), and eps_abs = eps_rel = 1e-10, i.e. the absolute and relative errors respectively.
//'
//' @param mu A real number specifying the location parameter.
//' @param sigma A positive real number specifying the scale parameter.
//' @param h,k Real numbers specifying the two shape parameters.
//' @param tau A real number between 0 and 1, corresponding to the CDF value at the point of truncation.
//' @param q1 The point (or vector for \code{kappa4Int2}) specifying the lower limit of the arc length integral.
//' @param q2 The point (or vector for \code{kappa4Int2}) specifying the upper limit of the arc length integral.
//' @param quantile Logical, TRUE/FALSE, whether \code{q1} and \code{q2} are quantiles, or actual points in the domain.
//'
//' @return kappa4Int: A list with the following components:
//' \itemize{
//' \item value: The resultant arc length.
//' \item abs.err: The absolute error between iterations.
//' subdivisions: Number of subdivisions used in the numerical approximation.
//' \item neval: Number of function evaluations used by the numerical approximation.
//' }
//'
//' @examples
//' library(alR)
//' mu <- 4
//' sigma <- 0.4
//' h <- -4
//' tau <- 1 ## no truncation
//' k <- kappa4tc(-4, 0, 1)$par
//' kappa4Int(mu, sigma, h, k, tau, 0.025, 0.975, TRUE)
//' p1 <- qkappa4(0.025, mu, sigma, h, k)
//' p2 <- qkappa4(0.975, mu, sigma, h, k)
//' kappa4Int(mu, sigma, h, k, tau, p1, p2, FALSE)
//'
//' @export
// [[Rcpp::export]]
List kappa4Int(double mu, double sigma, double h, double k, double tau, double q1, double q2, bool quantile)
{
double p1, p2;
Params param = {mu, sigma, h, k, tau};

if (quantile)
{
p1 = qkappa4(q1, mu, sigma, h, k);
p2 = qkappa4(q2, mu, sigma, h, k);
}
else
{
p1 = q1;
p2 = q2;
}

return Rcpp_integrate(acdfKappa4, &param, p1, p2);
}


//' @rdname kappa4Int
//' @return kappa4Int2: A vector having length equal to that of the vector of lower quantile bounds, containing the arc lengths requested for a four-parameter kappa cumulative distribution function.
//'
//' @examples
//' kappa4Int2(mu, sigma, h, k, tau, c(0.025, 0.5), c(0.5, 0.975), TRUE)
//' p12 <- qkappa4(0.5, mu, sigma, h, k)
//' kappa4Int2(mu, sigma, h, k, tau, c(p1, p12), c(p12, p2), FALSE)
//'
//' @export
// [[Rcpp::export]]
NumericVector kappa4Int2(double mu, double sigma, double h, double k, double tau, NumericVector q1, NumericVector q2, bool quantile)
{
int n = q1.size();
NumericVector result(n);

for (int i=0; i<n; i++)
{
result[i] = as<double>(kappa4Int(mu, sigma, h, k, tau, q1[i], q2[i], quantile)["value"]);
}

return result;
}


//' Sample arc length statistic.
//'
//' The arc length over a specified interval is calculated for use in non-linear estimation.
//'
//' @param x Numeric vector of independent outcomes.
//' @param y Numeric vector of dependent outcomes \eqn{y=F(x)}.
//' @param q1,q2 Quantiles (between 0 and 1) over which the arc length segment is to be computed.
//' @return kappa4IntApprox: The resultant arc length.
//' @param quantile Logical, TRUE/FALSE, whether \code{q1} and \code{q2} are quantiles, or actual points in the domain.
//' @examples
//' x <- rnorm(100)
//' y <- pnorm(x)
//' kappa4IntApprox(x, y, 0.025, 0.975, TRUE)
//' kappa4IntApprox(x, y, -1.96, 1.96, FALSE)
//'
//' @export
// [[Rcpp::export]]
double kappa4IntApprox(NumericVector x, NumericVector y, double q1, double q2, bool quantile)
{
int Q1, Q2;
LogicalVector select;

if (quantile)
{
Q1 = sum(ifelse(y < q1, 1, 0));
Q2 = sum(ifelse(y <= q2, 1, 0));
select = ifelse((x >= x[Q1]) & (x <= x[Q2]), true, false);
}
else
{
Q1 = q1;
Q2 = q2;
select = ifelse((x >= Q1) & (x <= Q2), true, false);
}

NumericVector x_seg = x[select];
NumericVector y_seg = y[select];

return sum(pow(pow(diff(x_seg), 2)+pow(diff(y_seg), 2), 0.5));
}


//' @rdname kappa4IntApprox
//' @return kappa4IntApprox2: A vector having length equal to that of the vector of lower quantile bounds, containing the discrete arc length segments over the specified intervals.
//'
//' @examples
//' kappa4IntApprox2(x, y, c(0.025, 0.5), c(0.5, 0.975), TRUE)
//' kappa4IntApprox2(x, y, c(-1.96, 0), c(0, 1.96), FALSE)
//'
//' @export
// [[Rcpp::export]]
NumericVector kappa4IntApprox2(NumericVector x, NumericVector y, NumericVector q1, NumericVector q2, bool quantile)
{
int i, n = q1.size();
NumericVector result(n);

for (i=0; i<n; i++)
{
result[i] = kappa4IntApprox(x, y, q1[i], q2[i], quantile);
}

return result;
}