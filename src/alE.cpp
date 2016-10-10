#include <Rcpp.h>
#include <R_ext/Applic.h>
using namespace Rcpp;
#include "alR.h"


typedef struct { NumericVector sampAL; NumericVector q1; NumericVector q2; } Aux;


double alEobj(int n, double *par, void *ex)
{
Aux *aux = (Aux *) ex;
NumericVector sampAL = aux->sampAL;
NumericVector q1 = aux->q1;
NumericVector q2 = aux->q2;

NumericVector theoAL = GaussInt2(par[0], par[1], q1, q2, false);

return sqrt(sum(pow(theoAL-sampAL, 2)));
}


//' Arc length estimation.
//'
//' A framework for arc length estimation.
//'
//'
//' \itemize{
//' \item Estimate distributional parameters using the method of arc lengths.
//' \item Simulate distributions for sample arc length statistics.
//' }
//'
//' This method is currently only implimented for the normal distribution.  The underlying C code for the Nelder-Mead method of the optim function is used for optimising the objective function.  The tolarence level is set at 1e-15, and a maximum number of 1000 iterations is allowed.  The maximum likelihood estimates are used as initial values for the Nelder-Mead algorithm.
//'
//' @param x A vector of sample values.
//' @param mu A real value specifying the mean of the normal distribution.
//' @param sigma A positive real number specifying the scale parameter of the normal distribution.
//' @param q1,q2 Vectors specifying the quantiles (or points if quantile=FALSE) over which arc length segments are to be computed.
//' @param quantile TRUE/FALSE whether q1 and q2 are quantiles, or elements of the domain of \code{x}.
//' @param dc TRUE/FALSE:  Should the discrete or continuous sample statistic be used.
//' @param type The type of bandwidth estimator for the underlying KDE; see \code{\link{bw}}.
//' @param n An integer specifying the sample size.
//' @param bootstraps An integer specifying the size of the parametric bootstrap.
//'
//' @return alE: A list with the following components (see \code{\link{optim}}):
//' \itemize{
//' \item par: The estimated parameters.
//' \item abstol: The absolute tolarence level (default 1e-15).
//' \item fail: An integer code indicating convergence.
//' \item fncount: Number of function evaluations.
//' }
//'
//' @examples
//' x <- rnorm(1000)
//' alE(x,0.025, 0.975, TRUE, -1)
//' alE(x,c(0.025, 0.5), c(0.5, 0.975), TRUE, -1)
//' alE(x,0.025, 0.975, FALSE, -1)
//' alE(x,c(0.025, 0.5), c(0.5, 0.975), FALSE, -1)
//'
//' @export
// [[Rcpp::export]]
List alE(NumericVector x, NumericVector q1, NumericVector q2, bool dc, double type)
{
double h = bw(x, type);
int n = q1.size();
NumericVector p1(n), p2(n), sampAL(n);

for (int i=0; i<n; i++)
{
p1[i] = qsamp(x, q1[i]);
p2[i] = qsamp(x, q2[i]);
}

if (dc)
{
sampAL = kdeGaussIntApprox2(x, h, p1, p2, false);
}
else
{
 sampAL = kdeGaussInt2(x, h, p1, p2, false);
}

Aux aux = {sampAL, p1, p2};
double xin[2] = {mean(x), sd(x)};

return Rcpp_nmmin(2, alEobj, xin, &aux);
}


//' @rdname alE
//' @return alEdist: A vector (matrix) of arc lengths over the specified interval(s), i.e. the simulated distribution for the chosen sample arc length statistic.
// [[Rcpp::export]]
NumericMatrix alEdist(int n, int bootstraps, double mu, double sigma, NumericVector q1, NumericVector q2, bool quantile, bool dc, double type)
{
NumericMatrix sampDist(bootstraps, q1.size());

for (int i=0; i<bootstraps; i++)
{
NumericVector x = rnorm(n, mu, sigma);
double h = bw(x, type);

if (dc)
{
sampDist(i, _) = kdeGaussIntApprox2(x, h, q1, q2, quantile);
}
else
{
 sampDist(i, _) = kdeGaussInt2(x, h, q1, q2, quantile);
}
}

return sampDist;
}