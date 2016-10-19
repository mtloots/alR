#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
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

//' @rdname alEfit
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


//' @rdname alEfit
//' @return alEfitdist: A matrix of parameter estimates resulting from the estimated arc lengths over the specified interval(s), i.e. the bootstrap distribution for the estimated parameters resulting from the chosen sample arc length statistic.
//' @examples
//' \dontrun{
//' alEfitdist(x, 0.025, 0.975, TRUE, -1, 100)
//' alEfitdist(x, 0.025, 0.975, FALSE, -1, 100)
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix alEfitdist(NumericVector x, NumericVector q1, NumericVector q2, bool dc, double type, int bootstraps)
{
NumericMatrix sampDist(bootstraps, 2);
NumericVector prob = NumericVector::create();
int n = x.size();

for (int i=0; i<bootstraps; i++)
{
NumericVector sample = RcppArmadillo::sample(x, n, true, prob);
sampDist(i, _) = as<NumericVector>(wrap(alE(sample, q1, q2, dc, type)[0]));
}

return sampDist;
}


//' @rdname alEtest
//' @return alEdist: A vector (matrix) of arc lengths over the specified interval(s), i.e. the simulated distribution for the chosen sample arc length statistic.
//' @examples
//' \dontrun{
//' alEdist(50, 100, 2, 3.5, 0.025, 0.975, TRUE, TRUE, -1)
//' alEdist(50, 100, 2, 3.5, c(0.025,0.5), c(0.5,0.975), TRUE, TRUE, -1)
//' alEdist(50, 100, 2, 3.5, 0.025, 0.975, TRUE, FALSE, -1)
//' alEdist(50, 100, 2, 3.5, c(0.025,0.5), c(0.5,0.975), TRUE, FALSE, -1)
//' alEdist(50, 100, 2, 3.5, qnorm(0.025,2,3.5), qnorm(0.975, 2, 3.5), FALSE, FALSE, -1)
//' alEdist(50, 100, 2, 3.5, c(qnorm(0.025, 2, 3.5),2), c(2,qnorm(0.975, 2, 3.5)), FALSE, FALSE, -1)
//' }
//' @export
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