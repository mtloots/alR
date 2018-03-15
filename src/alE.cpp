#include <Rcpp.h>
using namespace Rcpp;
#include "alR.h"


typedef struct { NumericVector sampAL; NumericVector q1; NumericVector q2; int distribution; } Aux;


double alEobj(int n, double *par, void *ex)
{
Aux *aux = (Aux *) ex;
NumericVector sampAL = aux->sampAL;
NumericVector q1 = aux->q1;
NumericVector q2 = aux->q2;
int distribution = aux->distribution;
NumericVector theoAL(q1.size());

switch (distribution)
{
case 1:
theoAL = GaussInt2(par[0], par[1], q1, q2, false);
break;
case 2:
theoAL = GPDInt2(par[0], std::abs(par[1]), par[2], q1, q2, false);
break;
default: throw exception("Unknown distribution selected.");
}

return pow(sum(pow(theoAL-sampAL, 2)), 0.5);
}


//' @rdname alEfit
//' @return alE: A list with the following components (see \code{\link{optim}}):
//' \itemize{
//' \item par: The estimated parameters.
//' \item abstol: The absolute tolerance level (default 1e-15).
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
List alE(NumericVector x, NumericVector q1, NumericVector q2, bool dc, double type, int distribution=1)
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

Aux aux = {sampAL, p1, p2, distribution};

switch (distribution)
{
case 1:
{
double xina[2] = {mean(x), sd(x)};
return Rcpp_nmmin(2, alEobj, xina, &aux);
}
case 2:
{
NumericVector temp = GPDqML(x);
double xinb[3] = {temp[0], temp[1], temp[2]};
return Rcpp_nmmin(3, alEobj, xinb, &aux);
}
default: throw exception("Unknown distribution selected.");
}
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
NumericMatrix alEfitdist(NumericVector x, NumericVector q1, NumericVector q2, bool dc, double type, int bootstraps, int distribution=1)
{
switch (distribution)
{
case 1:
{
NumericMatrix sampDist(bootstraps, 2);
int n = x.size();

for (int i=0; i<bootstraps; i++)
{
NumericVector sample = RcppSample(x, n);
sampDist(i, _) = as<NumericVector>(wrap(alE(sample, q1, q2, dc, type, distribution)["par"]));
}

return sampDist;
}
case 2:
{
NumericMatrix sampDist(bootstraps, 3);
int n = x.size();

for (int i=0; i<bootstraps; i++)
{
NumericVector sample = RcppSample(x, n);
sampDist(i, _) = as<NumericVector>(wrap(alE(sample, q1, q2, dc, type, distribution)["par"]));
}

return sampDist;
}
default:  throw exception("Unknown distribution selected.");
}
}


//' @rdname alEtest
//' @return alEdist: A vector (matrix) of arc lengths over the specified interval(s), i.e. the simulated distribution for the chosen sample arc length statistic.
//' @examples
//' \dontrun{
//' alEdist(50, 100, 2, 3.5, 0.025, 0.975, TRUE, TRUE, -1)
//' alEdist(50, 100, 2, 3.5, c(0.025,0.5), c(0.5,0.975), TRUE, TRUE, -1)
//' alEdist(50, 100, 2, 3.5, 0.025, 0.975, TRUE, FALSE, -1)
//' alEdist(50, 100, 2, 3.5, c(0.025,0.5), c(0.5,0.975), TRUE, FALSE, -1)
//' alEdist(50, 100, 2, 3.5, qnorm(0.025,2,3.5),
//' qnorm(0.975, 2, 3.5), FALSE, FALSE, -1)
//' alEdist(50, 100, 2, 3.5, c(qnorm(0.025, 2, 3.5),2),
//' c(2,qnorm(0.975, 2, 3.5)), FALSE, FALSE, -1)
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix alEdist(int n, int bootstraps, NumericVector mu, double sigma, NumericVector q1, NumericVector q2, bool quantile, bool dc, double type, int distribution=1)
{
NumericMatrix sampDist(bootstraps, q1.size());
NumericVector x(n);

for (int i=0; i<bootstraps; i++)
{
switch (distribution)
{
case 1:
x = rnorm(n, mu[0], sigma);
break;
case 2:
x = rGPD(n, mu[0], mu[1], mu[2]);
break;
default: throw exception("Unknown distribution selected.");
}

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