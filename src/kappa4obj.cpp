#include <Rcpp.h>
using namespace Rcpp;
#include "alR.h"


//' Sigmoidal curve fitting.
//'
//' Support functions for fitting four-parameter kappa sigmoidal curves.
//'
//' @param parms A numeric vector of parameters to be estimated.
//' @param  xvec A numeric vector of independent observations.
//' @param y A numeric vector of dependent observations.
//' @param x_min The minimum xvec value.
//' @param x_max The maximum xvec value.
//' @param q1,q2 Numeric vectors, for the lower and upper bounds of the intervals over which arc lengths are to be computed.
//' @param al_samp The sample arc length statistic.
//'
//' @return kappa4NLSobj: The nonlinear least squares objective function.
//' @export
//[[Rcpp::export]]
double kappa4NLSobj(NumericVector parms, NumericVector xvec, NumericVector y, double x_min, double x_max)
{
double x0;
int n = y.size();
NumericVector yHat(n);

if (parms[1]<=0)
{
x0 = x_min-(parms[0]/parms[2]);
for (int i=0; i<n; i++)
{
yHat[i] = pkappa4(xvec[i], x0, parms[0], parms[1], parms[2]);
}
}
else
{
x0 = x_min-parms[0]*(1-pow(parms[1], -parms[2]))/parms[2];
for (int i=0; i<n; i++)
{
yHat[i] = pkappa4(xvec[i], x0, parms[0], parms[1], parms[2]);
}
}

yHat = yHat/pkappa4(x_max, x0, parms[0], parms[1], parms[2]);

yHat = ifelse(is_nan(yHat), 1, yHat);

return pow(sum(pow(y-yHat, 2)), 0.5);
}


//' @rdname kappa4NLSobj
//' @return kappa4NLScon: A vector with three conditions evaluated.
//' @export
//[[Rcpp::export]]
NumericVector kappa4NLScon(NumericVector parms, NumericVector xvec, NumericVector y, double x_min, double x_max)
{
double x0;
NumericVector con(3);

con[0] = 0;

if (parms[1]<=0)
{
x0 = x_min-(parms[0]/parms[2]);
con[2] = (1.0/parms[1])-parms[2];
}
else
{
x0 = x_min-parms[0]*(1-pow(parms[1], -parms[2]))/parms[2];
con[2] = parms[2]-1;
if (parms[2] > 0)
{
con[0] = x0+parms[0]/parms[2]-x_max;
}
}

con[1] = dddkappa4(kappa4cond(x0, parms[0], parms[1], parms[2]), x0, parms[0], parms[1], parms[2]);

con = ifelse(is_na(con), 1, con);

return con;
}


//' @rdname kappa4NLSobj
//' @return kappa4NLShin: A vector specifying a single nonlinear inequality constraint.
//' @export
//[[Rcpp::export]]
NumericVector kappa4NLShin(NumericVector parms, NumericVector xvec, NumericVector y, double x_min, double x_max)
{
NumericVector c(1);

if (parms[1]<=0)
{
c[0] = parms[2]-(1.0/parms[1]);
}
else
{
c[0] = 1-parms[2];
}

c = ifelse(is_na(c), 1, c);

return c;
}


//' @rdname kappa4NLSobj
//' @return kappa4NLSheq: A vector specifying two nonlinear equality constraints.
//' @export
//[[Rcpp::export]]
NumericVector kappa4NLSheq(NumericVector parms, NumericVector xvec, NumericVector y, double x_min, double x_max)
{
double x0;
NumericVector c(2);

c[0] = 0;

if (parms[1]<=0)
{
x0 = x_min-(parms[0]/parms[2]);
}
else
{
x0 = x_min-parms[0]*(1-pow(parms[1], -parms[2]))/parms[2];
if (parms[2] > 0)
{
c[0] = x_max-x0-parms[0]/parms[2];
}
}

c[1] = dddkappa4(kappa4cond(x0, parms[0], parms[1], parms[2]), x0, parms[0], parms[1], parms[2]);

c = ifelse(is_na(c), 1, c);

return c;
}


//' @rdname kappa4NLSobj
//' @return kappa4ALobj: The arc length objective function.
//' @export
//[[Rcpp::export]]
double kappa4ALobj(NumericVector parms, NumericVector al_samp, double x_min, double x_max, NumericVector q1, NumericVector q2)
{
double k = parms[2];
double x0;
int n = q1.size();
NumericVector al_theo(n);

if (parms[1]<=0)
{
x0 = x_min-(parms[0]/k);
}
else
{
x0 = x_min-parms[0]*(1-pow(parms[1], -k))/k;
}

al_theo = kappa4Int2(x0, parms[0], parms[1], k, pkappa4(x_max, x0, parms[0], parms[1], k), q1, q2, false);

al_theo = ifelse(is_nan(al_theo), 0, al_theo);

return pow(sum(pow(al_theo-al_samp, 2)), 0.5);
}


//' @rdname kappa4NLSobj
//' @return kappa4ALcon: A vector with three conditions evaluated.
//' @export
//[[Rcpp::export]]
NumericVector kappa4ALcon(NumericVector parms, NumericVector al_samp, double x_min, double x_max, NumericVector q1, NumericVector q2)
{
double k = parms[2];
double x0;
NumericVector con(3);

con[0] = 0;

if (parms[1]<=0)
{
x0 = x_min-(parms[0]/k);
con[2] = (1.0/parms[1])-k;
}
else
{
x0 = x_min-parms[0]*(1-pow(parms[1], -k))/k;
con[2] = k-1;
if (k > 0)
{
con[0] = x0+parms[0]/k-x_max;
}
}

con[1] = dddkappa4(kappa4cond(x0, parms[0], parms[1], k), x0, parms[0], parms[1], k);

con = ifelse(is_na(con), 1, con);

return con;
}


//' @rdname kappa4NLSobj
//' @return kappa4ALhin: A vector specifying a single nonlinear inequality constraint.
//' @export
//[[Rcpp::export]]
NumericVector kappa4ALhin(NumericVector parms, NumericVector al_samp, double x_min, double x_max, NumericVector q1, NumericVector q2)
{
NumericVector c(1);
double k = parms[2];

if (parms[1]<=0)
{
c[0] = k-(1.0/parms[1]);
}
else
{
c[0] = 1-k;
}

c = ifelse(is_na(c), 1, c);

return c;
}


//' @rdname kappa4NLSobj
//' @return kappa4ALheq: A vector specifying two nonlinear equality constraints.
//' @export
//[[Rcpp::export]]
NumericVector kappa4ALheq(NumericVector parms, NumericVector al_samp, double x_min, double x_max, NumericVector q1, NumericVector q2)
{
double x0;
NumericVector c(2);
double k = parms[2];

c[0] = 0;

if (parms[1]<=0)
{
x0 = x_min-(parms[0]/k);
}
else
{
x0 = x_min-parms[0]*(1-pow(parms[1], -k))/k;
if (k > 0)
{
c[0] = x_max-x0-parms[0]/k;
}
}

c[1] = dddkappa4(kappa4cond(x0, parms[0], parms[1], k), x0, parms[0], parms[1], k);

c = ifelse(is_na(c), 1, c);

return c;
}