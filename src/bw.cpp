#include <Rcpp.h>
using namespace Rcpp;
#include "alR.h"


//' @rdname bw
//' @return Silverman: Bandwidth estimator based on Silverman's rule of thumb.
//'
//' @examples
//' set.seed(1)
//' x <- rnorm(100)
//' Silverman(x)
//'
//' @export
// [[Rcpp::export]]
double Silverman(NumericVector x)
{
return sd(x)*pow(4.0/(3*x.size()), 0.2);
}

//' @rdname bw
//' @return Silverman2: Bandwidth estimator based on Silverman's adapted rule of thumb.
//'
//' @examples
//' Silverman2(x)
//'
//' @export
//[[Rcpp::export]]
double Silverman2(NumericVector x)
{
double x25 = qsamp(x, 0.25);
double x75 = qsamp(x, 0.75);

double A = std::min(double(sd(x)), (x75-x25)/1.34);
return 0.9*A/pow(x.size(), 0.2);
}


//' Kernel density bandwidth estimators.
//'
//' Calculate the bandwidth estimator using various methods.
//'
//' @param x A vector of data points.
//' @param type One of the following options:
//' \itemize{
//' \item -1: Silverman's rule of thumb.
//' \item -2: Silverman's adapted rule of thumb.
//' \item >0: The real number is returned without any calculations.
//' }
//'
//' @return bw: Bandwidth estimator based on the selected method.
//'
//' @examples
//' bw(x, -1)
//' bw(x, -2)
//' bw(x, 0.5)
//'
//' @export
//[[Rcpp::export]]
double bw(NumericVector x, double type)
{
if(type == -1) return Silverman(x);
else if(type == -2) return Silverman2(x);
else return type;
}