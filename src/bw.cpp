#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]] 

//' Silverman's rule of thumb.
//'
//' Calculate the bandwidth estimator using Silverman's rule of thumb.
//'
//' @param x A vector of data points.
//'
//' @return Bandwidth estimator based on Silverman's rule of thumb.
//'
//' @examples
//' set.seed(1)
//' x <- rnorm(100)
//' Silverman(x)
//'
//' @export
// [[Rcpp::export]]
double Silverman(const arma::vec& x)
{
return arma::stddev(x)*pow(4.0/(3*x.n_rows), 0.2);
}

//' Silverman's adapted rule of thumb.
//'
//' Calculate the bandwidth estimator using Silverman's adapted rule of thumb.
//'
//' @param x A vector of data points.
//'
//' @return Bandwidth estimator based on Silverman's adapted rule of thumb.
//'
//' @examples
//' set.seed(1)
//' x <- rnorm(100)
//' Silverman2(x)
//'
//' @export
//[[Rcpp::export]]
double Silverman2(const arma::vec&  x)
{
int n = x.n_rows;
double y25, y75;
NumericVector y = NumericVector(x.begin(), x.end());

// Calculate 0.25th empirical quantile
int q25 = floor(n*0.25);
double g25 = n*0.25-q25;
std::nth_element(y.begin(), y.begin()+q25, y.end());
if (g25 == 0)
{
y25 = y[q25];
}
else
{
y25 = y[q25+1];
}

// Calculate 0.75th empirical quantile
int q75 = floor(n*0.75);
double g75 = n*0.75-q75;
std::nth_element(y.begin(), y.begin()+q75, y.end());
if (g75 == 0)
{
y75 = y[q75];
}
else
{
y75 = y[q75+1];
}

double A = std::min(arma::stddev(x), (y75-y25)/1.34);
return 0.9*A/pow(n, 0.2);
}

//' Kernel density bandwidth parameters.
//'
//' Calculate the bandwidth estimator using various methods.
//'
//' @param x A vector of data points.
//' @param type An integer:
//' \itemize{
//' \item 1: Silverman's rule of thumb.
//' \item 2: Silverman's adapted rule of thumb.
//' \item default: Silvermans rule of thumb (number 1).
//' }
//'
//' @return Bandwidth estimator based on Silverman's rule of thumb.
//'
//' @examples
//' set.seed(1)
//' x <- rnorm(100)
//' bw(x, 1)
//' bw(x, 2)
//'
//' @export
//[[Rcpp::export]]
double bw(const arma::vec& x, const int& type)
{
switch(type)
{
case 1:
return Silverman(x);
break;
case 2:
return Silverman2(x);
break;
default:
return Silverman(x);
}
}