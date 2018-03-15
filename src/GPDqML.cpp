#include <Rcpp.h>
using namespace Rcpp;
#include "alR.h"


typedef struct { NumericVector x; } Aux;


double nllGPD(int n, double *par, void *ex)
{
Aux *aux = (Aux *) ex;
NumericVector x = aux->x;
return -sum(Rcpp::log(dGPD(x, min(x), std::abs(par[0]), par[1])));
}


double mseGPD(int n, double *par, void *ex)
{
Aux *aux = (Aux *) ex;
NumericVector x = aux->x;
std::sort(x.begin(), x.end());

return -(sum(Rcpp::log(diff(pGPD(x, par[0], std::abs(par[1]), par[2])))))/(x.size()+1);
}


//' The qML Method for the Generalised Pareto Distribution (GPD).
//'
//' Estimation of the parameters of the generalised Pareto distribution (GPD) with location parameter \code{mu}, scale parameter \code{sigma}, and shape parameter \code{alpha}.  using the quasi-maximum-likelihood method.
//'
//' @param x A vector of data points.
//' @rdname GPDqML
//' @return GPDMLE:  A named vector containing the maximum likelihood estimates for the GPD fitted to the data \code{x}.
//' @examples
//' x <- rGPD(100, 0, 1, 2)
//' GPDMLE(x)
//' @export
// [[Rcpp::export]]
NumericVector GPDMLE(NumericVector x)
{
Aux aux = {x};
double xin[2] = {mean(x), 0};
NumericVector mle = Rcpp_nmmin(2, nllGPD, xin, &aux)["par"];
NumericVector par(3);

par[0]=min(x);
par[1]=mle[0];
par[2]=mle[1];

par.names() = CharacterVector::create("mu", "sigma", "alpha");

return par;
}


//' @rdname GPDqML
//' @return GPDqML:  A named vector containing the estimated parameters for the GPD fitted to the data \code{x}, using quasi-maximum likelihood.
//' @examples
//' GPDqML(x)
//' @export
// [[Rcpp::export]]
NumericVector GPDqML(NumericVector x)
{
double mu=min(x), sigma, alpha, z;
NumericVector par(3);

NumericVector a = Rcpp::log(1-(x-mu)/(max(x)-mu));
alpha = -1.0/(x.size()-1)*sum(as<NumericVector>(a[!is_infinite(a)]));

sigma = alpha*(max(x)-mu);

z = 1-(sum(pow(x, 2))/x.size())/(2*pow(mean(x), 2));

if (alpha<0.75 && z<0.2)
{
par = GPDMLE(x);
}
else
{
par[0] = mu;
par[1] = sigma;
par[2] = alpha;
}

par.names() = CharacterVector::create("mu", "sigma", "alpha");

return par;
}


//' @rdname GPDqML
//' @return GPDMSE:  A named vector containing the estimated parameters for the GPD fitted to the data \code{x}, using maximum spacings estimation (MSE).
//' @examples
//' \dontrun{
//' GPDMSE(x)
//' }
//' @export
// [[Rcpp::export]]
NumericVector GPDMSE(NumericVector x)
{
Aux aux = {x};
NumericVector start = GPDqML(x);
double xin[3] = {start[0], start[1], start[2]};
NumericVector mse = Rcpp_nmmin(3, mseGPD, xin, &aux)["par"];

mse.names() = CharacterVector::create("mu", "sigma", "alpha");

return mse;
}