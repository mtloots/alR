#ifndef alR_H
#define alR_H

typedef double brents_fn(double x, void *ex);
typedef void integr_fn(double *t, int n, void *ex);
typedef double optimfn(int n, double *par, void *ex);

Rcpp::List brents_fun(brents_fn f, void *ex, double lower, double upper, double tol, unsigned int max_iter);

double bw(NumericVector x, double type);

Rcpp::List pairSort(Rcpp::NumericVector x);

double pkappa4(double x, double mu, double sigma, double h, double k);

double dkappa4(double x, double mu, double sigma, double h, double k);

double qkappa4(double x, double mu, double sigma, double h, double k);

double dddkappa4(double x, double mu, double sigma, double h, double k);

double kappa4cond(double mu, double sigma, double h, double k);

NumericVector kappa4Int2(double mu, double sigma, double h, double k, double tau, NumericVector q1, NumericVector q2, bool quantile);

Rcpp::List kappa4tc(double h, double mu, double sigma);


double dkdeGauss(double x, Rcpp::NumericVector mu, double h);

Rcpp::NumericVector GaussInt2(double mu, double sigma, Rcpp::NumericVector q1, Rcpp::NumericVector q2, bool quantile);

Rcpp::NumericVector kdeGaussInt2(Rcpp::NumericVector mu, double h, Rcpp::NumericVector q1, Rcpp::NumericVector q2, bool quantile);

Rcpp::NumericVector kdeGaussIntApprox2(Rcpp::NumericVector mu, double h, Rcpp::NumericVector q1, Rcpp::NumericVector q2, bool quantile);

Rcpp::NumericVector lulu(Rcpp::NumericVector x);

Rcpp::List qkdeGauss(double x, Rcpp::NumericVector mu, double h);

double qlin(NumericVector x, NumericVector y, double q);

double qsamp(Rcpp::NumericVector x, double q);

Rcpp::List Rcpp_integrate(integr_fn f, void *ex, double lower, double upper, int subdiv = 100, double eps_abs = 1e-10, double eps_rel = 1e-10);

Rcpp::List Rcpp_nmmin(int n, optimfn fn, double *xin, void *ex, double Fmin = 0.0, double abstol=1e-15, double intol=1e-15, double alpha=1.0, double beta=0.5, double gamma=2.0, int trace=0, int maxit=1000);

Rcpp::NumericVector RcppSample(Rcpp::NumericVector sample, int n);

#endif