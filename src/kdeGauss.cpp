#include <Rcpp.h>
#include <R_ext/Applic.h>
using namespace Rcpp;


typedef struct { NumericVector mu; double h; } Params;
typedef double brents_fn(double x, void *ex);
typedef struct { NumericVector mu; double h; double q; } Qarams;


List Rcpp_integrate_b(integr_fn f, void *ex, double bound, int inf = -1, int subdiv = 100, double eps_abs = 1e-10, double eps_rel = 1e-10)
{
int lenw = 4 * subdiv;
int *iwork = new int[subdiv];
double *work = new double[lenw];
double value;
double abs_err;
int subdiv_used;
int neval;
int error;

Rdqagi(f, ex, &bound, &inf, &eps_abs, &eps_rel, &value, &abs_err, &neval, &error, &subdiv, &lenw, &subdiv_used, iwork, work);

delete [] iwork;
delete [] work;

return List::create(_["value"] = value,
_["abs.err"] = abs_err,
_["subdivisions"] = subdiv_used,
_["neval"] = neval);
}


List brents_fun(brents_fn f, void *ex, double lower, double upper, double tol, unsigned int max_iter)
{
int msg = 0;
unsigned int iterations = 0;
double a = lower;
double b = upper;
double fa = f(a, ex); // calculated now to save function calls
double fb = f(b, ex); // calculated now to save function calls
double fs = 0; // initialize 
 
if (fa*fb >= 0)
{
msg = -1; //Signs of f(lower_bound) and f(upper_bound) must be opposites
}
 
if (std::abs(fa) < std::abs(b)) // if magnitude of f(lower_bound) is less than magnitude of f(upper_bound)
{
std::swap(a,b);
std::swap(fa,fb);
}
 
double c = a; // c now equals the largest magnitude of the lower and upper bounds
double fc = fa; // precompute function evalutation for point c by assigning it the same value as fa
bool mflag = true; // boolean flag used to evaluate if statement later on
double s = 0; // Our Root that will be returned
double d = 0; // Only used if mflag is unset (mflag == false)
 
if (msg == 0)
{
for (unsigned int iter = 1; iter < max_iter; ++iter)
{
// stop if converged on root or error is less than tolerance
if (std::abs(b-a) < tol)
{
iterations = iter-1;
iter = max_iter;
} // end if
 
if (fa != fc && fb != fc)
{
// use inverse quadratic interpolation
s =  ( a * fb * fc / ((fa - fb) * (fa - fc)) )
+ ( b * fa * fc / ((fb - fa) * (fb - fc)) )
+ ( c * fa * fb / ((fc - fa) * (fc - fb)) );
}
else
{
// secant method
s = b - fb * (b - a) / (fb - fa);
}
 
// checks to see whether we can use the faster converging quadratic && secant methods or if we need to use bisection
if (( (s < (3 * a + b) * 0.25) || (s > b) ) ||
( mflag && (std::abs(s-b) >= (std::abs(b-c) * 0.5)) ) ||
( !mflag && (std::abs(s-b) >= (std::abs(c-d) * 0.5)) ) ||
( mflag && (std::abs(b-c) < tol) ) ||
( !mflag && (std::abs(c-d) < tol)))
{
// bisection method
s = (a+b)*0.5;
 
mflag = true;
}
else
{
mflag = false;
}
 
fs = f(s, ex); // calculate fs
d = c; // first time d is being used (wasnt used on first iteration because mflag was set)
c = b; // set c equal to upper bound
fc = fb; // set f(c) = f(b)
 
if ( fa * fs < 0)// fa and fs have opposite signs
{
b = s;
fb = fs; // set f(b) = f(s)
}
else
{
a = s;
fa = fs; // set f(a) = f(s)
}
 
if (std::abs(fa) < std::abs(fb)) // if magnitude of fa is less than magnitude of fb
{
std::swap(a,b); // swap a and b
std::swap(fa,fb); // make sure f(a) and f(b) are correct after swap
}
 
} // end for
 } // end if (msg == 0)

return List::create(_["result"] = s,
_["value"] = f(s, ex),
_["iterations"] = iterations,
_["msg"] = msg);
 } // end brents_fun
 

//' Gaussian kernel density estimator.
//'
//' Estimate a density function using a kernel density estimator with a Gaussian kernel.
//'
//' The cumulative distribution function is calculated using the numerical integration C code implimented for R's integrate functions, i.e. using Rdqagi.  For this approximation, subdiv = 100 (100 subdivisions), and eps_abs = eps_rel = 1e-10, i.e. the absolute and relative errors respectively.
//'
//' The quantiles of the Gaussian kernel density estimator are calculated using Brent's method.  This method requires an interval in which a solution is saught.  The objective funcion for which a zero is saught is \code{\link{dkdeGauss}}-\code{x}, where \code{x} is the quantile saught.  The first interval in which a solution is searched for, corresponds to the range of \code{mu}, and is expanded in multiples thereof in consequtive steps.  The maximum number of iterations is set at 1000, and the accuracy saught between iterations, is set at 1e-10.
//'
//' @param x A data point, or quantile, at which the kernel density estimator should be evaluated.
//' @param mu A vector of data points on which the kernel density estimator is based.
//' @param h The kernel density estimator bandwidth.
//' @rdname kdeGauss
//' @return pkdeGauss: The estimated value of the density function at the point x.
//' @examples
//' library(alR)
//' x <- rnorm(100)
//' h_x <- Silverman(x)
//' pkdeGauss(0, x, h_x)
//' @export
// [[Rcpp::export]]
double pkdeGauss(double x, NumericVector mu, double h)
{
return (1.0/(mu.size()*sqrt(2*PI)*h))*sum(exp(-0.5*pow((x-mu)/h, 2)));
}


void pdfkdeGauss(double *t, int n, void *ex)
{
Params *param = (Params *) ex;
NumericVector mu = param->mu;
double h = param->h;

for(int i = 0; i < n; i++)
{
t[i] = pkdeGauss(t[i], mu, h);
}
} 


//' @rdname kdeGauss
//' @examples
//' dkdeGauss(0, x, h_x)
//' @return dkdeGauss: A list with the following components:
//' \itemize{
//' \item value: The estimated value of the cumulative distribution function at the point \code{x}.
//' \item abs.err: The absolute error between iterations.
//' subdivisions: Number of subdivisions used in the numerical approximation.
//' \item neval: Number of function evaluations used by the numerical approximation.
//' }
//' @export
// [[Rcpp::export]]
List dkdeGauss(double x, NumericVector mu, double h)
{
Params param = {mu, h};

return Rcpp_integrate_b(pdfkdeGauss, &param, x);
}


double qkde(double x, void *ex)
{
Qarams *param = (Qarams *) ex;
NumericVector mu = param->mu;
double h = param->h;
double q = param->q;

return double(dkdeGauss(x, mu, h)[0])-q;
}


//' @rdname kdeGauss
//' @examples
//' qkdeGauss(0.5, x, h_x)
//' @return qkdeGauss: A list with the following components:
//' \itemize{
//' \item result: The \code{x}th quantile of the Gaussian kernel density estimator.
//' \item value: The value of the cumulative distribution function of the Gaussian kernel density estimator at the \code{x}th quantile.
//' \item obj.fun: The value of the objective function resulting from Brent's method; should be less than 1e-10.
//' \item iterations: Number of iterations for Brent's method in order to achieve the desired accuracy.
//' \item steps: Number of range expansions of the search boundaries for Brent's method.
//' }
//' @export
// [[Rcpp::export]]
List qkdeGauss(double x, NumericVector mu, double h)
{
Qarams param = {mu, h, x};
double minmu = min(mu);
double maxmu = max(mu);
int steps = 0;
double lower, upper;
List bf;

for (int i=1; i<=1000; i++)
{
steps++;
lower = steps*minmu;
upper = steps*maxmu;
bf = brents_fun(qkde, &param, lower, upper, 1e-10, 1000);
if (bf[3] == 0)
{
i = 1001;
}
}

return List::create(_["result"] = bf[0],
_["value"] = dkdeGauss(bf[0], mu, h)[0],
_["obj.fun"] = bf[1],
_["iterations"] = bf[2],
_["steps"] = steps);
}