#include <Rcpp.h>
#include <R_ext/Applic.h>
using namespace Rcpp;


typedef double brents_fn(double x, void *ex);


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


//' Point corresponding to quantile.
//' Calculate the \eqn{x} point corresponding to a quantile \eqn{F(x)} using linear interpolation.
//'
//' @param x A numeric vector, specifying the \eqn{x} values.
//' @param y A numeric vector, specifying the \eqn{F(x)} values.
//' @param q A real number between 0 and 1 inclusive, specifying the desired quantile.
//'
//' @return The interpolated quantile, \eqn{x}, corresponding to \eqn{q=F(x)}.
//' @examples
//' x <- rnorm(100)
//' y <- pnorm(x)
//' qlin(x, y, 0.5)
//' @export
// [[Rcpp::export]]
double qlin(NumericVector x, NumericVector y, double q)
{
int pos2 = sum(ifelse(y<=q, 1, 0));

if (pos2 == 0)
{
pos2 = 1;
}

int pos1 = pos2-1;

std::nth_element(y.begin(), y.begin()+pos2, y.end());
std::nth_element(y.begin(), y.begin()+pos1, y.end());

double x1 = x[pos1];
double x2 = x[pos2];
double y1 = y[pos1];
double y2 = y[pos2];

double m = (y2-y1)/(x2-x1);
double c = y1-m*x1;

return (q-c)/m;
}


//' Empirical sample quantile.
//' Calculate empirical sample quantile.
//'
//' @param x A numeric vector, specifying the sample for which the quantile is to be calculated.
//' @param q A real number between 0 and 1 inclusive, specifying the desired quantile.
//'
//' @return The empirical quantile of the provided sample.
//' @examples
//' x<-rnorm(100)
//' qsamp(x, 0.5)
//' @export
// [[Rcpp::export]]
double qsamp(NumericVector x, double q)
{
int n = x.size();
int Q = floor(n*q);
double g = n*q-double(Q);
if (g == 0) Q -= 1; //Q = Q, but c++ uses 0-based indexing
else Q = Q; //Q -= 1 0-based indexing

std::nth_element(x.begin(), x.begin()+Q, x.end());

return x[Q];
}


List Rcpp_integrate(integr_fn f, void *ex, double lower, double upper, int subdiv = 100, double eps_abs = 1e-10, double eps_rel = 1e-10)
{
int lenw = 4 * subdiv;
int *iwork = new int[subdiv];
double *work = new double[lenw];
double value;
double abs_err;
int subdiv_used;
int neval;
int error;

Rdqags(f, ex, &lower, &upper, &eps_abs, &eps_rel, &value, &abs_err, &neval, &error, &subdiv, &lenw, &subdiv_used, iwork, work);

delete [] iwork;
delete [] work;

return List::create(_["value"] = value,
_["abs.err"] = abs_err,
_["subdivisions"] = subdiv_used,
_["neval"] = neval);
}


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


List Rcpp_nmmin(int n, optimfn fn, double *xin, void *ex, double Fmin = 0.0, double abstol=1e-15, double intol=1e-15, double alpha=1.0, double beta=0.5, double gamma=2.0, int trace=0, int maxit=1000)
{
double x;
double *p;
int fail;
int fncount;
NumericVector X(n);

nmmin(n, xin, &x, &Fmin, fn, &fail, abstol, intol, ex, alpha, beta, gamma, trace, &fncount, maxit);

p = &x;

for (int i=0; i<n; i++)
{
X[i] = *(p+i);
}

return List::create(_["par"] = X,
_["abstol"] = abstol,
_["fail"] = fail,
_["fncount"] = fncount);
}

NumericVector RcppSample(NumericVector sample, int n)
{
return sample[round(runif(n)*(sample.size()-1.0), 0)];
}


//' Sorted vector index.
//'
//' The sorted vector is returned along with the original index of the vector it belonged to.
//'
//' @param x A real-valued vector to be sorted.
//'
//' @return A list with two components:
//' \itemize{
//' \item sorted: The sorted version of \code{x}.
//' \item index: The index of the \eqn{i^{th}} element in \code{x}.
//' }
//'
//' @examples
//' pairSort(c(5, 2, 6))
//' @export
//[[Rcpp::export]]
List pairSort(NumericVector x)
{
int n = x.size();
NumericVector x_sorted(n);
IntegerVector x_index(n);

std::vector<double> arr(x.begin(), x.end());
std::vector<std::pair<double,int> >V;

for (int i=0; i<x.size(); i++)
{
std::pair<double,int>P = std::make_pair(arr[i], i);
V.push_back(P);
}

std::sort(V.begin(), V.end());

for (int i=0; i<x.size(); i++)
{
x_sorted[i] = V[i].first;
x_index[i] = V[i].second;
}

return List::create(_["sorted"] = x_sorted,
_["index"] = x_index);
}


//' LULU smoother.
//'
//' Performs LULU smoothing of the provided vector.
//'
//' @param x A real-valued vector.
//'
//' @return The LULU-smoothed version of \code{x}.
//'
//' @examples
//' x <- rnorm(10)
//' lulu(x)
//' @export
//[[Rcpp::export]]
NumericVector lulu(NumericVector x)
{
int n = x.size();
NumericVector x_out(n);

x.insert(0,0);
x.insert(n+1,0);

for (int i=1; i<=n; i++)
{
if ((x[i] < x[i-1]) & (x[i] < x[i+1]))
{
x_out[i-1] = std::min(x[i-1], x[i+1]);
}
else if ((x[i] > x[i-1]) & (x[i] > x[i+1]))
{
x_out[i-1] = std::max(x[i-1], x[i+1]);
}
else
{
x_out[i-1] = x[i];
}
}

return x_out;
}