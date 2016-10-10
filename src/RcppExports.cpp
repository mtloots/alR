// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// alE
List alE(NumericVector x, NumericVector q1, NumericVector q2, bool dc, double type);
RcppExport SEXP alR_alE(SEXP xSEXP, SEXP q1SEXP, SEXP q2SEXP, SEXP dcSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q1(q1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q2(q2SEXP);
    Rcpp::traits::input_parameter< bool >::type dc(dcSEXP);
    Rcpp::traits::input_parameter< double >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(alE(x, q1, q2, dc, type));
    return rcpp_result_gen;
END_RCPP
}
// alEdist
NumericMatrix alEdist(int n, int bootstraps, double mu, double sigma, NumericVector q1, NumericVector q2, bool quantile, bool dc, double type);
RcppExport SEXP alR_alEdist(SEXP nSEXP, SEXP bootstrapsSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP q1SEXP, SEXP q2SEXP, SEXP quantileSEXP, SEXP dcSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type bootstraps(bootstrapsSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q1(q1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q2(q2SEXP);
    Rcpp::traits::input_parameter< bool >::type quantile(quantileSEXP);
    Rcpp::traits::input_parameter< bool >::type dc(dcSEXP);
    Rcpp::traits::input_parameter< double >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(alEdist(n, bootstraps, mu, sigma, q1, q2, quantile, dc, type));
    return rcpp_result_gen;
END_RCPP
}
// matrixSqrt
NumericMatrix matrixSqrt(NumericMatrix x);
RcppExport SEXP alR_matrixSqrt(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(matrixSqrt(x));
    return rcpp_result_gen;
END_RCPP
}
// alrKDE
double alrKDE(const arma::vec& beta, const arma::mat& gamma, const arma::vec& aly, const arma::vec& q1, const arma::vec& q2, const int& type);
RcppExport SEXP alR_alrKDE(SEXP betaSEXP, SEXP gammaSEXP, SEXP alySEXP, SEXP q1SEXP, SEXP q2SEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type aly(alySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type q1(q1SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type q2(q2SEXP);
    Rcpp::traits::input_parameter< const int& >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(alrKDE(beta, gamma, aly, q1, q2, type));
    return rcpp_result_gen;
END_RCPP
}
// Silverman
double Silverman(NumericVector x);
RcppExport SEXP alR_Silverman(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(Silverman(x));
    return rcpp_result_gen;
END_RCPP
}
// Silverman2
double Silverman2(NumericVector x);
RcppExport SEXP alR_Silverman2(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(Silverman2(x));
    return rcpp_result_gen;
END_RCPP
}
// bw
double bw(NumericVector x, double type);
RcppExport SEXP alR_bw(SEXP xSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(bw(x, type));
    return rcpp_result_gen;
END_RCPP
}
// GaussInt
List GaussInt(double mu, double sigma, double q1, double q2, bool quantile);
RcppExport SEXP alR_GaussInt(SEXP muSEXP, SEXP sigmaSEXP, SEXP q1SEXP, SEXP q2SEXP, SEXP quantileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type q1(q1SEXP);
    Rcpp::traits::input_parameter< double >::type q2(q2SEXP);
    Rcpp::traits::input_parameter< bool >::type quantile(quantileSEXP);
    rcpp_result_gen = Rcpp::wrap(GaussInt(mu, sigma, q1, q2, quantile));
    return rcpp_result_gen;
END_RCPP
}
// GaussInt2
NumericVector GaussInt2(double mu, double sigma, NumericVector q1, NumericVector q2, bool quantile);
RcppExport SEXP alR_GaussInt2(SEXP muSEXP, SEXP sigmaSEXP, SEXP q1SEXP, SEXP q2SEXP, SEXP quantileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q1(q1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q2(q2SEXP);
    Rcpp::traits::input_parameter< bool >::type quantile(quantileSEXP);
    rcpp_result_gen = Rcpp::wrap(GaussInt2(mu, sigma, q1, q2, quantile));
    return rcpp_result_gen;
END_RCPP
}
// dkdeGauss
double dkdeGauss(double x, NumericVector mu, double h);
RcppExport SEXP alR_dkdeGauss(SEXP xSEXP, SEXP muSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(dkdeGauss(x, mu, h));
    return rcpp_result_gen;
END_RCPP
}
// pkdeGauss
double pkdeGauss(double x, NumericVector mu, double h);
RcppExport SEXP alR_pkdeGauss(SEXP xSEXP, SEXP muSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(pkdeGauss(x, mu, h));
    return rcpp_result_gen;
END_RCPP
}
// qkdeGauss
List qkdeGauss(double x, NumericVector mu, double h);
RcppExport SEXP alR_qkdeGauss(SEXP xSEXP, SEXP muSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(qkdeGauss(x, mu, h));
    return rcpp_result_gen;
END_RCPP
}
// kdeGaussInt
List kdeGaussInt(NumericVector mu, double h, double q1, double q2, bool quantile);
RcppExport SEXP alR_kdeGaussInt(SEXP muSEXP, SEXP hSEXP, SEXP q1SEXP, SEXP q2SEXP, SEXP quantileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type q1(q1SEXP);
    Rcpp::traits::input_parameter< double >::type q2(q2SEXP);
    Rcpp::traits::input_parameter< bool >::type quantile(quantileSEXP);
    rcpp_result_gen = Rcpp::wrap(kdeGaussInt(mu, h, q1, q2, quantile));
    return rcpp_result_gen;
END_RCPP
}
// kdeGaussInt2
NumericVector kdeGaussInt2(NumericVector mu, double h, NumericVector q1, NumericVector q2, bool quantile);
RcppExport SEXP alR_kdeGaussInt2(SEXP muSEXP, SEXP hSEXP, SEXP q1SEXP, SEXP q2SEXP, SEXP quantileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q1(q1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q2(q2SEXP);
    Rcpp::traits::input_parameter< bool >::type quantile(quantileSEXP);
    rcpp_result_gen = Rcpp::wrap(kdeGaussInt2(mu, h, q1, q2, quantile));
    return rcpp_result_gen;
END_RCPP
}
// kdeGaussIntApprox
double kdeGaussIntApprox(NumericVector mu, double h, double q1, double q2, bool quantile);
RcppExport SEXP alR_kdeGaussIntApprox(SEXP muSEXP, SEXP hSEXP, SEXP q1SEXP, SEXP q2SEXP, SEXP quantileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type q1(q1SEXP);
    Rcpp::traits::input_parameter< double >::type q2(q2SEXP);
    Rcpp::traits::input_parameter< bool >::type quantile(quantileSEXP);
    rcpp_result_gen = Rcpp::wrap(kdeGaussIntApprox(mu, h, q1, q2, quantile));
    return rcpp_result_gen;
END_RCPP
}
// kdeGaussIntApprox2
NumericVector kdeGaussIntApprox2(NumericVector mu, double h, NumericVector q1, NumericVector q2, bool quantile);
RcppExport SEXP alR_kdeGaussIntApprox2(SEXP muSEXP, SEXP hSEXP, SEXP q1SEXP, SEXP q2SEXP, SEXP quantileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q1(q1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q2(q2SEXP);
    Rcpp::traits::input_parameter< bool >::type quantile(quantileSEXP);
    rcpp_result_gen = Rcpp::wrap(kdeGaussIntApprox2(mu, h, q1, q2, quantile));
    return rcpp_result_gen;
END_RCPP
}
// momKDE
double momKDE(const arma::vec& beta, const arma::mat& gamma, const arma::vec& momy, Function kdeGaussMom, const int& type);
RcppExport SEXP alR_momKDE(SEXP betaSEXP, SEXP gammaSEXP, SEXP momySEXP, SEXP kdeGaussMomSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type momy(momySEXP);
    Rcpp::traits::input_parameter< Function >::type kdeGaussMom(kdeGaussMomSEXP);
    Rcpp::traits::input_parameter< const int& >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(momKDE(beta, gamma, momy, kdeGaussMom, type));
    return rcpp_result_gen;
END_RCPP
}
// qsamp
double qsamp(NumericVector x, double q);
RcppExport SEXP alR_qsamp(SEXP xSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(qsamp(x, q));
    return rcpp_result_gen;
END_RCPP
}
