// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// alE
List alE(NumericVector x, NumericVector q1, NumericVector q2, bool dc, double type, int distribution);
RcppExport SEXP _alR_alE(SEXP xSEXP, SEXP q1SEXP, SEXP q2SEXP, SEXP dcSEXP, SEXP typeSEXP, SEXP distributionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q1(q1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q2(q2SEXP);
    Rcpp::traits::input_parameter< bool >::type dc(dcSEXP);
    Rcpp::traits::input_parameter< double >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type distribution(distributionSEXP);
    rcpp_result_gen = Rcpp::wrap(alE(x, q1, q2, dc, type, distribution));
    return rcpp_result_gen;
END_RCPP
}
// alEfitdist
NumericMatrix alEfitdist(NumericVector x, NumericVector q1, NumericVector q2, bool dc, double type, int bootstraps, int distribution);
RcppExport SEXP _alR_alEfitdist(SEXP xSEXP, SEXP q1SEXP, SEXP q2SEXP, SEXP dcSEXP, SEXP typeSEXP, SEXP bootstrapsSEXP, SEXP distributionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q1(q1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q2(q2SEXP);
    Rcpp::traits::input_parameter< bool >::type dc(dcSEXP);
    Rcpp::traits::input_parameter< double >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type bootstraps(bootstrapsSEXP);
    Rcpp::traits::input_parameter< int >::type distribution(distributionSEXP);
    rcpp_result_gen = Rcpp::wrap(alEfitdist(x, q1, q2, dc, type, bootstraps, distribution));
    return rcpp_result_gen;
END_RCPP
}
// alEdist
NumericMatrix alEdist(int n, int bootstraps, NumericVector mu, double sigma, NumericVector q1, NumericVector q2, bool quantile, bool dc, double type, int distribution);
RcppExport SEXP _alR_alEdist(SEXP nSEXP, SEXP bootstrapsSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP q1SEXP, SEXP q2SEXP, SEXP quantileSEXP, SEXP dcSEXP, SEXP typeSEXP, SEXP distributionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type bootstraps(bootstrapsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q1(q1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q2(q2SEXP);
    Rcpp::traits::input_parameter< bool >::type quantile(quantileSEXP);
    Rcpp::traits::input_parameter< bool >::type dc(dcSEXP);
    Rcpp::traits::input_parameter< double >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type distribution(distributionSEXP);
    rcpp_result_gen = Rcpp::wrap(alEdist(n, bootstraps, mu, sigma, q1, q2, quantile, dc, type, distribution));
    return rcpp_result_gen;
END_RCPP
}
// alrKDE
double alrKDE(const arma::vec& beta, const arma::mat& gamma, const arma::vec& aly, const arma::vec& q1, const arma::vec& q2, const int& type);
RcppExport SEXP _alR_alrKDE(SEXP betaSEXP, SEXP gammaSEXP, SEXP alySEXP, SEXP q1SEXP, SEXP q2SEXP, SEXP typeSEXP) {
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
RcppExport SEXP _alR_Silverman(SEXP xSEXP) {
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
RcppExport SEXP _alR_Silverman2(SEXP xSEXP) {
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
RcppExport SEXP _alR_bw(SEXP xSEXP, SEXP typeSEXP) {
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
RcppExport SEXP _alR_GaussInt(SEXP muSEXP, SEXP sigmaSEXP, SEXP q1SEXP, SEXP q2SEXP, SEXP quantileSEXP) {
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
RcppExport SEXP _alR_GaussInt2(SEXP muSEXP, SEXP sigmaSEXP, SEXP q1SEXP, SEXP q2SEXP, SEXP quantileSEXP) {
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
// dGPD
NumericVector dGPD(NumericVector x, double mu, double sigma, double alpha);
RcppExport SEXP _alR_dGPD(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(dGPD(x, mu, sigma, alpha));
    return rcpp_result_gen;
END_RCPP
}
// pGPD
NumericVector pGPD(NumericVector q, double mu, double sigma, double alpha);
RcppExport SEXP _alR_pGPD(SEXP qSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(pGPD(q, mu, sigma, alpha));
    return rcpp_result_gen;
END_RCPP
}
// qGPD
double qGPD(double p, double mu, double sigma, double alpha);
RcppExport SEXP _alR_qGPD(SEXP pSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(qGPD(p, mu, sigma, alpha));
    return rcpp_result_gen;
END_RCPP
}
// rGPD
NumericVector rGPD(const int n, double mu, double sigma, double alpha);
RcppExport SEXP _alR_rGPD(SEXP nSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(rGPD(n, mu, sigma, alpha));
    return rcpp_result_gen;
END_RCPP
}
// GPDInt
List GPDInt(double mu, double sigma, double alpha, double q1, double q2, bool quantile);
RcppExport SEXP _alR_GPDInt(SEXP muSEXP, SEXP sigmaSEXP, SEXP alphaSEXP, SEXP q1SEXP, SEXP q2SEXP, SEXP quantileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type q1(q1SEXP);
    Rcpp::traits::input_parameter< double >::type q2(q2SEXP);
    Rcpp::traits::input_parameter< bool >::type quantile(quantileSEXP);
    rcpp_result_gen = Rcpp::wrap(GPDInt(mu, sigma, alpha, q1, q2, quantile));
    return rcpp_result_gen;
END_RCPP
}
// GPDInt2
NumericVector GPDInt2(double mu, double sigma, double alpha, NumericVector q1, NumericVector q2, bool quantile);
RcppExport SEXP _alR_GPDInt2(SEXP muSEXP, SEXP sigmaSEXP, SEXP alphaSEXP, SEXP q1SEXP, SEXP q2SEXP, SEXP quantileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q1(q1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q2(q2SEXP);
    Rcpp::traits::input_parameter< bool >::type quantile(quantileSEXP);
    rcpp_result_gen = Rcpp::wrap(GPDInt2(mu, sigma, alpha, q1, q2, quantile));
    return rcpp_result_gen;
END_RCPP
}
// GPDMLE
NumericVector GPDMLE(NumericVector x);
RcppExport SEXP _alR_GPDMLE(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(GPDMLE(x));
    return rcpp_result_gen;
END_RCPP
}
// GPDqML
NumericVector GPDqML(NumericVector x);
RcppExport SEXP _alR_GPDqML(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(GPDqML(x));
    return rcpp_result_gen;
END_RCPP
}
// GPDMSE
NumericVector GPDMSE(NumericVector x);
RcppExport SEXP _alR_GPDMSE(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(GPDMSE(x));
    return rcpp_result_gen;
END_RCPP
}
// pkappa4
double pkappa4(double x, double mu, double sigma, double h, double k);
RcppExport SEXP _alR_pkappa4(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP hSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(pkappa4(x, mu, sigma, h, k));
    return rcpp_result_gen;
END_RCPP
}
// dkappa4
double dkappa4(double x, double mu, double sigma, double h, double k);
RcppExport SEXP _alR_dkappa4(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP hSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(dkappa4(x, mu, sigma, h, k));
    return rcpp_result_gen;
END_RCPP
}
// qkappa4
double qkappa4(double x, double mu, double sigma, double h, double k);
RcppExport SEXP _alR_qkappa4(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP hSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(qkappa4(x, mu, sigma, h, k));
    return rcpp_result_gen;
END_RCPP
}
// rkappa4
NumericVector rkappa4(int n, double mu, double sigma, double h, double k);
RcppExport SEXP _alR_rkappa4(SEXP nSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP hSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(rkappa4(n, mu, sigma, h, k));
    return rcpp_result_gen;
END_RCPP
}
// dddkappa4
double dddkappa4(double x, double mu, double sigma, double h, double k);
RcppExport SEXP _alR_dddkappa4(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP hSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(dddkappa4(x, mu, sigma, h, k));
    return rcpp_result_gen;
END_RCPP
}
// kappa4cond
double kappa4cond(double mu, double sigma, double h, double k);
RcppExport SEXP _alR_kappa4cond(SEXP muSEXP, SEXP sigmaSEXP, SEXP hSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(kappa4cond(mu, sigma, h, k));
    return rcpp_result_gen;
END_RCPP
}
// kappa4tc
List kappa4tc(double h, double mu, double sigma);
RcppExport SEXP _alR_kappa4tc(SEXP hSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(kappa4tc(h, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// kappa4Int
List kappa4Int(double mu, double sigma, double h, double k, double tau, double q1, double q2, bool quantile);
RcppExport SEXP _alR_kappa4Int(SEXP muSEXP, SEXP sigmaSEXP, SEXP hSEXP, SEXP kSEXP, SEXP tauSEXP, SEXP q1SEXP, SEXP q2SEXP, SEXP quantileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type q1(q1SEXP);
    Rcpp::traits::input_parameter< double >::type q2(q2SEXP);
    Rcpp::traits::input_parameter< bool >::type quantile(quantileSEXP);
    rcpp_result_gen = Rcpp::wrap(kappa4Int(mu, sigma, h, k, tau, q1, q2, quantile));
    return rcpp_result_gen;
END_RCPP
}
// kappa4Int2
NumericVector kappa4Int2(double mu, double sigma, double h, double k, double tau, NumericVector q1, NumericVector q2, bool quantile);
RcppExport SEXP _alR_kappa4Int2(SEXP muSEXP, SEXP sigmaSEXP, SEXP hSEXP, SEXP kSEXP, SEXP tauSEXP, SEXP q1SEXP, SEXP q2SEXP, SEXP quantileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q1(q1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q2(q2SEXP);
    Rcpp::traits::input_parameter< bool >::type quantile(quantileSEXP);
    rcpp_result_gen = Rcpp::wrap(kappa4Int2(mu, sigma, h, k, tau, q1, q2, quantile));
    return rcpp_result_gen;
END_RCPP
}
// kappa4IntApprox
double kappa4IntApprox(NumericVector x, NumericVector y, double q1, double q2, bool quantile);
RcppExport SEXP _alR_kappa4IntApprox(SEXP xSEXP, SEXP ySEXP, SEXP q1SEXP, SEXP q2SEXP, SEXP quantileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type q1(q1SEXP);
    Rcpp::traits::input_parameter< double >::type q2(q2SEXP);
    Rcpp::traits::input_parameter< bool >::type quantile(quantileSEXP);
    rcpp_result_gen = Rcpp::wrap(kappa4IntApprox(x, y, q1, q2, quantile));
    return rcpp_result_gen;
END_RCPP
}
// kappa4IntApprox2
NumericVector kappa4IntApprox2(NumericVector x, NumericVector y, NumericVector q1, NumericVector q2, bool quantile);
RcppExport SEXP _alR_kappa4IntApprox2(SEXP xSEXP, SEXP ySEXP, SEXP q1SEXP, SEXP q2SEXP, SEXP quantileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q1(q1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q2(q2SEXP);
    Rcpp::traits::input_parameter< bool >::type quantile(quantileSEXP);
    rcpp_result_gen = Rcpp::wrap(kappa4IntApprox2(x, y, q1, q2, quantile));
    return rcpp_result_gen;
END_RCPP
}
// kappa4NLSobj
double kappa4NLSobj(NumericVector parms, NumericVector xvec, NumericVector y, double x_min, double x_max);
RcppExport SEXP _alR_kappa4NLSobj(SEXP parmsSEXP, SEXP xvecSEXP, SEXP ySEXP, SEXP x_minSEXP, SEXP x_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type parms(parmsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xvec(xvecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type x_min(x_minSEXP);
    Rcpp::traits::input_parameter< double >::type x_max(x_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(kappa4NLSobj(parms, xvec, y, x_min, x_max));
    return rcpp_result_gen;
END_RCPP
}
// kappa4NLScon
NumericVector kappa4NLScon(NumericVector parms, NumericVector xvec, NumericVector y, double x_min, double x_max);
RcppExport SEXP _alR_kappa4NLScon(SEXP parmsSEXP, SEXP xvecSEXP, SEXP ySEXP, SEXP x_minSEXP, SEXP x_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type parms(parmsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xvec(xvecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type x_min(x_minSEXP);
    Rcpp::traits::input_parameter< double >::type x_max(x_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(kappa4NLScon(parms, xvec, y, x_min, x_max));
    return rcpp_result_gen;
END_RCPP
}
// kappa4NLShin
NumericVector kappa4NLShin(NumericVector parms, NumericVector xvec, NumericVector y, double x_min, double x_max);
RcppExport SEXP _alR_kappa4NLShin(SEXP parmsSEXP, SEXP xvecSEXP, SEXP ySEXP, SEXP x_minSEXP, SEXP x_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type parms(parmsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xvec(xvecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type x_min(x_minSEXP);
    Rcpp::traits::input_parameter< double >::type x_max(x_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(kappa4NLShin(parms, xvec, y, x_min, x_max));
    return rcpp_result_gen;
END_RCPP
}
// kappa4NLSheq
NumericVector kappa4NLSheq(NumericVector parms, NumericVector xvec, NumericVector y, double x_min, double x_max);
RcppExport SEXP _alR_kappa4NLSheq(SEXP parmsSEXP, SEXP xvecSEXP, SEXP ySEXP, SEXP x_minSEXP, SEXP x_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type parms(parmsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xvec(xvecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type x_min(x_minSEXP);
    Rcpp::traits::input_parameter< double >::type x_max(x_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(kappa4NLSheq(parms, xvec, y, x_min, x_max));
    return rcpp_result_gen;
END_RCPP
}
// kappa4ALobj
double kappa4ALobj(NumericVector parms, NumericVector al_samp, double x_min, double x_max, NumericVector q1, NumericVector q2);
RcppExport SEXP _alR_kappa4ALobj(SEXP parmsSEXP, SEXP al_sampSEXP, SEXP x_minSEXP, SEXP x_maxSEXP, SEXP q1SEXP, SEXP q2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type parms(parmsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type al_samp(al_sampSEXP);
    Rcpp::traits::input_parameter< double >::type x_min(x_minSEXP);
    Rcpp::traits::input_parameter< double >::type x_max(x_maxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q1(q1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q2(q2SEXP);
    rcpp_result_gen = Rcpp::wrap(kappa4ALobj(parms, al_samp, x_min, x_max, q1, q2));
    return rcpp_result_gen;
END_RCPP
}
// kappa4ALcon
NumericVector kappa4ALcon(NumericVector parms, NumericVector al_samp, double x_min, double x_max, NumericVector q1, NumericVector q2);
RcppExport SEXP _alR_kappa4ALcon(SEXP parmsSEXP, SEXP al_sampSEXP, SEXP x_minSEXP, SEXP x_maxSEXP, SEXP q1SEXP, SEXP q2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type parms(parmsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type al_samp(al_sampSEXP);
    Rcpp::traits::input_parameter< double >::type x_min(x_minSEXP);
    Rcpp::traits::input_parameter< double >::type x_max(x_maxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q1(q1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q2(q2SEXP);
    rcpp_result_gen = Rcpp::wrap(kappa4ALcon(parms, al_samp, x_min, x_max, q1, q2));
    return rcpp_result_gen;
END_RCPP
}
// kappa4ALhin
NumericVector kappa4ALhin(NumericVector parms, NumericVector al_samp, double x_min, double x_max, NumericVector q1, NumericVector q2);
RcppExport SEXP _alR_kappa4ALhin(SEXP parmsSEXP, SEXP al_sampSEXP, SEXP x_minSEXP, SEXP x_maxSEXP, SEXP q1SEXP, SEXP q2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type parms(parmsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type al_samp(al_sampSEXP);
    Rcpp::traits::input_parameter< double >::type x_min(x_minSEXP);
    Rcpp::traits::input_parameter< double >::type x_max(x_maxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q1(q1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q2(q2SEXP);
    rcpp_result_gen = Rcpp::wrap(kappa4ALhin(parms, al_samp, x_min, x_max, q1, q2));
    return rcpp_result_gen;
END_RCPP
}
// kappa4ALheq
NumericVector kappa4ALheq(NumericVector parms, NumericVector al_samp, double x_min, double x_max, NumericVector q1, NumericVector q2);
RcppExport SEXP _alR_kappa4ALheq(SEXP parmsSEXP, SEXP al_sampSEXP, SEXP x_minSEXP, SEXP x_maxSEXP, SEXP q1SEXP, SEXP q2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type parms(parmsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type al_samp(al_sampSEXP);
    Rcpp::traits::input_parameter< double >::type x_min(x_minSEXP);
    Rcpp::traits::input_parameter< double >::type x_max(x_maxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q1(q1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q2(q2SEXP);
    rcpp_result_gen = Rcpp::wrap(kappa4ALheq(parms, al_samp, x_min, x_max, q1, q2));
    return rcpp_result_gen;
END_RCPP
}
// dkdeGauss
double dkdeGauss(double x, NumericVector mu, double h);
RcppExport SEXP _alR_dkdeGauss(SEXP xSEXP, SEXP muSEXP, SEXP hSEXP) {
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
RcppExport SEXP _alR_pkdeGauss(SEXP xSEXP, SEXP muSEXP, SEXP hSEXP) {
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
RcppExport SEXP _alR_qkdeGauss(SEXP xSEXP, SEXP muSEXP, SEXP hSEXP) {
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
RcppExport SEXP _alR_kdeGaussInt(SEXP muSEXP, SEXP hSEXP, SEXP q1SEXP, SEXP q2SEXP, SEXP quantileSEXP) {
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
RcppExport SEXP _alR_kdeGaussInt2(SEXP muSEXP, SEXP hSEXP, SEXP q1SEXP, SEXP q2SEXP, SEXP quantileSEXP) {
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
RcppExport SEXP _alR_kdeGaussIntApprox(SEXP muSEXP, SEXP hSEXP, SEXP q1SEXP, SEXP q2SEXP, SEXP quantileSEXP) {
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
RcppExport SEXP _alR_kdeGaussIntApprox2(SEXP muSEXP, SEXP hSEXP, SEXP q1SEXP, SEXP q2SEXP, SEXP quantileSEXP) {
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
RcppExport SEXP _alR_momKDE(SEXP betaSEXP, SEXP gammaSEXP, SEXP momySEXP, SEXP kdeGaussMomSEXP, SEXP typeSEXP) {
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
// qlin
double qlin(NumericVector x, NumericVector y, double q);
RcppExport SEXP _alR_qlin(SEXP xSEXP, SEXP ySEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(qlin(x, y, q));
    return rcpp_result_gen;
END_RCPP
}
// qsamp
double qsamp(NumericVector x, double q);
RcppExport SEXP _alR_qsamp(SEXP xSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(qsamp(x, q));
    return rcpp_result_gen;
END_RCPP
}
// pairSort
List pairSort(NumericVector x);
RcppExport SEXP _alR_pairSort(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(pairSort(x));
    return rcpp_result_gen;
END_RCPP
}
// lulu
NumericVector lulu(NumericVector x);
RcppExport SEXP _alR_lulu(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(lulu(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_alR_alE", (DL_FUNC) &_alR_alE, 6},
    {"_alR_alEfitdist", (DL_FUNC) &_alR_alEfitdist, 7},
    {"_alR_alEdist", (DL_FUNC) &_alR_alEdist, 10},
    {"_alR_alrKDE", (DL_FUNC) &_alR_alrKDE, 6},
    {"_alR_Silverman", (DL_FUNC) &_alR_Silverman, 1},
    {"_alR_Silverman2", (DL_FUNC) &_alR_Silverman2, 1},
    {"_alR_bw", (DL_FUNC) &_alR_bw, 2},
    {"_alR_GaussInt", (DL_FUNC) &_alR_GaussInt, 5},
    {"_alR_GaussInt2", (DL_FUNC) &_alR_GaussInt2, 5},
    {"_alR_dGPD", (DL_FUNC) &_alR_dGPD, 4},
    {"_alR_pGPD", (DL_FUNC) &_alR_pGPD, 4},
    {"_alR_qGPD", (DL_FUNC) &_alR_qGPD, 4},
    {"_alR_rGPD", (DL_FUNC) &_alR_rGPD, 4},
    {"_alR_GPDInt", (DL_FUNC) &_alR_GPDInt, 6},
    {"_alR_GPDInt2", (DL_FUNC) &_alR_GPDInt2, 6},
    {"_alR_GPDMLE", (DL_FUNC) &_alR_GPDMLE, 1},
    {"_alR_GPDqML", (DL_FUNC) &_alR_GPDqML, 1},
    {"_alR_GPDMSE", (DL_FUNC) &_alR_GPDMSE, 1},
    {"_alR_pkappa4", (DL_FUNC) &_alR_pkappa4, 5},
    {"_alR_dkappa4", (DL_FUNC) &_alR_dkappa4, 5},
    {"_alR_qkappa4", (DL_FUNC) &_alR_qkappa4, 5},
    {"_alR_rkappa4", (DL_FUNC) &_alR_rkappa4, 5},
    {"_alR_dddkappa4", (DL_FUNC) &_alR_dddkappa4, 5},
    {"_alR_kappa4cond", (DL_FUNC) &_alR_kappa4cond, 4},
    {"_alR_kappa4tc", (DL_FUNC) &_alR_kappa4tc, 3},
    {"_alR_kappa4Int", (DL_FUNC) &_alR_kappa4Int, 8},
    {"_alR_kappa4Int2", (DL_FUNC) &_alR_kappa4Int2, 8},
    {"_alR_kappa4IntApprox", (DL_FUNC) &_alR_kappa4IntApprox, 5},
    {"_alR_kappa4IntApprox2", (DL_FUNC) &_alR_kappa4IntApprox2, 5},
    {"_alR_kappa4NLSobj", (DL_FUNC) &_alR_kappa4NLSobj, 5},
    {"_alR_kappa4NLScon", (DL_FUNC) &_alR_kappa4NLScon, 5},
    {"_alR_kappa4NLShin", (DL_FUNC) &_alR_kappa4NLShin, 5},
    {"_alR_kappa4NLSheq", (DL_FUNC) &_alR_kappa4NLSheq, 5},
    {"_alR_kappa4ALobj", (DL_FUNC) &_alR_kappa4ALobj, 6},
    {"_alR_kappa4ALcon", (DL_FUNC) &_alR_kappa4ALcon, 6},
    {"_alR_kappa4ALhin", (DL_FUNC) &_alR_kappa4ALhin, 6},
    {"_alR_kappa4ALheq", (DL_FUNC) &_alR_kappa4ALheq, 6},
    {"_alR_dkdeGauss", (DL_FUNC) &_alR_dkdeGauss, 3},
    {"_alR_pkdeGauss", (DL_FUNC) &_alR_pkdeGauss, 3},
    {"_alR_qkdeGauss", (DL_FUNC) &_alR_qkdeGauss, 3},
    {"_alR_kdeGaussInt", (DL_FUNC) &_alR_kdeGaussInt, 5},
    {"_alR_kdeGaussInt2", (DL_FUNC) &_alR_kdeGaussInt2, 5},
    {"_alR_kdeGaussIntApprox", (DL_FUNC) &_alR_kdeGaussIntApprox, 5},
    {"_alR_kdeGaussIntApprox2", (DL_FUNC) &_alR_kdeGaussIntApprox2, 5},
    {"_alR_momKDE", (DL_FUNC) &_alR_momKDE, 5},
    {"_alR_qlin", (DL_FUNC) &_alR_qlin, 3},
    {"_alR_qsamp", (DL_FUNC) &_alR_qsamp, 2},
    {"_alR_pairSort", (DL_FUNC) &_alR_pairSort, 1},
    {"_alR_lulu", (DL_FUNC) &_alR_lulu, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_alR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
