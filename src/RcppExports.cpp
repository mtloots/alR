// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// alrKDE
double alrKDE(const arma::vec& beta, const arma::mat& gamma, const arma::vec& aly, const arma::vec& q1, const arma::vec& q2);
RcppExport SEXP alR_alrKDE(SEXP betaSEXP, SEXP gammaSEXP, SEXP alySEXP, SEXP q1SEXP, SEXP q2SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type aly(alySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type q1(q1SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type q2(q2SEXP);
    __result = Rcpp::wrap(alrKDE(beta, gamma, aly, q1, q2));
    return __result;
END_RCPP
}
// dkdeGauss
double dkdeGauss(double x, NumericVector mu, double h);
RcppExport SEXP alR_dkdeGauss(SEXP xSEXP, SEXP muSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    __result = Rcpp::wrap(dkdeGauss(x, mu, h));
    return __result;
END_RCPP
}
// pkdeGauss
List pkdeGauss(double x, NumericVector mu, double h);
RcppExport SEXP alR_pkdeGauss(SEXP xSEXP, SEXP muSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    __result = Rcpp::wrap(pkdeGauss(x, mu, h));
    return __result;
END_RCPP
}
// qkdeGauss
List qkdeGauss(double x, NumericVector mu, double h);
RcppExport SEXP alR_qkdeGauss(SEXP xSEXP, SEXP muSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    __result = Rcpp::wrap(qkdeGauss(x, mu, h));
    return __result;
END_RCPP
}
// kdeGaussInt
List kdeGaussInt(NumericVector mu, double h, double q1, double q2, bool quantile);
RcppExport SEXP alR_kdeGaussInt(SEXP muSEXP, SEXP hSEXP, SEXP q1SEXP, SEXP q2SEXP, SEXP quantileSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type q1(q1SEXP);
    Rcpp::traits::input_parameter< double >::type q2(q2SEXP);
    Rcpp::traits::input_parameter< bool >::type quantile(quantileSEXP);
    __result = Rcpp::wrap(kdeGaussInt(mu, h, q1, q2, quantile));
    return __result;
END_RCPP
}
// kdeGaussInt2
NumericVector kdeGaussInt2(NumericVector mu, double h, NumericVector q1, NumericVector q2, bool quantile);
RcppExport SEXP alR_kdeGaussInt2(SEXP muSEXP, SEXP hSEXP, SEXP q1SEXP, SEXP q2SEXP, SEXP quantileSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q1(q1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q2(q2SEXP);
    Rcpp::traits::input_parameter< bool >::type quantile(quantileSEXP);
    __result = Rcpp::wrap(kdeGaussInt2(mu, h, q1, q2, quantile));
    return __result;
END_RCPP
}
// momKDE
double momKDE(const arma::vec& beta, const arma::mat& gamma, const arma::vec& momy, Function kdeGaussMom);
RcppExport SEXP alR_momKDE(SEXP betaSEXP, SEXP gammaSEXP, SEXP momySEXP, SEXP kdeGaussMomSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type momy(momySEXP);
    Rcpp::traits::input_parameter< Function >::type kdeGaussMom(kdeGaussMomSEXP);
    __result = Rcpp::wrap(momKDE(beta, gamma, momy, kdeGaussMom));
    return __result;
END_RCPP
}
// Silverman
double Silverman(const arma::vec& x);
RcppExport SEXP alR_Silverman(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    __result = Rcpp::wrap(Silverman(x));
    return __result;
END_RCPP
}
