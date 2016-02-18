// invPWA_R.cpp
// Keith O'Hara
// 02/18/2016

#include "invPWA_R.h"
using namespace Rcpp;

SEXP invPWA_R(SEXP a_R, SEXP B_R, SEXP C_R)
{
    try {
        arma::vec a = as<arma::vec>(a_R);
        arma::mat B = as<arma::mat>(B_R);
        arma::mat C = as<arma::mat>(C_R);
        
        arma::vec vals = invPWA(a,B,C);
        
        return Rcpp::List::create(Rcpp::Named("vals") = vals);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "C++ exception (unknown reason)" );
    }
    return R_NilValue;
}