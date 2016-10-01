/*################################################################################
  ##
  ##   Copyright (C) 2015 - 2016 the TraME Team:
  ##      Alfred Galichon
  ##      Keith O'Hara
  ##      Simon Weber
  ##
  ##   This file is part of TraME.
  ##
  ##   TraME is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   TraME is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ##   You should have received a copy of the GNU General Public License
  ##   along with TraME. If not, see <http://www.gnu.org/licenses/>.
  ##
  ################################################################################*/

/*
 * logit class module
 *
 * Keith O'Hara
 * 09/06/2016
 */

//#define TRAME_RCPP_ARMADILLO

#include "trame.hpp"
#include "trame_R_modules.hpp"

// wrapper functions to catch errors and handle memory pointers (which Rcpp can't do)
SEXP logit_R::G_R(arma::vec n)
{
    try {
        double val_out = this->G(n);
        //
        return Rcpp::List::create(Rcpp::Named("val") = val_out);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
	}
    return R_NilValue;
}

SEXP logit_R::G_R(arma::vec n, arma::mat U_inp)
{
    try {
        arma::mat mu_out;
        double val_out = this->G(n, U_inp, mu_out);
        //
        return Rcpp::List::create(Rcpp::Named("val") = val_out, Rcpp::Named("mu") = mu_out);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
	}
    return R_NilValue;
}

SEXP logit_R::Gstar_R(arma::vec n)
{
    try {
        double val_out = this->Gstar(n);
        //
        return Rcpp::List::create(Rcpp::Named("val") = val_out);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
	}
    return R_NilValue;
}

SEXP logit_R::Gstar_R(arma::vec n, arma::mat mu_inp)
{
    try {
        arma::mat U_out;
        double val_out = this->Gstar(n, mu_inp, U_out);
        //
        return Rcpp::List::create(Rcpp::Named("val") = val_out, Rcpp::Named("U") = U_out);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
	}
    return R_NilValue;
}

SEXP logit_R::Gbar_R(arma::mat U_bar, arma::mat mu_bar, arma::vec n)
{
    try {
        arma::mat U_out, mu_out;
        double val_out = this->Gbar(U_bar, mu_bar, n, U_out, mu_out);

        return Rcpp::List::create(Rcpp::Named("val") = val_out, Rcpp::Named("U") = U_out, Rcpp::Named("mu") = mu_out);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
	}
    return R_NilValue;
}

empirical_R logit_R::simul_R(int nbDraws)
{
    trame::empirical emp_obj = this->simul(&nbDraws,NULL);

    empirical_R emp_R_obj = static_cast<empirical_R&>(emp_obj);

    return emp_R_obj;
}

RCPP_EXPOSED_CLASS(empirical_R)
RCPP_EXPOSED_CLASS(logit_R)

RCPP_MODULE(logit_module)
{
    using namespace Rcpp ;

    // function overloading requires some trickery
    void (trame::logit::*build_1)(int, int) = &trame::logit::build ;
    void (trame::logit::*build_2)(int, int, double, bool) = &trame::logit::build ;

    SEXP (logit_R::*G_1)(arma::vec) = &logit_R::G_R ;
    SEXP (logit_R::*G_2)(arma::vec, arma::mat) = &logit_R::G_R ;

    SEXP (logit_R::*Gstar_1)(arma::vec) = &logit_R::Gstar_R ;
    SEXP (logit_R::*Gstar_2)(arma::vec, arma::mat) = &logit_R::Gstar_R ;

    // now we can declare the class
    class_<trame::logit>( "logit" )
        .default_constructor()

        // basic objects
        .field( "nbX", &trame::logit::nbX )
        .field( "nbY", &trame::logit::nbY )

        .field( "nbParams", &trame::logit::nbParams )
        .field( "sigma", &trame::logit::sigma )
        .field( "outsideOption", &trame::logit::outsideOption )

        .field( "U", &trame::logit::U )
        .field( "mu", &trame::logit::mu )

        .field( "U_sol", &trame::logit::U_sol )
        .field( "mu_sol", &trame::logit::mu_sol )

        // read only objects
        //.field_readonly( "", &trame::logit:: )

        // member functions
        .method( "build", build_1 )
        .method( "build", build_2 )
    ;

    class_<logit_R>( "logit_R" )
        .derives<trame::logit>( "logit" )
        .default_constructor()

        .method( "G", G_1 )
        .method( "G", G_2 )
        .method( "Gstar", Gstar_1 )
        .method( "Gstar", Gstar_2 )
        .method( "Gbar", &logit_R::Gbar_R )
        .method( "simul", &logit_R::simul_R )
    ;
}
