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
 * RUSC class module
 *
 * Keith O'Hara
 * 08/08/2016
 */

//#define TRAME_RCPP_ARMADILLO

#include "trame.hpp"
#include "trame_R_modules.hpp"

// wrapper functions to catch errors and handle memory pointers (which Rcpp can't do)
SEXP rusc_R::G_R(arma::vec n)
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

SEXP rusc_R::G_R(arma::vec n, arma::mat U_inp)
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

SEXP rusc_R::Gx_R(arma::mat U_x_inp, int x)
{
    try {
        arma::mat mu_x_out;
        double val_x_out = this->Gx(U_x_inp, mu_x_out, x);
        //
        return Rcpp::List::create(Rcpp::Named("val_x") = val_x_out, Rcpp::Named("mu_x") = mu_x_out);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
	}
    return R_NilValue;
}

SEXP rusc_R::Gstar_R(arma::vec n)
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

SEXP rusc_R::Gstar_R(arma::vec n, arma::mat mu_inp)
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

SEXP rusc_R::Gstarx_R(arma::mat mu_x_inp, int x)
{   
    try {
        arma::mat U_x_out;
        double val_x_out = this->Gstarx(mu_x_inp, U_x_out, x);
        //
        return Rcpp::List::create(Rcpp::Named("val_x") = val_x_out, Rcpp::Named("U_x") = U_x_out);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
	}
    return R_NilValue;
}

SEXP rusc_R::Gbar_R(arma::mat U_bar, arma::mat mu_bar, arma::vec n)
{
    try {
        arma::mat U_out, mu_out;
        double val_out = this->Gbar(U_bar, mu_bar, n, U_out, mu_out);
        //
        return Rcpp::List::create(Rcpp::Named("val") = val_out, Rcpp::Named("U") = U_out, Rcpp::Named("mu") = mu_out);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
	}
    return R_NilValue;
}

empirical_R rusc_R::simul_R(int nbDraws)
{
    trame::empirical emp_obj = this->simul(&nbDraws,NULL);

    empirical_R emp_R_obj = static_cast<empirical_R&>(emp_obj);

    return emp_R_obj;
}

RCPP_EXPOSED_CLASS(empirical_R)
RCPP_EXPOSED_CLASS(rusc_R)

RCPP_MODULE(rusc_module)
{
    using namespace Rcpp ;

    // function overloading requires some trickery
    SEXP (rusc_R::*G_1)(arma::vec) = &rusc_R::G_R ;
    SEXP (rusc_R::*G_2)(arma::vec, arma::mat) = &rusc_R::G_R ;

    SEXP (rusc_R::*Gstar_1)(arma::vec) = &rusc_R::Gstar_R ;
    SEXP (rusc_R::*Gstar_2)(arma::vec, arma::mat) = &rusc_R::Gstar_R ;
  
    // now we can declare the class
    class_<trame::rusc>( "rusc" )
        .default_constructor()

        // basic objects
        .field( "nbX", &trame::rusc::nbX )
        .field( "nbY", &trame::rusc::nbY )

        .field( "nbParams", &trame::rusc::nbParams )
        .field( "outsideOption", &trame::rusc::outsideOption )

        .field( "zeta", &trame::rusc::zeta )

        .field( "U", &trame::rusc::U )
        .field( "mu", &trame::rusc::mu )

        .field( "U_sol", &trame::rusc::U_sol )
        .field( "mu_sol", &trame::rusc::mu_sol )

        // read only objects
        .field_readonly( "aux_ord", &trame::rusc::aux_ord )

        // member functions
        .method( "build", &trame::rusc::build )
    ;

    class_<rusc_R>( "rusc_R" )
        .derives<trame::rusc>( "rusc" )
        .default_constructor()

        .method( "G", G_1 )
        .method( "G", G_2 )
        .method( "Gx", &rusc_R::Gx_R )
        .method( "Gstar", Gstar_1 )
        .method( "Gstar", Gstar_2 )
        .method( "Gstarx", &rusc_R::Gstarx_R )
        .method( "Gbar", &rusc_R::Gbar_R )
        .method( "simul", &rusc_R::simul_R )
    ;
}
