/*################################################################################
  ##
  ##   Copyright (C) 2015 - 2017 the TraME Team:
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
 * none class module
 *
 * Keith O'Hara
 * 08/08/2016
 */

#include "trameR.hpp"

RCPP_EXPOSED_CLASS(empirical_R)
RCPP_EXPOSED_CLASS(none_R)

RCPP_MODULE(none_module)
{
    using namespace Rcpp ;

    // function overloading requires some trickery
    SEXP (none_R::*G_1)(arma::vec) = &none_R::G_R ;
    SEXP (none_R::*G_2)(arma::vec, arma::mat) = &none_R::G_R ;

    SEXP (none_R::*Gstar_1)(arma::vec) = &none_R::Gstar_R ;
    SEXP (none_R::*Gstar_2)(arma::vec, arma::mat) = &none_R::Gstar_R ;
  
    // now we can declare the class
    class_<trame::none>( "none" )
        .default_constructor()

        // basic objects
        .field( "nbX", &trame::none::nbX )
        .field( "nbY", &trame::none::nbY )

        .field( "nbParams", &trame::none::nbParams )

        .field( "mu", &trame::none::mu )
        .field( "U", &trame::none::U )

        .field( "mu_sol", &trame::none::mu )
        .field( "U_sol", &trame::none::U )

        // read only objects
        //.field_readonly( "", &trame::none:: )

        // member functions
        .method( "build", &trame::none::build )
    ;

    class_<none_R>( "none_R" )
        .derives<trame::none>( "none" )
        .default_constructor()

        .method( "G", G_1 )
        .method( "G", G_2 )
        .method( "Gx", &none_R::Gx_R )
        .method( "Gstar", Gstar_1 )
        .method( "Gstar", Gstar_2 )
        .method( "Gstarx", &none_R::Gstarx_R )
        .method( "Gbar", &none_R::Gbar_R )
        .method( "simul", &none_R::simul_R )
    ;
}

// wrapper functions to catch errors and handle memory pointers
SEXP none_R::G_R(arma::vec n)
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

SEXP none_R::G_R(arma::vec n, arma::mat U_inp)
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

SEXP none_R::Gx_R(arma::mat U_x_inp, int x)
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

SEXP none_R::Gstar_R(arma::vec n)
{   
    Rprintf("Gstar not yet defined for no arums case.\n");
    //
    return R_NilValue;
}

SEXP none_R::Gstar_R(arma::vec n, arma::mat mu_inp)
{   
    Rprintf("Gstar not yet defined for no arums case.\n");
    //
    return R_NilValue;
}

SEXP none_R::Gstarx_R(arma::mat mu_x_inp, int x)
{   
    Rprintf("Gstar not yet defined for no arums case.\n");
    //
    return R_NilValue;
}

SEXP none_R::Gbar_R(arma::mat U_bar, arma::mat mu_bar, arma::vec n)
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

empirical_R none_R::simul_R(int nbDraws)
{
    trame::empirical emp_obj = this->simul(&nbDraws,NULL);

    empirical_R emp_R_obj = static_cast<empirical_R&>(emp_obj);

    return emp_R_obj;
}
