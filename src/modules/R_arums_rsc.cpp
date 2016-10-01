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
 * RSC class module
 *
 * Keith O'Hara
 * 08/08/2016
 */

//#define TRAME_RCPP_ARMADILLO

#include "trame.hpp"
#include "trame_R_modules.hpp"

// wrapper functions to catch errors and handle memory pointers (which Rcpp can't do)
SEXP rsc_R::G_R(arma::vec n)
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

SEXP rsc_R::G_R(arma::vec n, arma::mat U_inp)
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

SEXP rsc_R::Gx_R(arma::mat U_x_inp, int x)
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

SEXP rsc_R::Gstar_R(arma::vec n)
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

SEXP rsc_R::Gstar_R(arma::vec n, arma::mat mu_inp)
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

SEXP rsc_R::Gstarx_R(arma::mat mu_x_inp, int x)
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

SEXP rsc_R::Gbar_R(arma::mat U_bar, arma::mat mu_bar, arma::vec n)
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

empirical_R rsc_R::simul_R(int nbDraws)
{
    trame::empirical emp_obj = this->simul(&nbDraws,NULL);

    empirical_R emp_R_obj = static_cast<empirical_R&>(emp_obj);

    return emp_R_obj;
}

RCPP_EXPOSED_CLASS(empirical_R)
RCPP_EXPOSED_CLASS(rsc_R)

RCPP_MODULE(rsc_module)
{
    using namespace Rcpp ;

    // function overloading requires some trickery
    SEXP (rsc_R::*G_1)(arma::vec) = &rsc_R::G_R ;
    SEXP (rsc_R::*G_2)(arma::vec, arma::mat) = &rsc_R::G_R ;

    SEXP (rsc_R::*Gstar_1)(arma::vec) = &rsc_R::Gstar_R ;
    SEXP (rsc_R::*Gstar_2)(arma::vec, arma::mat) = &rsc_R::Gstar_R ;
  
    // now we can declare the class
    class_<trame::rsc>( "rsc" )
        .default_constructor()

        // basic objects
        .field( "nbX", &trame::rsc::nbX )
        .field( "nbY", &trame::rsc::nbY )

        .field( "nbParams", &trame::rsc::nbParams )
        .field( "outsideOption", &trame::rsc::outsideOption )

        .field( "zeta", &trame::rsc::zeta )

        .field( "U", &trame::rsc::U )
        .field( "mu", &trame::rsc::mu )

        .field( "U_sol", &trame::rsc::U_sol )
        .field( "mu_sol", &trame::rsc::mu_sol )

        // read only objects
        .field_readonly( "aux_ord", &trame::rsc::aux_ord )

        .field_readonly( "aux_Influence_lhs", &trame::rsc::aux_Influence_lhs )
        .field_readonly( "aux_Influence_rhs", &trame::rsc::aux_Influence_rhs )

        .field_readonly( "aux_DinvPsigma", &trame::rsc::aux_DinvPsigma )
        .field_readonly( "aux_Psigma", &trame::rsc::aux_Psigma )

        // member functions
        .method( "build", &trame::rsc::build )
        .method( "build_beta", &trame::rsc::build_beta )
    ;

    class_<rsc_R>( "rsc_R" )
        .derives<trame::rsc>( "rsc" )
        .default_constructor()

        .method( "G", G_1 )
        .method( "G", G_2 )
        .method( "Gx", &rsc_R::Gx_R )
        .method( "Gstar", Gstar_1 )
        .method( "Gstar", Gstar_2 )
        .method( "Gstarx", &rsc_R::Gstarx_R )
        .method( "Gbar", &rsc_R::Gbar_R )
        .method( "simul", &rsc_R::simul_R )
    ;
}
