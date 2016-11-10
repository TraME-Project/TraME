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
 * transfers class module
 *
 * Keith O'Hara
 * 11/08/2016
 */

//#define TRAME_RCPP_ARMADILLO

#include "trame.hpp"
#include "../R_modules/trame_R_modules.hpp"
#include "trame_R_markets.hpp"

// wrapper functions to catch errors and handle memory pointers (which Rcpp can't do)
void transfers_R::build_ETU_R(arma::mat alpha_ETU, arma::mat gamma_ETU, arma::mat tau_ETU)
{
    try {
        this->build_ETU(alpha_ETU,gamma_ETU,tau_ETU);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
	}
}

void transfers_R::build_LTU_R(arma::mat lambda_LTU, arma::mat phi_LTU)
{
    try {
        this->build_LTU(lambda_LTU,phi_LTU);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
	}
}

void transfers_R::build_NTU_R(arma::mat alpha_NTU, arma::mat gamma_NTU)
{
    try {
        this->build_NTU(alpha_NTU,gamma_NTU);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
	}
}

void transfers_R::build_TU_R(arma::mat phi_TU)
{
    try {
        this->build_TU(phi_TU);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
	}
}

void transfers_R::trans_R()
{
    try {
        this->trans();
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
	}
}

SEXP transfers_R::Psi_R(arma::mat U, arma::mat V)
{
    try {
        arma::mat psi_out = this->Psi(U,V);
        //
        return Rcpp::List::create(Rcpp::Named("psi") = psi_out);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
	}
    return R_NilValue;
}

SEXP transfers_R::Psi_R(arma::mat U, arma::mat V, Rcpp::IntegerVector x_ind, Rcpp::IntegerVector y_ind)
{
    try {
        int x_ind_size = x_ind.size();
        int y_ind_size = y_ind.size();

        arma::mat psi_out;
        //
        // default case to mirror NULL
        if (x_ind_size == 0 && y_ind_size == 0) {
            psi_out = this->Psi(U,V);
        }
        //
        // correct for zero indexing vs R indexing
        arma::uvec x_ind_uvec, y_ind_uvec;

        if (x_ind_size != 0) {
            x_ind_uvec = Rcpp::as<arma::uvec>(x_ind) - 1;
        }
        if (y_ind_size != 0) {
            y_ind_uvec = Rcpp::as<arma::uvec>(y_ind) - 1;
        }
        //
        if (x_ind_size != 0 && y_ind_size == 0) {
            psi_out = this->Psi(U,V,&x_ind_uvec,NULL);
        } else if (x_ind_size == 0 && y_ind_size != 0) {
            psi_out = this->Psi(U,V,NULL,&y_ind_uvec);
        } else {
            psi_out = this->Psi(U,V,&x_ind_uvec,&y_ind_uvec);
        }
        //
        return Rcpp::List::create(Rcpp::Named("psi") = psi_out);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
	}
    return R_NilValue;
}

RCPP_MODULE(transfers_module)
{
    using namespace Rcpp ;

    // function overloading requires some trickery
    SEXP (transfers_R::*Psi_1)(arma::mat, arma::mat) = &transfers_R::Psi_R ;
    SEXP (transfers_R::*Psi_2)(arma::mat, arma::mat, Rcpp::IntegerVector, Rcpp::IntegerVector) = &transfers_R::Psi_R ;

    // now we can declare the class
    class_<trame::transfers>( "transfers" )
        .default_constructor()

        // basic objects
        .field( "ETU", &trame::transfers::ETU )
        .field( "LTU", &trame::transfers::LTU )
        .field( "NTU", &trame::transfers::NTU )
        .field( "TU", &trame::transfers::TU )

        .field( "nbX", &trame::transfers::nbX )
        .field( "nbY", &trame::transfers::nbY )
        .field( "nbParams", &trame::transfers::nbParams )

        .field( "phi", &trame::transfers::phi )
        
        .field( "alpha", &trame::transfers::alpha )
        .field( "gamma", &trame::transfers::gamma )
        .field( "lambda", &trame::transfers::lambda )
        .field( "tau", &trame::transfers::tau )

        // read only objects
        //.field_readonly( "", &trame::transfers:: )

        // member functions
    ;

    class_<transfers_R>( "transfers_R" )
        .derives<trame::transfers>( "transfers" )
        .default_constructor()

        .method( "build_ETU", &transfers_R::build_ETU_R )
        .method( "build_LTU", &transfers_R::build_LTU_R )
        .method( "build_NTU", &transfers_R::build_NTU_R )
        .method( "build_TU", &transfers_R::build_TU_R )

        .method( "trans", &transfers_R::trans_R )

        .method( "Psi", Psi_1 )
        .method( "Psi", Psi_2 )
    ;
}
