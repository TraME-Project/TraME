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
 * mfe<mmf> class module
 *
 * Keith O'Hara
 * 10/20/2016
 */

#include "trameR.hpp"

RCPP_MODULE(mfe_mmf_module)
{
    using namespace Rcpp ;

    // function overloading requires some trickery
    void (mfe_mmf_R::*build_ETU_1)(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, arma::mat tau_inp, bool need_norm_inp) = &mfe_mmf_R::build_ETU_R ;
    void (mfe_mmf_R::*build_ETU_2)(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, arma::mat tau_inp, double sigma_inp, bool need_norm_inp) = &mfe_mmf_R::build_ETU_R ;

    void (mfe_mmf_R::*build_LTU_1)(arma::vec n_inp, arma::vec m_inp, arma::mat lambda_inp, arma::mat phi_inp, bool need_norm_inp) = &mfe_mmf_R::build_LTU_R ;
    void (mfe_mmf_R::*build_LTU_2)(arma::vec n_inp, arma::vec m_inp, arma::mat lambda_inp, arma::mat phi_inp, double sigma_inp, bool need_norm_inp) = &mfe_mmf_R::build_LTU_R ;

    void (mfe_mmf_R::*build_NTU_1)(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, bool need_norm_inp) = &mfe_mmf_R::build_NTU_R ;
    void (mfe_mmf_R::*build_NTU_2)(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, double sigma_inp, bool need_norm_inp) = &mfe_mmf_R::build_NTU_R ;

    void (mfe_mmf_R::*build_TU_1)(arma::vec n_inp, arma::vec m_inp, arma::mat phi_inp, bool need_norm_inp) = &mfe_mmf_R::build_TU_R ;
    void (mfe_mmf_R::*build_TU_2)(arma::vec n_inp, arma::vec m_inp, arma::mat phi_inp, double sigma_inp, bool need_norm_inp) = &mfe_mmf_R::build_TU_R ;
  
    // now we can declare the class
    class_<trame::mfe<trame::mmf>>( "mfe_mmf" )
        .default_constructor()

        // basic objects
        .field( "ETU", &trame::mfe<trame::mmf>::ETU )
        .field( "LTU", &trame::mfe<trame::mmf>::LTU )
        .field( "NTU", &trame::mfe<trame::mmf>::NTU )
        .field( "TU", &trame::mfe<trame::mmf>::TU )

        .field( "need_norm", &trame::mfe<trame::mmf>::need_norm )
        .field( "outsideOption", &trame::mfe<trame::mmf>::outsideOption )

        .field( "nbX", &trame::mfe<trame::mmf>::nbX )
        .field( "nbY", &trame::mfe<trame::mmf>::nbY )

        // member functions
        .method( "trans", &trame::mfe<trame::mmf>::trans )
    ;

    class_<mfe_mmf_R>( "mfe_mmf_R" )
        .derives<trame::mfe<trame::mmf>>( "mfe_mmf" )
        .default_constructor()

        .method( "build_ETU", build_ETU_1 )
        .method( "build_ETU", build_ETU_2 )
        .method( "build_LTU", build_LTU_1 )
        .method( "build_LTU", build_LTU_2 )
        .method( "build_NTU", build_NTU_1 )
        .method( "build_NTU", build_NTU_2 )
        .method( "build_TU", build_TU_1 )
        .method( "build_TU", build_TU_2 )
        .method( "solve", &mfe_mmf_R::solve_R )
    ;
}

// wrapper functions to catch errors and handle memory pointers
void mfe_mmf_R::build_ETU_R(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, arma::mat tau_inp, bool need_norm_inp)
{
    try {
        this->build_ETU(n_inp, m_inp, alpha_inp, gamma_inp, tau_inp, NULL, need_norm_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
}

void mfe_mmf_R::build_ETU_R(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, arma::mat tau_inp, double sigma_inp, bool need_norm_inp)
{
    try {
        this->build_ETU(n_inp, m_inp, alpha_inp, gamma_inp, tau_inp, &sigma_inp, need_norm_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
}

void mfe_mmf_R::build_LTU_R(arma::vec n_inp, arma::vec m_inp, arma::mat lambda_inp, arma::mat phi_inp, bool need_norm_inp)
{
    try {
        this->build_LTU(n_inp, m_inp, lambda_inp, phi_inp, NULL, need_norm_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
}

void mfe_mmf_R::build_LTU_R(arma::vec n_inp, arma::vec m_inp, arma::mat lambda_inp, arma::mat phi_inp, double sigma_inp, bool need_norm_inp)
{
    try {
        this->build_LTU(n_inp, m_inp, lambda_inp, phi_inp, &sigma_inp, need_norm_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
}

void mfe_mmf_R::build_NTU_R(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, bool need_norm_inp)
{
    try {
        this->build_NTU(n_inp, m_inp, alpha_inp, gamma_inp, NULL, need_norm_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
}

void mfe_mmf_R::build_NTU_R(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, double sigma_inp, bool need_norm_inp)
{
    try {
        this->build_NTU(n_inp, m_inp, alpha_inp, gamma_inp, &sigma_inp, need_norm_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
}

void mfe_mmf_R::build_TU_R(arma::vec n_inp, arma::vec m_inp, arma::mat phi_inp, bool need_norm_inp)
{
    try {
        this->build_TU(n_inp, m_inp, phi_inp, NULL, need_norm_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
}

void mfe_mmf_R::build_TU_R(arma::vec n_inp, arma::vec m_inp, arma::mat phi_inp, double sigma_inp, bool need_norm_inp)
{
    try {
        this->build_TU(n_inp, m_inp, phi_inp, &sigma_inp, need_norm_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
}

SEXP mfe_mmf_R::solve_R()
{
    try {
        arma::mat mu_sol;
        bool success = this->solve(mu_sol);
        //
        return Rcpp::List::create(Rcpp::Named("mu") = mu_sol, Rcpp::Named("success") = success);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
    return R_NilValue;
}
