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
 * dse<logit> class module
 *
 * Keith O'Hara
 * 10/23/2016
 */

//#define TRAME_RCPP_ARMADILLO

#include "trame.hpp"
#include "../R_modules/trame_R_modules.hpp"
#include "trame_R_markets.hpp"

void dse_logit_R::build_LTU_R(arma::vec n_inp, arma::vec m_inp, arma::mat lambda_inp, arma::mat phi_inp, bool need_norm_inp)
{
    try {
        int nbX = n_inp.n_elem;
        int nbY = m_inp.n_elem;

        trame::logit logit_1(nbX,nbY), logit_2(nbY,nbX);
        this->build_LTU(n_inp,m_inp,lambda_inp,phi_inp,logit_1,logit_2,need_norm_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
}

void dse_logit_R::build_LTU_R(arma::vec n_inp, arma::vec m_inp, arma::mat lambda_inp, arma::mat phi_inp, logit_R arums_G_inp, logit_R arums_H_inp, bool need_norm_inp)
{
    try {
        trame::logit logit_1 = static_cast<trame::logit&>(arums_G_inp);
        trame::logit logit_2 = static_cast<trame::logit&>(arums_H_inp);

        this->build_LTU(n_inp,m_inp,lambda_inp,phi_inp,logit_1,logit_2,need_norm_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
}

void dse_logit_R::build_NTU_R(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, bool need_norm_inp)
{
    try {
        int nbX = n_inp.n_elem;
        int nbY = m_inp.n_elem;

        trame::logit logit_1(nbX,nbY), logit_2(nbY,nbX);
        this->build_NTU(n_inp,m_inp,alpha_inp,gamma_inp,logit_1,logit_2,need_norm_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
}

void dse_logit_R::build_NTU_R(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, logit_R arums_G_inp, logit_R arums_H_inp, bool need_norm_inp)
{
    try {
        trame::logit logit_1 = static_cast<trame::logit&>(arums_G_inp);
        trame::logit logit_2 = static_cast<trame::logit&>(arums_H_inp);

        this->build_NTU(n_inp,m_inp,alpha_inp,gamma_inp,logit_1,logit_2,need_norm_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
}

void dse_logit_R::build_TU_R(arma::vec n_inp, arma::vec m_inp, arma::mat phi_inp, bool need_norm_inp)
{
    try {
        int nbX = n_inp.n_elem;
        int nbY = m_inp.n_elem;

        trame::logit logit_1(nbX,nbY), logit_2(nbY,nbX);
        this->build_TU(n_inp,m_inp,phi_inp,logit_1,logit_2,need_norm_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
}

void dse_logit_R::build_TU_R(arma::vec n_inp, arma::vec m_inp, arma::mat phi_inp, logit_R arums_G_inp, logit_R arums_H_inp, bool need_norm_inp)
{
    try {
        trame::logit logit_1 = static_cast<trame::logit&>(arums_G_inp);
        trame::logit logit_2 = static_cast<trame::logit&>(arums_H_inp);

        this->build_TU(n_inp,m_inp,phi_inp,logit_1,logit_2,need_norm_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
}

SEXP dse_logit_R::solve_R()
{
    try {
        arma::mat mu_sol;
        bool success = this->solve(mu_sol, (char*) "darum");
        //
        return Rcpp::List::create(Rcpp::Named("mu") = mu_sol, Rcpp::Named("success") = success);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
    return R_NilValue;
}

SEXP dse_logit_R::solve_R(Rcpp::CharacterVector solver_inp)
{
    try {
        arma::mat mu_sol;
        //char* solver = solver_inp[0];
        bool success = this->solve(mu_sol, solver_inp[0]);
        //
        return Rcpp::List::create(Rcpp::Named("mu") = mu_sol, Rcpp::Named("success") = success);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
    return R_NilValue;
}

logit_R dse_logit_R::get_arums_G()
{
    logit_R arums_obj_out = static_cast<logit_R&>(arums_G);

    return arums_obj_out;
}

void dse_logit_R::set_arums_G(logit_R arums_G_inp)
{
    try {
        arums_G = static_cast<trame::logit&>(arums_G_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
}

logit_R dse_logit_R::get_arums_H()
{
    logit_R arums_obj_out = static_cast<logit_R&>(arums_H);

    return arums_obj_out;
}

void dse_logit_R::set_arums_H(logit_R arums_H_inp)
{
    try {
        arums_H = static_cast<trame::logit&>(arums_H_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
}

void dse_logit_R::set_arums(logit_R arums_G_inp, logit_R arums_H_inp)
{
    try {
        arums_G = static_cast<trame::logit&>(arums_G_inp);
        arums_H = static_cast<trame::logit&>(arums_H_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
}

transfers_R dse_logit_R::get_transfers_R()
{
    transfers_R trans_obj_out = static_cast<transfers_R&>(trans_obj);

    return trans_obj_out;
}

void dse_logit_R::set_transfers_R(transfers_R trans_obj_inp)
{
    try {
        trans_obj = static_cast<trame::transfers&>(trans_obj_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
}

RCPP_EXPOSED_CLASS(logit_R)
RCPP_EXPOSED_CLASS(transfers_R)
RCPP_EXPOSED_CLASS(dse_logit_R)

RCPP_MODULE(dse_logit_module)
{
    using namespace Rcpp ;

    // function overloading requires some trickery
    void (dse_logit_R::*build_LTU_1)(arma::vec n_inp, arma::vec m_inp, arma::mat lambda_inp, arma::mat phi_inp, bool need_norm_inp) = &dse_logit_R::build_LTU_R ;
    void (dse_logit_R::*build_LTU_2)(arma::vec n_inp, arma::vec m_inp, arma::mat lambda_inp, arma::mat phi_inp, logit_R arums_G_inp, logit_R arums_H_inp, bool need_norm_inp) = &dse_logit_R::build_LTU_R ;
    
    void (dse_logit_R::*build_NTU_1)(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, bool need_norm_inp) = &dse_logit_R::build_NTU_R ;
    void (dse_logit_R::*build_NTU_2)(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, logit_R arums_G_inp, logit_R arums_H_inp, bool need_norm_inp) = &dse_logit_R::build_NTU_R ;
    
    void (dse_logit_R::*build_TU_1)(arma::vec n_inp, arma::vec m_inp, arma::mat phi_inp, bool need_norm_inp) = &dse_logit_R::build_TU_R ;
    void (dse_logit_R::*build_TU_2)(arma::vec n_inp, arma::vec m_inp, arma::mat phi_inp, logit_R arums_G_inp, logit_R arums_H_inp, bool need_norm_inp) = &dse_logit_R::build_TU_R ;
    
    //SEXP (dse_logit_R::*solve_R_1)() = &dse_logit_R::solve_R ;
    SEXP (dse_logit_R::*solve_R_2)(Rcpp::CharacterVector solver_inp) = &dse_logit_R::solve_R ;

    // now we can declare the class
    class_<trame::dse<trame::logit>>( "dse_logit" )
        .default_constructor()

        // basic objects
        .field( "LTU", &trame::dse<trame::logit>::LTU )
        .field( "NTU", &trame::dse<trame::logit>::NTU )
        .field( "TU", &trame::dse<trame::logit>::TU )

        .field( "need_norm", &trame::dse<trame::logit>::need_norm )
        .field( "outsideOption", &trame::dse<trame::logit>::outsideOption )

        .field( "nbX", &trame::dse<trame::logit>::nbX )
        .field( "nbY", &trame::dse<trame::logit>::nbY )

        // member functions
        .method( "trans", &trame::dse<trame::logit>::trans )
    ;

    class_<dse_logit_R>( "dse_logit_R" )
        .derives<trame::dse<trame::logit>>( "dse_logit" )
        .default_constructor()

        .method( "build_LTU", build_LTU_1 )
        .method( "build_LTU", build_LTU_2 )
        .method( "build_NTU", build_NTU_1 )
        .method( "build_NTU", build_NTU_2 )
        .method( "build_TU", build_TU_1 )
        .method( "build_TU", build_TU_2 )

        //.method( "solve", solve_R_1 )
        .method( "solve", solve_R_2 )

        .method( "get_arums_G", &dse_logit_R::get_arums_G )
        .method( "set_arums_G", &dse_logit_R::set_arums_G )
        .method( "get_arums_H", &dse_logit_R::get_arums_H )
        .method( "set_arums_H", &dse_logit_R::set_arums_H )
        .method( "set_arums", &dse_logit_R::set_arums )

        .method( "get_transfers", &dse_logit_R::get_transfers_R )
        .method( "set_transfers", &dse_logit_R::set_transfers_R )
    ;
}
