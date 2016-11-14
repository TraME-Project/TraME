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
 * dse<empirical> class module
 *
 * Keith O'Hara
 * 11/10/2016
 *
 * This version:
 * 11/13/2016
 */

#include "trameR.hpp"

RCPP_EXPOSED_CLASS(empirical_R)
RCPP_EXPOSED_CLASS(logit_R)
RCPP_EXPOSED_CLASS(none_R)
RCPP_EXPOSED_CLASS(probit_R)
RCPP_EXPOSED_CLASS(rsc_R)
RCPP_EXPOSED_CLASS(rusc_R)
RCPP_EXPOSED_CLASS(transfers_R)
RCPP_EXPOSED_CLASS(dse_empirical_R)

RCPP_MODULE(dse_empirical_module)
{
    using namespace Rcpp ;

    // function overloading requires some trickery
    void (dse_empirical_R::*build_LTU_1)(arma::vec n_inp, arma::vec m_inp, arma::mat lambda_inp, arma::mat phi_inp, bool need_norm_inp) = &dse_empirical_R::build_LTU_R ;
    void (dse_empirical_R::*build_LTU_2)(arma::vec n_inp, arma::vec m_inp, arma::mat lambda_inp, arma::mat phi_inp, SEXP arums_G_inp, SEXP arums_H_inp, bool need_norm_inp) = &dse_empirical_R::build_LTU_R ;
    
    void (dse_empirical_R::*build_NTU_1)(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, bool need_norm_inp) = &dse_empirical_R::build_NTU_R ;
    void (dse_empirical_R::*build_NTU_2)(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, SEXP arums_G_inp, SEXP arums_H_inp, bool need_norm_inp) = &dse_empirical_R::build_NTU_R;
    
    void (dse_empirical_R::*build_TU_1)(arma::vec n_inp, arma::vec m_inp, arma::mat phi_inp, bool need_norm_inp) = &dse_empirical_R::build_TU_R ;
    void (dse_empirical_R::*build_TU_2)(arma::vec n_inp, arma::vec m_inp, arma::mat phi_inp, SEXP arums_G_inp, SEXP arums_H_inp, bool need_norm_inp) = &dse_empirical_R::build_TU_R ;

    //SEXP (dse_empirical_R::*solve_R_1)() = &dse_empirical_R::solve_R ;
    SEXP (dse_empirical_R::*solve_R_2)(Rcpp::CharacterVector solver_inp) = &dse_empirical_R::solve_R ;

    // now we can declare the class
    class_<trame::dse<trame::empirical>>( "dse_empirical" )
        .default_constructor()

        // basic objects
        .field( "LTU", &trame::dse<trame::empirical>::LTU )
        .field( "NTU", &trame::dse<trame::empirical>::NTU )
        .field( "TU", &trame::dse<trame::empirical>::TU )

        .field( "need_norm", &trame::dse<trame::empirical>::need_norm )
        .field( "outsideOption", &trame::dse<trame::empirical>::outsideOption )

        .field( "nbX", &trame::dse<trame::empirical>::nbX )
        .field( "nbY", &trame::dse<trame::empirical>::nbY )

        // member functions
        .method( "trans", &trame::dse<trame::empirical>::trans )
    ;

    class_<dse_empirical_R>( "dse_empirical_R" )
        .derives<trame::dse<trame::empirical>>( "dse_empirical" )
        .default_constructor()

        .method( "build_LTU", build_LTU_1 )
        .method( "build_LTU", build_LTU_2 )
        .method( "build_NTU", build_NTU_1 )
        .method( "build_NTU", build_NTU_2 )
        .method( "build_TU", build_TU_1 )
        .method( "build_TU", build_TU_2 )

        //.method( "solve", solve_R_1 )
        .method( "solve", solve_R_2 )

        .method( "get_arums_G", &dse_empirical_R::get_arums_G )
        .method( "set_arums_G", &dse_empirical_R::set_arums_G )
        .method( "get_arums_H", &dse_empirical_R::get_arums_H )
        .method( "set_arums_H", &dse_empirical_R::set_arums_H )
        .method( "set_arums", &dse_empirical_R::set_arums )

        .method( "get_transfers", &dse_empirical_R::get_transfers_R )
        .method( "set_transfers", &dse_empirical_R::set_transfers_R )
    ;
}

// wrapper functions to catch errors and handle memory pointers
void dse_empirical_R::build_LTU_R(arma::vec n_inp, arma::vec m_inp, arma::mat lambda_inp, arma::mat phi_inp, bool need_norm_inp)
{
    try {
        this->build_LTU(n_inp,m_inp,lambda_inp,phi_inp,need_norm_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
}

void dse_empirical_R::build_LTU_R(arma::vec n_inp, arma::vec m_inp, arma::mat lambda_inp, arma::mat phi_inp, SEXP arums_G_inp, SEXP arums_H_inp, bool need_norm_inp)
{
    try {
        if (Rf_inherits(arums_G_inp, "Rcpp_logit_R")) {
            this->build_LTU_R_int<logit_R>(n_inp, m_inp, lambda_inp, phi_inp, arums_G_inp, arums_H_inp, need_norm_inp);
        } 
        else if (Rf_inherits(arums_G_inp, "Rcpp_none_R")) {
            this->build_LTU_R_int<none_R>(n_inp, m_inp, lambda_inp, phi_inp, arums_G_inp, arums_H_inp, need_norm_inp);
        }
        else if (Rf_inherits(arums_G_inp, "Rcpp_probit_R")) {
            this->build_LTU_R_int<probit_R>(n_inp, m_inp, lambda_inp, phi_inp, arums_G_inp, arums_H_inp, need_norm_inp);
        }
        else if (Rf_inherits(arums_G_inp, "Rcpp_rsc_R")) {
            this->build_LTU_R_int<rsc_R>(n_inp, m_inp, lambda_inp, phi_inp, arums_G_inp, arums_H_inp, need_norm_inp);
        }
        else if (Rf_inherits(arums_G_inp, "Rcpp_rusc_R")) {
            this->build_LTU_R_int<rusc_R>(n_inp, m_inp, lambda_inp, phi_inp, arums_G_inp, arums_H_inp, need_norm_inp);
        }
        else {
            ::Rf_error( "trame: unrecognized arums type" );
        }
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
}

void dse_empirical_R::build_NTU_R(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, bool need_norm_inp)
{
    try {
        this->build_NTU(n_inp,m_inp,alpha_inp,gamma_inp,need_norm_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
}

void dse_empirical_R::build_NTU_R(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, SEXP arums_G_inp, SEXP arums_H_inp, bool need_norm_inp)
{
    try {
        if (Rf_inherits(arums_G_inp, "Rcpp_logit_R")) {
            this->build_NTU_R_int<logit_R>(n_inp, m_inp, alpha_inp, gamma_inp, arums_G_inp, arums_H_inp, need_norm_inp);
        } 
        else if (Rf_inherits(arums_G_inp, "Rcpp_none_R")) {
            this->build_NTU_R_int<none_R>(n_inp, m_inp, alpha_inp, gamma_inp, arums_G_inp, arums_H_inp, need_norm_inp);
        }
        else if (Rf_inherits(arums_G_inp, "Rcpp_probit_R")) {
            this->build_NTU_R_int<probit_R>(n_inp, m_inp, alpha_inp, gamma_inp, arums_G_inp, arums_H_inp, need_norm_inp);
        }
        else if (Rf_inherits(arums_G_inp, "Rcpp_rsc_R")) {
            this->build_NTU_R_int<rsc_R>(n_inp, m_inp, alpha_inp, gamma_inp, arums_G_inp, arums_H_inp, need_norm_inp);
        }
        else if (Rf_inherits(arums_G_inp, "Rcpp_rusc_R")) {
            this->build_NTU_R_int<rusc_R>(n_inp, m_inp, alpha_inp, gamma_inp, arums_G_inp, arums_H_inp, need_norm_inp);
        }
        else {
            ::Rf_error( "trame: unrecognized arums type" );
        }
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
}

void dse_empirical_R::build_TU_R(arma::vec n_inp, arma::vec m_inp, arma::mat phi_inp, bool need_norm_inp)
{
    try {
        this->build_TU(n_inp,m_inp,phi_inp,need_norm_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
}

/*void (dse_empirical_R::*build_TU_logit) (arma::vec n_inp, arma::vec m_inp, arma::mat phi_inp, SEXP arums_G_inp, SEXP arums_H_inp, bool need_norm_inp) = &dse_empirical_R::build_TU_R<logit_R>;
void (dse_empirical_R::*build_TU_none)  (arma::vec n_inp, arma::vec m_inp, arma::mat phi_inp, SEXP arums_G_inp, SEXP arums_H_inp, bool need_norm_inp) = &dse_empirical_R::build_TU_R<none_R>;
void (dse_empirical_R::*build_TU_probit)(arma::vec n_inp, arma::vec m_inp, arma::mat phi_inp, SEXP arums_G_inp, SEXP arums_H_inp, bool need_norm_inp) = &dse_empirical_R::build_TU_R<probit_R>;
void (dse_empirical_R::*build_TU_rsc)   (arma::vec n_inp, arma::vec m_inp, arma::mat phi_inp, SEXP arums_G_inp, SEXP arums_H_inp, bool need_norm_inp) = &dse_empirical_R::build_TU_R<rsc_R>;
void (dse_empirical_R::*build_TU_rusc)  (arma::vec n_inp, arma::vec m_inp, arma::mat phi_inp, SEXP arums_G_inp, SEXP arums_H_inp, bool need_norm_inp) = &dse_empirical_R::build_TU_R<rusc_R>;*/

void dse_empirical_R::build_TU_R(arma::vec n_inp, arma::vec m_inp, arma::mat phi_inp, SEXP arums_G_inp, SEXP arums_H_inp, bool need_norm_inp)
{
    try {
        if (Rf_inherits(arums_G_inp, "Rcpp_logit_R")) {
            this->build_TU_R_int<logit_R>(n_inp, m_inp, phi_inp, arums_G_inp, arums_H_inp, need_norm_inp);
        } 
        else if (Rf_inherits(arums_G_inp, "Rcpp_none_R")) {
            this->build_TU_R_int<none_R>(n_inp, m_inp, phi_inp, arums_G_inp, arums_H_inp, need_norm_inp);
        }
        else if (Rf_inherits(arums_G_inp, "Rcpp_probit_R")) {
            this->build_TU_R_int<probit_R>(n_inp, m_inp, phi_inp, arums_G_inp, arums_H_inp, need_norm_inp);
        }
        else if (Rf_inherits(arums_G_inp, "Rcpp_rsc_R")) {
            this->build_TU_R_int<rsc_R>(n_inp, m_inp, phi_inp, arums_G_inp, arums_H_inp, need_norm_inp);
        }
        else if (Rf_inherits(arums_G_inp, "Rcpp_rusc_R")) {
            this->build_TU_R_int<rusc_R>(n_inp, m_inp, phi_inp, arums_G_inp, arums_H_inp, need_norm_inp);
        }
        else {
            ::Rf_error( "trame: unrecognized arums type" );
        }
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
}

SEXP dse_empirical_R::solve_R()
{
    try {
        arma::mat mu_sol;
        bool success = this->solve_R_int(mu_sol, (char*) "darum");
        //
        return Rcpp::List::create(Rcpp::Named("mu") = mu_sol, Rcpp::Named("success") = success);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
    return R_NilValue;
}

SEXP dse_empirical_R::solve_R(Rcpp::CharacterVector solver_inp)
{
    try {
        arma::mat mu_sol;
        //char* solver = solver_inp[0];
        bool success = this->solve_R_int(mu_sol, solver_inp[0]);
        //
        return Rcpp::List::create(Rcpp::Named("mu") = mu_sol, Rcpp::Named("success") = success);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
    return R_NilValue;
}

// this is just a copy of the dse<empirical>::solve function
bool dse_empirical_R::solve_R_int(arma::mat& mu_sol, const char* solver)
{
    bool res = false;
    const char sig = (solver != NULL) ? solver[0] : char(0);
    
    if (solver) { // not NULL
        if (sig=='c') {
            res = cupids_lp(*this,mu_sol);
        }
        if (sig=='d') {
            res = darum(*this,mu_sol);
        }
        if (sig=='e') {
            res = eap_nash(*this,mu_sol);
        }
        if (sig=='j') {
            res = jacobi(*this,mu_sol);
        }
        if (sig=='m') {
            res = max_welfare(*this,mu_sol);
        }
        if (sig=='o') {
            res = oap_lp(*this,mu_sol);
        }
        // default
        if (sig=='n') {
            if (NTU) {
                res = darum(*this,mu_sol);
            } else if (TU) {
                res = max_welfare(*this,mu_sol);
            } else {
                res = jacobi(*this,mu_sol);
            }
        }
    } else {
        if (NTU) {
            res = darum(*this,mu_sol);
        } else if (TU) {
            res = max_welfare(*this,mu_sol);
        } else {
            res = jacobi(*this,mu_sol);
        }
    }
    //
    return res;
}

empirical_R dse_empirical_R::get_arums_G()
{
    empirical_R arums_obj_out = static_cast<empirical_R&>(arums_G);

    return arums_obj_out;
}

void dse_empirical_R::set_arums_G(empirical_R arums_G_inp)
{
    try {
        arums_G = static_cast<trame::empirical&>(arums_G_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
}

empirical_R dse_empirical_R::get_arums_H()
{
    empirical_R arums_obj_out = static_cast<empirical_R&>(arums_H);

    return arums_obj_out;
}

void dse_empirical_R::set_arums_H(empirical_R arums_H_inp)
{
    try {
        arums_H = static_cast<trame::empirical&>(arums_H_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
}

void dse_empirical_R::set_arums(empirical_R arums_G_inp, empirical_R arums_H_inp)
{
    try {
        arums_G = static_cast<trame::empirical&>(arums_G_inp);
        arums_H = static_cast<trame::empirical&>(arums_H_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
}

transfers_R dse_empirical_R::get_transfers_R()
{
    transfers_R trans_obj_out = static_cast<transfers_R&>(trans_obj);

    return trans_obj_out;
}

void dse_empirical_R::set_transfers_R(transfers_R trans_obj_inp)
{
    try {
        trans_obj = static_cast<trame::transfers&>(trans_obj_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
}
