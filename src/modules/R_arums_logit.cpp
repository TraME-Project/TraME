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

// test functions; delete later
Rcpp::List logit_R::test_1(double x_inp, double y_inp)
{
    double ret = x_inp + y_inp;

    return Rcpp::List::create(Rcpp::Named("val") = ret);
}

arma::mat logit_R::test_2(int x_inp, int y_inp)
{
    arma::mat ret(x_inp, y_inp);
    arma::mat mat_1 = arma::randu(x_inp, y_inp);
    arma::mat mat_2 = arma::randu(x_inp, y_inp);
    ret = mat_1 + mat_2;

    return ret;
}

SEXP logit_R::test_3(int x_inp, int y_inp)
{
    try {
        Rprintf("input 1: %d\n",x_inp);
        Rprintf("input 2: %d\n",y_inp);
        arma::mat ret = test_mat_add_2(x_inp,y_inp);
        int ret_elems = ret.n_elem;
        int ret_rows = ret.n_rows;
        int ret_cols = ret.n_cols;
        //double ret_11 = ret(0,0);
        Rprintf("ret elem: %d\n",ret_elems);
        Rprintf("ret rows: %d\n",ret_rows);
        Rprintf("ret cols: %d\n",ret_cols);
        return R_NilValue;
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
	}
    return R_NilValue;
}

// wrapper function as Rcpp can't handle memory pointers
Rcpp::List logit_R::G_R(arma::vec n, arma::mat U_inp)
{
    arma::mat mu_out;

    double val_out = this->G(n, U_inp, mu_out);

    return Rcpp::List::create(Rcpp::Named("val") = val_out, Rcpp::Named("mu") = mu_out);
}

Rcpp::List logit_R::Gstar_R(arma::vec n, arma::mat mu_inp)
{
    arma::mat U_out;

    double val_out = this->Gstar(n, mu_inp, U_out);

    return Rcpp::List::create(Rcpp::Named("val") = val_out, Rcpp::Named("U") = U_out);
}

Rcpp::List logit_R::Gbar_R(arma::mat U_bar, arma::mat mu_bar, arma::vec n)
{
    arma::mat U_out, mu_out;

    double val_out = this->Gbar(U_bar, mu_bar, n, U_out, mu_out);

    return Rcpp::List::create(Rcpp::Named("val") = val_out, Rcpp::Named("U") = U_out, Rcpp::Named("mu") = mu_out);
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
    double (trame::logit::*G_1)(arma::vec) = &trame::logit::G ;
    Rcpp::List (logit_R::*G_2)(arma::vec, arma::mat) = &logit_R::G_R ;

    double (trame::logit::*Gstar_1)(arma::vec) = &trame::logit::Gstar ;
    Rcpp::List (logit_R::*Gstar_2)(arma::vec, arma::mat) = &logit_R::Gstar_R ;

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
        .method( "test_add", &trame::logit::test_add )
        .method( "test_create", &trame::logit::test_create )
        .method( "test_mat_add", &trame::logit::test_mat_add )
        .method( "test_mat_add_2", &trame::logit::test_mat_add_2 )
        .method( "build", build_1 )
        .method( "build", build_2 )
        .method( "G", G_1 )
        .method( "Gstar", Gstar_1 )
        //.method( "Gstarx", &trame::logit::Gstarx )
    ;

    class_<logit_R>( "logit_R" )
        .derives<trame::logit>( "logit" )
        .default_constructor()

        .method( "test_1", &logit_R::test_1 )
        .method( "test_2", &logit_R::test_2 )
        .method( "test_3", &logit_R::test_3 )
        .method( "G", G_2 )
        .method( "Gstar", Gstar_2 )
        .method( "Gbar", &logit_R::Gbar_R )
        .method( "simul", &logit_R::simul_R )
    ;
}
