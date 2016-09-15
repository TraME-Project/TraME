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

#include <RcppArmadillo.h>

#include "trame.hpp"

// derived class to provide wrappers to some functions
class logit_R : public trame::logit
{
    public:
        Rcpp::List G_R(arma::vec n, arma::mat U_inp);
        Rcpp::List Gbar_R(arma::mat U_bar, arma::mat mu_bar, arma::vec n);
};

// wrapper function as Rcpp can't handle memory pointers
Rcpp::List logit_R::G_R(arma::vec n, arma::mat U_inp)
{
    arma::mat mu_out;

    double val_out = this->G(n, U_inp, mu_out);

    return Rcpp::List::create(Rcpp::Named("val") = val_out, Rcpp::Named("mu") = mu_out);
}

Rcpp::List logit_R::Gbar_R(arma::mat U_bar, arma::mat mu_bar, arma::vec n)
{
    arma::mat U_out, mu_out;

    double val_out = this->Gbar(U_bar, mu_bar, n, U_out, mu_out);

    return Rcpp::List::create(Rcpp::Named("val") = val_out, Rcpp::Named("U") = U_out, Rcpp::Named("mu") = mu_out);
}

RCPP_MODULE(logit_module)
{
    using namespace Rcpp ;

    // function overloading requires some trickery
    double (trame::logit::*G_1)(arma::vec) = &trame::logit::G ;
    Rcpp::List (logit_R::*G_2)(arma::vec, arma::mat) = &logit_R::G_R ;

    double (trame::logit::*Gstar_1)(arma::vec) = &trame::logit::Gstar ;
  
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
        .method( "build", &trame::logit::build )
        .method( "G", G_1 )
        .method( "Gstar", Gstar_1 )
        //.method( "Gstarx", &trame::logit::Gstarx )
    ;

    class_<logit_R>( "logit_R" )
        .derives<trame::logit>( "logit" )
        .default_constructor()

        .method( "G", G_2 )
        .method( "Gbar", &logit_R::Gbar_R )
    ;
}
