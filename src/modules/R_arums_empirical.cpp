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
 * empirical class module
 *
 * Keith O'Hara
 * 09/06/2016
 */

#include <RcppArmadillo.h>

#include "trame.hpp"

// derived class to provide wrappers to some functions
class empirical_R : public trame::empirical
{
    public:
        Rcpp::List G_R(arma::vec n, arma::mat U_inp);
        Rcpp::List Gbar_R(arma::mat U_bar, arma::mat mu_bar, arma::vec n);
};

// wrapper function as Rcpp can't handle memory pointers
Rcpp::List empirical_R::G_R(arma::vec n, arma::mat U_inp)
{
    arma::mat mu_out;

    double val_out = this->G(n, U_inp, mu_out);

    return Rcpp::List::create(Rcpp::Named("val") = val_out, Rcpp::Named("mu") = mu_out);
}

Rcpp::List empirical_R::Gbar_R(arma::mat U_bar, arma::mat mu_bar, arma::vec n)
{
    arma::mat U_out, mu_out;

    double val_out = this->Gbar(U_bar, mu_bar, n, U_out, mu_out);

    return Rcpp::List::create(Rcpp::Named("val") = val_out, Rcpp::Named("U") = U_out, Rcpp::Named("mu") = mu_out);
}

RCPP_MODULE(empirical_module)
{
    using namespace Rcpp ;

    // function overloading requires some trickery
    double (trame::empirical::*G_1)(arma::vec) = &trame::empirical::G ;
    Rcpp::List (empirical_R::*G_2)(arma::vec, arma::mat) = &empirical_R::G_R ;

    double (trame::empirical::*Gstar_1)(arma::vec) = &trame::empirical::Gstar ;
    
    // now we can declare the class
    class_<trame::empirical>( "empirical" )
        .default_constructor()

        // basic objects
        .field( "nbX", &trame::empirical::nbX )
        .field( "nbY", &trame::empirical::nbY )

        .field( "nbParams", &trame::empirical::nbParams )
        .field( "aux_nbDraws", &trame::empirical::aux_nbDraws )
        .field( "nbOptions", &trame::empirical::nbOptions )

        .field( "xHomogenous", &trame::empirical::xHomogenous )
        .field( "outsideOption", &trame::empirical::outsideOption )

        .field( "atoms", &trame::empirical::atoms )

        // read only objects
        //.field_readonly( "k_Gstar", &empirical::k_Gstar )

        // member functions
        .method( "build", &trame::empirical::build )
        .method( "G", G_1 )
        .method( "Gx", &trame::empirical::Gx )
        .method( "Gstar", Gstar_1 )
        //.method( "Gstarx", &empirical::Gstarx )
    ;

    class_<empirical_R>( "empirical_R" )
        .derives<trame::empirical>( "empirical" )
        .default_constructor()

        .method( "G", G_2 )
        .method( "Gbar", &empirical_R::Gbar_R )
    ;
}
