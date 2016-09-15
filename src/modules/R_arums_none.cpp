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
 * none class module
 *
 * Keith O'Hara
 * 08/08/2016
 */


#include <RcppArmadillo.h>

#include "trame.hpp"

// derived class to provide wrappers to some functions
class none_R : public trame::none
{
    public:
        Rcpp::List G_R(arma::vec n, arma::mat U_inp);
        Rcpp::List Gx_R(arma::mat U_x_inp, int x);
        double Gstar_R(arma::vec n);
        double Gstarx_R(arma::mat mu_x_inp, int x);
        Rcpp::List Gbar_R(arma::mat U_bar, arma::mat mu_bar, arma::vec n);
};

// wrapper function as Rcpp can't handle memory pointers
Rcpp::List none_R::G_R(arma::vec n, arma::mat U_inp)
{
    arma::mat mu_out;

    double val_out = this->G(n, U_inp, mu_out);

    return Rcpp::List::create(Rcpp::Named("val") = val_out, Rcpp::Named("mu") = mu_out);
}

Rcpp::List none_R::Gx_R(arma::mat U_x_inp, int x)
{
    arma::mat mu_x_out;

    double val_x_out = this->Gx(U_x_inp, mu_x_out, x);

    return Rcpp::List::create(Rcpp::Named("val_x") = val_x_out, Rcpp::Named("mu_x") = mu_x_out);
}

double none_R::Gstar_R(arma::vec n)
{   
    Rprintf("Gstar not yet defined for no arums case.\n");

    return 0.0;
}

double none_R::Gstarx_R(arma::mat mu_x_inp, int x)
{   
    Rprintf("Gstar not yet defined for no arums case.\n");

    return 0.0;
}

Rcpp::List none_R::Gbar_R(arma::mat U_bar, arma::mat mu_bar, arma::vec n)
{
    arma::mat U_out, mu_out;

    double val_out = this->Gbar(U_bar, mu_bar, n, U_out, mu_out);

    return Rcpp::List::create(Rcpp::Named("val") = val_out, Rcpp::Named("U") = U_out, Rcpp::Named("mu") = mu_out);
}

RCPP_MODULE(none_module)
{
    using namespace Rcpp ;

    // function overloading requires some trickery
    double (trame::none::*G_1)(arma::vec) = &trame::none::G ;
    Rcpp::List (none_R::*G_2)(arma::vec, arma::mat) = &none_R::G_R ;
  
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
        .method( "G", G_1 )
        //.method( "Gstar", &trame::none::Gstar ) // not defined
        //.method( "Gbarx", &trame::none::Gbarx )
    ;

    class_<none_R>( "none_R" )
        .derives<trame::none>( "none" )
        .default_constructor()

        .method( "G", G_2 )
        .method( "Gx", &none_R::Gx_R )
        .method( "Gstar", &none_R::Gstar_R ) // not defined
        .method( "Gstarx", &none_R::Gstar_R ) // not defined
        .method( "Gbar", &none_R::Gbar_R )
    ;
}
