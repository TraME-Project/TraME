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

//#define TRAME_RCPP_ARMADILLO

#include "trame.hpp"
#include "trame_R_modules.hpp"

// wrapper function as Rcpp can't handle memory pointers
Rcpp::List empirical_R::G_R(arma::vec n, arma::mat U_inp)
{
    arma::mat mu_out;

    double val_out = this->G(n, U_inp, mu_out);

    return Rcpp::List::create(Rcpp::Named("val") = val_out, Rcpp::Named("mu") = mu_out);
}

Rcpp::List empirical_R::Gx_R(arma::mat U_x_inp, int x)
{
    arma::mat mu_x_out;

    double val_x_out = this->Gx(U_x_inp, mu_x_out, x);

    return Rcpp::List::create(Rcpp::Named("val_x") = val_x_out, Rcpp::Named("mu_x") = mu_x_out);
}

Rcpp::List empirical_R::Gstar_R(arma::vec n, arma::mat mu_inp)
{   
    arma::mat U_out;

    double val_out = this->Gstar(n, mu_inp, U_out);

    return Rcpp::List::create(Rcpp::Named("val") = val_out, Rcpp::Named("U") = U_out);
}

Rcpp::List empirical_R::Gstarx_R(arma::mat mu_x_inp, int x)
{   
    arma::mat U_x_out;

    double val_x_out = this->Gstarx(mu_x_inp, U_x_out, x);

    return Rcpp::List::create(Rcpp::Named("val_x") = val_x_out, Rcpp::Named("U_x") = U_x_out);
}

Rcpp::List empirical_R::Gbar_R(arma::mat U_bar, arma::mat mu_bar, arma::vec n)
{
    arma::mat U_out, mu_out;

    double val_out = this->Gbar(U_bar, mu_bar, n, U_out, mu_out);

    return Rcpp::List::create(Rcpp::Named("val") = val_out, Rcpp::Named("U") = U_out, Rcpp::Named("mu") = mu_out);
}

RCPP_EXPOSED_CLASS(empirical_R)

RCPP_MODULE(empirical_module)
{
    using namespace Rcpp ;

    // function overloading requires some trickery
    double (trame::empirical::*G_1)(arma::vec) = &trame::empirical::G ;
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

        .field( "U", &trame::empirical::U )
        .field( "mu", &trame::empirical::mu )

        .field( "U_sol", &trame::empirical::U_sol )
        .field( "mu_sol", &trame::empirical::mu_sol )

        // read only objects
        //.field_readonly( "k_Gstar", &empirical::k_Gstar )

        // member functions
        .method( "build", &trame::empirical::build )
        .method( "G", G_1 )
        .method( "Gstar", Gstar_1 )
    ;

    class_<empirical_R>( "empirical_R" )
        .derives<trame::empirical>( "empirical" )
        .default_constructor()

        .method( "G", &empirical_R::G_R )
        .method( "Gx", &empirical_R::Gx_R )
        .method( "Gstar", &empirical_R::Gstar_R )
        .method( "Gstarx", &empirical_R::Gstarx_R )
        .method( "Gbar", &empirical_R::Gbar_R )
    ;
}
