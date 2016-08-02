/*################################################################################
  ##
  ##   Copyright (C) 2015 - 2016 Alfred Galichon and the TraME Team
  ##
  ##   This file is part of the R package TraME.
  ##
  ##   The R package TraME is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   The R package TraME is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ##   You should have received a copy of the GNU General Public License
  ##   along with TraME. If not, see <http://www.gnu.org/licenses/>.
  ##
  ################################################################################*/

#include <RcppArmadillo.h>

#include "../headers/arums_empirical.hpp"

// derived class to provide wrappers to some functions
class empirical_R : public empirical
{
    public:
        Rcpp::List Gbar_R(arma::mat U_bar, arma::mat mu_bar, arma::vec n);
};

// wrapper function as Rcpp can't handle memory pointers
Rcpp::List empirical_R::Gbar_R(arma::mat U_bar, arma::mat mu_bar, arma::vec n)
{
    arma::mat U_out, mu_out;

    Gbar(U_bar, mu_bar, n, U_out, mu_out);

    return Rcpp::List::create(Rcpp::Named("U") = U_out, Rcpp::Named("mu") = mu_out);
}

RCPP_MODULE(empirical_module)
{
    using namespace Rcpp ;

    // function overloading requires some trickery
    double (empirical::*Gstar_1)(arma::vec) = &empirical::Gstar ;
    
    // now we can declare the class
    class_<empirical>( "R_empirical" )
        .default_constructor()

        // basic objects
        .field( "nbX", &empirical::nbX )
        .field( "nbY", &empirical::nbY )

        .field( "nbParams", &empirical::nbParams )
        .field( "aux_nbDraws", &empirical::aux_nbDraws )
        .field( "nbOptions", &empirical::nbOptions )

        .field( "xHomogenous", &empirical::xHomogenous )
        .field( "outsideOption", &empirical::outsideOption )

        .field( "atoms", &empirical::atoms )

        // read only objects
        //.field_readonly( "k_Gstar", &empirical::k_Gstar )

        // member functions
        .method( "build", &empirical::build )
        .method( "G", &empirical::G )
        .method( "Gx", &empirical::Gx )
        .method( "Gstar", Gstar_1 )
        .method( "Gstarx", &empirical::Gstarx )
        .method( "Gbar", &empirical::Gbar )
        .method( "Gbarx", &empirical::Gbarx )
    ;
}
