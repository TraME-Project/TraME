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

#include "../aux/trame_aux.hpp"
#include "../aux/trame_structs.hpp"

#include "../headers/arums_empirical.hpp"

#include "../headers/arums_none.hpp"

// derived class to provide wrappers to some functions
class none_R : public none
{
    public:
        Rcpp::List Gbar_R(arma::mat U_bar, arma::mat mu_bar, arma::vec n);
};

// wrapper function as Rcpp can't handle memory pointers
Rcpp::List none_R::Gbar_R(arma::mat U_bar, arma::mat mu_bar, arma::vec n)
{
    arma::mat U_out, mu_out;

    Gbar(U_bar, mu_bar, n, U_out, mu_out);

    return Rcpp::List::create(Rcpp::Named("U") = U_out, Rcpp::Named("mu") = mu_out);
}

RCPP_MODULE(none_module)
{
    using namespace Rcpp ;
  
    class_<none>( "none" )
        .default_constructor()

        // basic objects
        .field( "nbX", &none::nbX )
        .field( "nbY", &none::nbY )

        .field( "nbParams", &none::nbParams )

        .field( "mu", &none::mu )
        .field( "U", &none::U )

        .field( "mu_sol", &none::mu )
        .field( "U_sol", &none::U )

        // read only objects
        //.field_readonly( "", &none:: )

        // member functions
        .method( "build", &none::build )
        .method( "G", &none::G )
        .method( "Gx", &none::Gx )
        //.method( "Gstar", &none::Gstar ) // not defined
        .method( "Gbarx", &none::Gbarx )
    ;

    class_<none_R>( "none_R" )
        .derives<none>( "none" )
        .default_constructor()

        .method( "Gbar", &none_R::Gbar_R )
    ;
}
