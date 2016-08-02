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

#include "../headers/arums_RUSC.hpp"

// derived class to provide wrappers to some functions
class RUSC_R : public RUSC
{
    public:
        Rcpp::List Gbar_R(arma::mat U_bar, arma::mat mu_bar, arma::vec n);
};

// wrapper function as Rcpp can't handle memory pointers
Rcpp::List RUSC_R::Gbar_R(arma::mat U_bar, arma::mat mu_bar, arma::vec n)
{
    arma::mat U_out, mu_out;

    Gbar(U_bar, mu_bar, n, U_out, mu_out);

    return Rcpp::List::create(Rcpp::Named("U") = U_out, Rcpp::Named("mu") = mu_out);
}

RCPP_MODULE(RUSC_module)
{
    using namespace Rcpp ;

    // function overloading requires some trickery
    double (RUSC::*Gstar_1)(arma::vec) = &RUSC::Gstar ;
  
    // now we can declare the class
    class_<RUSC>( "RUSC" )
        .default_constructor()

        // basic objects
        .field( "nbX", &RUSC::nbX )
        .field( "nbY", &RUSC::nbY )

        .field( "nbParams", &RUSC::nbParams )
        .field( "outsideOption", &RUSC::outsideOption )

        .field( "mu", &RUSC::mu )
        .field( "U", &RUSC::U )

        .field( "zeta", &RUSC::zeta )

        .field( "mu_sol", &RUSC::mu_sol )
        .field( "U_sol", &RUSC::U_sol )

        // read only objects
        .field_readonly( "aux_ord", &RUSC::aux_ord )

        .field_readonly( "aux_Influence_lhs", &RUSC::aux_Influence_lhs )
        .field_readonly( "aux_Influence_rhs", &RUSC::aux_Influence_rhs )

        .field_readonly( "aux_DinvPsigma", &RUSC::aux_DinvPsigma )
        .field_readonly( "aux_Psigma", &RUSC::aux_Psigma )

        // member functions
        .method( "build", &RUSC::build )
        .method( "G", &RUSC::G )
        .method( "Gstar", Gstar_1 )
    ;

    class_<RUSC_R>( "RUSC_R" )
        .derives<RUSC>( "RUSC" )
        .default_constructor()

        .method( "Gbar", Gbar_R )
    ;
}
