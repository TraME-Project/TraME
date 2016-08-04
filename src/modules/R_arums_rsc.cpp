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

#include "../nlopt/generic_nlopt.hpp"

#include "../aux/trame_aux.hpp"
//#include "../aux/trame_structs.hpp"

#include "../headers/arums_empirical.hpp"

#include "../headers/arums_RSC.hpp"

// derived class to provide wrappers to some functions
class RSC_R : public RSC
{
    public:
        Rcpp::List Gbar_R(arma::mat U_bar, arma::mat mu_bar, arma::vec n);
};

// wrapper function as Rcpp can't handle memory pointers
Rcpp::List RSC_R::Gbar_R(arma::mat U_bar, arma::mat mu_bar, arma::vec n)
{
    arma::mat U_out, mu_out;

    Gbar(U_bar, mu_bar, n, U_out, mu_out);

    return Rcpp::List::create(Rcpp::Named("U") = U_out, Rcpp::Named("mu") = mu_out);
}

RCPP_MODULE(RSC_module)
{
    using namespace Rcpp ;

    // function overloading requires some trickery
    double (RSC::*Gstar_1)(arma::vec) = &RSC::Gstar ;
  
    // now we can declare the class
    class_<RSC>( "RSC" )
        .default_constructor()

        // basic objects
        .field( "nbX", &RSC::nbX )
        .field( "nbY", &RSC::nbY )

        .field( "nbParams", &RSC::nbParams )
        .field( "outsideOption", &RSC::outsideOption )

        .field( "mu", &RSC::mu )
        .field( "U", &RSC::U )

        .field( "zeta", &RSC::zeta )

        .field( "mu_sol", &RSC::mu_sol )
        .field( "U_sol", &RSC::U_sol )

        // read only objects
        .field_readonly( "aux_ord", &RSC::aux_ord )

        .field_readonly( "aux_Influence_lhs", &RSC::aux_Influence_lhs )
        .field_readonly( "aux_Influence_rhs", &RSC::aux_Influence_rhs )

        .field_readonly( "aux_DinvPsigma", &RSC::aux_DinvPsigma )
        .field_readonly( "aux_Psigma", &RSC::aux_Psigma )

        // member functions
        .method( "build", &RSC::build )
        .method( "build_beta", &RSC::build_beta )
        .method( "G", &RSC::G )
        .method( "Gstar", Gstar_1 )
    ;

    class_<RSC_R>( "RSC_R" )
        .derives<RSC>( "RSC" )
        .default_constructor()

        .method( "Gbar", &RSC_R::Gbar_R )
    ;
}
