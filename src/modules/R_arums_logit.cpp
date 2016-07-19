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

#include "headers/arums_logit.hpp"

RCPP_MODULE(logit_module)
{
    using namespace Rcpp ;

    // function overloading requires some trickery
    void (logit::*G_1)(arma::vec) = &logit::G ;
    void (logit::*G_2)(arma::mat,arma::vec) = &logit::G ;

    void (logit::*Gstar_1)(arma::vec) = &logit::Gstar ;
    void (logit::*Gstar_2)(arma::mat,arma::vec) = &logit::Gstar ;
  
    // now we can declare the class
    class_<logit>( "R_logit" )

    .default_constructor()

    // basic objects
    .field( "nbX", &logit::nbX )
    .field( "nbY", &logit::nbY )

    .field( "nbParams", &logit::nbParams )
    .field( "sigma", &logit::sigma )
    .field( "outsideOption", &logit::outsideOption )

    .field( "mu", &logit::mu )
    .field( "U", &logit::U )

    .field( "mu_sol", &logit::mu )
    .field( "U_sol", &logit::U )

    // read only objects
    //.field_readonly( "", &logit:: )

    // member functions
    .method( "build", &logit::build )
    .method( "G", G_1 )
    .method( "G", G_2 )
    .method( "Gstar", Gstar_1 )
    .method( "Gstar", GStar_2 )
    .method( "Gstarx", &logit::Gstarx )
    .method( "Gbar", &logit::Gbar )
    .method( "Gbarx", &logit::Gbarx )
    ;
}
