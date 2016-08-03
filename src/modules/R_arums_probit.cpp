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

#include "../headers/arums_probit.hpp"

RCPP_MODULE(probit_module)
{
    using namespace Rcpp ;

    void (probit::*unifCorrelCovMatrices_1)() = &probit::unifCorrelCovMatrices;
    arma::cube (probit::*unifCorrelCovMatrices_2)(double) = &probit::unifCorrelCovMatrices ;
  
    // now we can declare the class
    class_<probit>( "probit" )
        .default_constructor()

        // basic objects
        .field( "nbX", &probit::nbX )
        .field( "nbY", &probit::nbY )

        .field( "nbParams", &probit::nbParams )
        .field( "aux_nbOptions", &probit::aux_nbOptions )
        .field( "outsideOption", &probit::outsideOption )

        .field( "rho", &probit::rho )

        .field( "Covar", &probit::Covar )

        // read only objects
        //.field_readonly( "", &probit:: )

        // member functions
        .method( "build", &probit::build )
        .method( "unifCorrelCovMatrices", unifCorrelCovMatrices_1 )
        .method( "unifCorrelCovMatrices", unifCorrelCovMatrices_2 )
    ;
}
