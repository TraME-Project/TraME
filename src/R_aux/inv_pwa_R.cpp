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
 * inverse PWA function
 *
 * Keith O'Hara
 * 09/06/2016
 */

//#define TRAME_RCPP_ARMADILLO

#include "inv_pwa_R.hpp"
using namespace Rcpp;

SEXP inv_pwa_R(SEXP a_R, SEXP B_R, SEXP C_R, SEXP k_R)
{
    try {
        arma::vec a = as<arma::vec>(a_R);
        arma::mat B = as<arma::mat>(B_R);
        arma::mat C = as<arma::mat>(C_R);
        double k = as<double>(k_R);
        
        arma::vec ret = trame::inv_pwa(a,B,C,k);
        //
        return Rcpp::List::create(Rcpp::Named("vals") = ret);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
    return R_NilValue;
}
