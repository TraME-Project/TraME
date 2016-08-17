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
 * ipfp for logit
 *
 * Keith O'Hara
 * 08/16/2016
 */

#ifndef _trame_ipfp_HPP
#define _trame_ipfp_HPP

bool ipfp(mfe market, bool xFirst, double* tol_inp, arma::vec* by_start, arma::mat& mu, arma::vec& mux0, arma::vec& mu0y, arma::mat& U, arma::mat& V, arma::vec& u, arma::vec& v);

#endif
