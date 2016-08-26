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
 * OAP LP
 *
 * Keith O'Hara
 * 08/16/2016
 */

#ifndef _trame_oap_lp_HPP
#define _trame_oap_lp_HPP

template<typename Ta>
bool oap_lp(dse<Ta> market, bool xFirst, arma::mat& mu, arma::vec& mux0, arma::vec& mu0y, arma::vec& u, arma::vec& v, double& val, arma::mat& residuals);

#include "oap_lp.tpp"

#endif

