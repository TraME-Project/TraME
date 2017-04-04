/*################################################################################
  ##
  ##   Copyright (C) 2015 - 2017 the TraME Team:
  ##      Alfred Galichon
  ##      Keith O'Hara
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
 * Generalized logit transform
 *
 * Keith O'Hara
 * 11/28/2014
 */

#ifndef _logit_transform_HPP
#define _logit_transform_HPP

arma::vec logit_trans(const arma::vec& pars, const arma::vec& lower_bounds, const arma::vec& upper_bounds);
arma::vec logit_trans(const arma::vec& pars);
double logit_trans(const double& pars, const double& lower_bounds, const double& upper_bounds);

arma::vec logit_inv_trans(const arma::vec& pars_trans, const arma::vec& lower_bounds, const arma::vec& upper_bounds);
arma::vec logit_inv_trans(const arma::vec& pars_trans);
double logit_inv_trans(const double& pars_trans, const double& lower_bounds, const double& upper_bounds);

#endif
