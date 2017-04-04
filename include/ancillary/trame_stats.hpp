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
 * Stats functions
 *
 * Keith O'Hara
 * 08/08/2016
 */

#ifndef _trame_stats_HPP
#define _trame_stats_HPP

double pbeta (double x, double* fn_pars);
arma::vec pbeta (arma::vec x, double* fn_pars);

double qbeta (double x, double* fn_pars);
arma::vec qbeta (arma::vec x, double* fn_pars);

double iqbeta (double x, double* fn_pars);
arma::vec iqbeta (arma::vec x, double* fn_pars);

double dbeta (double x, double* fn_pars);
arma::vec dbeta (arma::vec x, double* fn_pars);

#endif
