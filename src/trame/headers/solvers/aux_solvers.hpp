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
 * auxiliary functions for equilibrium solvers
 *
 * Keith O'Hara
 * 08/23/2016
 */

#ifndef _trame_aux_solvers_HPP
#define _trame_aux_solvers_HPP

arma::mat u_from_vs(transfers trans_obj, arma::mat v, double* tol_inp, arma::mat* subdiff);
arma::mat v_from_us(transfers trans_obj, arma::mat u, double* tol_inp, arma::mat* subdiff);
arma::mat update_v(transfers trans_obj, arma::mat v, arma::vec n, arma::vec m, bool xFirst);

template<typename Ta>
arma::mat w_upper_bound(dse<Ta> market);

#include "aux_solvers.tpp"

#endif
