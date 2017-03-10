/*################################################################################
  ##
  ##   Copyright (C) 2015 - 2017 the TraME Team:
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
 * auxiliary functions
 *
 * Keith O'Hara
 * 08/08/2016
 */

#ifndef _trame_aux_HPP
#define _trame_aux_HPP

inline arma::uvec which_max(const arma::mat& X, int which_dim);
inline arma::vec unit_vec(int j, int n);
inline arma::uvec uvec_linspace(int a, int b);

inline arma::mat elem_add(const arma::mat& mat_1, const arma::mat& mat_2);
inline arma::mat elem_sub(const arma::mat& mat_1, const arma::mat& mat_2);
inline arma::mat elem_prod(const arma::mat& mat_1, const arma::mat& mat_2);
inline arma::mat elem_div(const arma::mat& mat_1, const arma::mat& mat_2);
inline double elem_min(const arma::mat& mat_1);
inline double elem_max(const arma::mat& mat_1);
inline arma::mat elem_max(const arma::mat& mat_1, const double& comp_val);

inline arma::mat cube_sum(const arma::cube& cube_inp, int which_dim);

#include "trame_aux.ipp"

#endif
