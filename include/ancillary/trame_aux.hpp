/*################################################################################
  ##
  ##   Copyright (C) 2015-2018 the TraME Team:
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
 * auxiliary functions
 *
 * Keith O'Hara
 * 08/08/2016
 *
 * This version:
 * 02/03/2018
 */

#ifndef _trame_aux_HPP
#define _trame_aux_HPP

double lse(const arma::mat& X);

arma::uvec which_max(const arma::mat& X, const uint_t which_dim);
arma::vec unit_vec(const uint_t j, const size_t n);
arma::uvec uvec_linspace(const int a, const int b);
int* uword_to_int(const arma::uword* var_inp, const size_t n_elem);

arma::mat byrow(const arma::mat& X, const size_t n_rows, const size_t n_cols);

arma::mat elem_add(const arma::mat& mat_1, const arma::mat& mat_2);
arma::mat elem_sub(const arma::mat& mat_1, const arma::mat& mat_2);
arma::mat elem_prod(const arma::mat& mat_1, const arma::mat& mat_2);
arma::mat elem_div(const arma::mat& mat_1, const arma::mat& mat_2);

double elem_min(const arma::mat& mat_1);
arma::mat elem_min(const arma::mat& mat_1, const arma::mat& mat_2);
arma::mat elem_min(const arma::mat& mat_1, const double comp_val);
arma::mat elem_min(const double comp_val, const arma::mat& mat_1);
double elem_min(const double comp_val_1, const double comp_val_2);

double elem_max(const arma::mat& mat_1);
arma::mat elem_max(const arma::mat& mat_1, const arma::mat& mat_2);
arma::mat elem_max(const arma::mat& mat_1, const double comp_val);
arma::mat elem_max(const double comp_val, const arma::mat& mat_1);
double elem_max(const double comp_val_1, const double comp_val_2);

arma::mat cube_sum(const arma::cube& cube_inp, const uint_t which_dim);
arma::mat cube_to_mat(const arma::cube& cube_inp);

#include "trame_aux.ipp"

#endif
