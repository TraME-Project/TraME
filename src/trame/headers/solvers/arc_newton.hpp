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
 * Arc Newton for TU_general
 *
 * Keith O'Hara
 * 01/17/2016
 */

#ifndef _trame_arc_newton_HPP
#define _trame_arc_newton_HPP

// internal function
template<typename Ta>
bool arc_newton_int(const dse<Ta>& market, arma::mat* mu_out, arma::vec* mu_x0_out, arma::vec* mu_0y_out, arma::mat* U_out, arma::mat* V_out, double* val_out, const double* tol_inp, const int* max_iter_inp);

// wrappers
template<typename Ta>
bool arc_newton(const dse<Ta>& market, arma::mat& mu_out);

template<typename Ta>
bool arc_newton(const dse<Ta>& market, arma::mat& mu_out, const double& tol_inp);

template<typename Ta>
bool arc_newton(const dse<Ta>& market, arma::mat& mu_out, const int& max_iter_inp);

template<typename Ta>
bool arc_newton(const dse<Ta>& market, arma::mat& mu_out, const double& tol_inp, const int& max_iter_inp);

template<typename Ta>
bool arc_newton(const dse<Ta>& market, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::mat& U_out, arma::mat& V_out, double& val_out, const double* tol_inp, const int* max_iter_inp);

// optimization-related functions

bool arc_newton_optim(arma::vec& init_out_vals, std::function<arma::vec (const arma::vec& vals_inp, void* opt_data)> opt_objfn, void* opt_data);

template<typename Ta>
arma::vec arc_newton_opt_objfn(const arma::vec& vals_inp, void *opt_data);

#include "arc_newton.tpp"

#endif
