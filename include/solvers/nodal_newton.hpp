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
 * Nodal Newton for MFE
 *
 * Keith O'Hara
 * 01/17/2016
 *
 * This version:
 * 02/15/2017
 */

#ifndef _trame_nodal_newton_HPP
#define _trame_nodal_newton_HPP

// internal function
template<typename Tm>
bool nodal_newton_int(const mfe<Tm>& market, arma::mat* mu_out, arma::vec* mu_x0_out, arma::vec* mu_0y_out, arma::mat* U_out, arma::mat* V_out, double* val_out, const double* tol_inp, const int* max_iter_inp);

// wrappers
template<typename Tm>
bool nodal_newton(const mfe<Tm>& market, arma::mat& mu_out);

template<typename Tm>
bool nodal_newton(const mfe<Tm>& market, arma::mat& mu_out, const double& tol_inp);

template<typename Tm>
bool nodal_newton(const mfe<Tm>& market, arma::mat& mu_out, const int& max_iter_inp);

template<typename Tm>
bool nodal_newton(const mfe<Tm>& market, arma::mat& mu_out, const double& tol_inp, const int& max_iter_inp);

template<typename Tm>
bool nodal_newton(const mfe<Tm>& market, arma::mat& mu_out, arma::mat& U_out, arma::mat& V_out);

template<typename Tm>
bool nodal_newton(const mfe<Tm>& market, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::mat& U_out, arma::mat& V_out, double& val_out, const double* tol_inp, const int* max_iter_inp);

// optimization-related functions

bool nodal_newton_optim(arma::vec& init_out_vals, std::function<arma::vec (const arma::vec& vals_inp, void* opt_data)> opt_objfn, void* opt_data);

bool nodal_newton_optim(arma::vec& init_out_vals, std::function<arma::vec (const arma::vec& vals_inp, void* opt_data)> opt_objfn, void* opt_data,
                               std::function<arma::mat (const arma::vec& vals_inp, void* jacob_data)> jacob_objfn, void* jacob_data);

template<typename Tm>
arma::vec nodal_newton_opt_objfn(const arma::vec& vals_inp, void *opt_data);

template<typename Tm>
arma::mat nodal_newton_jacobian(const arma::vec& vals_inp, void *jacob_data);

#include "nodal_newton.tpp"

#endif
