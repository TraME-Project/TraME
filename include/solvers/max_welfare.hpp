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
 * max_welfare for general ARUMs with TU
 *
 * Keith O'Hara
 * 08/16/2016
 *
 * This version:
 * 02/04/2018
 */

#ifndef _trame_max_welfare_HPP
#define _trame_max_welfare_HPP

// internal functions
template<typename Tg, typename Th, typename Tt, typename std::enable_if<!std::is_same<Tt,transfers::tu>::value>::type* = nullptr>
bool max_welfare_int(const dse<Tg,Th,Tt>& market, arma::mat* mu_out, arma::vec* mu_x0_out, arma::vec* mu_0y_out, arma::mat* U_out, arma::mat* V_out, 
                     double* val_out, const double err_tol = 1E-06, const uint_t max_iter = 2000);

template<typename Tg, typename Th, typename Tt, typename std::enable_if<std::is_same<Tt,transfers::tu>::value>::type* = nullptr>
bool max_welfare_int(const dse<Tg,Th,Tt>& market, arma::mat* mu_out, arma::vec* mu_x0_out, arma::vec* mu_0y_out, arma::mat* U_out, arma::mat* V_out, 
                     double* val_out, const double err_tol = 1E-06, const uint_t max_iter = 2000);

// wrappers
template<typename Tg, typename Th, typename Tt>
bool max_welfare(const dse<Tg,Th,Tt>& market, arma::mat& mu_out);

template<typename Tg, typename Th, typename Tt>
bool max_welfare(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, const double err_tol_inp, const uint_t max_iter_inp);

template<typename Tg, typename Th, typename Tt>
bool max_welfare(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, arma::mat& U_out, arma::mat& V_out);

template<typename Tg, typename Th, typename Tt>
bool max_welfare(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::mat& U_out, arma::mat& V_out,
                 double& val_out, const double err_tol_inp, const uint_t max_iter_inp);

// optimization-related functions

bool max_welfare_optim(arma::vec& init_out_vals, std::function<double (const arma::vec& vals_inp, arma::vec* grad, void* opt_data)> opt_objfn, void* opt_data, optim::algo_settings* settings_inp, const uint_t optim_method);

template<typename Tg, typename Th, typename Tt>
double max_welfare_opt_objfn(const arma::vec& vals_inp, arma::vec* grad, void *opt_data);

#include "max_welfare.tpp"

#endif
