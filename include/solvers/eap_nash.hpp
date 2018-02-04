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
 * Equilibrium Assignment Problem (EAP) - Nash
 *
 * Keith O'Hara
 * 08/16/2016
 *
 * This version:
 * 02/04/2018
 */

#ifndef _trame_eap_nash_HPP
#define _trame_eap_nash_HPP

// internal function
template<typename Tg, typename Th, typename Tt>
bool eap_nash_int(const dse<Tg,Th,Tt>& market, arma::mat* mu_out, arma::vec* mu_x0_out, arma::vec* mu_0y_out, arma::mat* u_out, arma::mat* v_out,
                  const bool x_first = true, const double err_tol = 1E-10, const uint_t max_iter = 5000);

// wrappers
template<typename Tg, typename Th, typename Tt>
bool eap_nash(const dse<Tg,Th,Tt>& market, arma::mat& mu_out);

template<typename Tg, typename Th, typename Tt>
bool eap_nash(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, const bool x_first_inp);

template<typename Tg, typename Th, typename Tt>
bool eap_nash(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, const double err_tol_inp, const uint_t max_iter_inp);

template<typename Tg, typename Th, typename Tt>
bool eap_nash(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, const bool x_first_inp, const double err_tol_inp, const uint_t max_iter_inp);

template<typename Tg, typename Th, typename Tt>
bool eap_nash(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, arma::mat& u_out, arma::mat& v_out);

template<typename Tg, typename Th, typename Tt>
bool eap_nash(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::mat& u_out, arma::mat& v_out,
              const bool x_first_inp, const double err_tol_inp, const uint_t max_iter_inp);

// internal functions

template<typename Tt>
arma::mat u_from_vs(const Tt& trans_obj, const arma::mat& v, const double* tol_inp, arma::mat* subdiff);

template<typename Tt>
arma::mat v_from_us(const Tt& trans_obj, const arma::mat& u, const double* tol_inp, arma::mat* subdiff);

template<typename Tt>
arma::mat update_v(const Tt& trans_obj, const arma::mat& v, const arma::vec& n, const arma::vec& m, const bool x_first);

#include "eap_nash.tpp"

#endif
