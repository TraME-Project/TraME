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
 * Deferred Acceptance (DA) algorithm for general ARUMs with NTU
 *
 * Keith O'Hara
 * 08/16/2016
 *
 * This version:
 * 02/04/2018
 */

#ifndef _trame_darum_HPP
#define _trame_darum_HPP

// internal functions
template<typename Tg, typename Th, typename Tt, typename std::enable_if<!std::is_same<Tt,transfers::ntu>::value>::type* = nullptr>
bool darum_int(const dse<Tg,Th,Tt>& market, arma::mat* mu_out, arma::vec* mu_x0_out, arma::vec* mu_0y_out, arma::mat* U_out, arma::mat* V_out, 
               const double err_tol = 1E-10, const uint_t max_iter = 5000);

template<typename Tg, typename Th, typename Tt, typename std::enable_if<std::is_same<Tt,transfers::ntu>::value>::type* = nullptr>
bool darum_int(const dse<Tg,Th,Tt>& market, arma::mat* mu_out, arma::vec* mu_x0_out, arma::vec* mu_0y_out, arma::mat* U_out, arma::mat* V_out, 
               const double err_tol = 1E-10, const uint_t max_iter = 5000);

// wrappers
template<typename Tg, typename Th, typename Tt>
bool darum(const dse<Tg,Th,Tt>& market, arma::mat& mu_out);

template<typename Tg, typename Th, typename Tt>
bool darum(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, const double tol_inp, const uint_t max_iter_inp);

template<typename Tg, typename Th, typename Tt>
bool darum(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, arma::mat& U_out, arma::mat& V_out);

template<typename Tg, typename Th, typename Tt>
bool darum(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::mat& U_out, arma::mat& V_out,
           const double tol_inp, const uint_t max_iter_inp);

#include "darum.tpp"

#endif
