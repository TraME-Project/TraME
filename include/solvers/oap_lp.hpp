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
 * Optimal Assignment Problem (OAP) LP solver for TU case
 *
 * Keith O'Hara
 * 08/16/2016
 *
 * This version:
 * 02/04/2018
 */

#ifndef _trame_oap_lp_HPP
#define _trame_oap_lp_HPP

// internal functions
template<typename Tg, typename Th, typename Tt, typename std::enable_if<!std::is_same<Tt,transfers::tu>::value>::type* = nullptr>
bool oap_lp_int(const dse<Tg,Th,Tt>& market, arma::mat* mu_out, arma::vec* mu_x0_out, arma::vec* mu_0y_out, arma::vec* u_out, arma::vec* v_out, 
                double* val_out, arma::mat* residuals_out, const bool x_first = true);

template<typename Tg, typename Th, typename Tt, typename std::enable_if<std::is_same<Tt,transfers::tu>::value>::type* = nullptr>
bool oap_lp_int(const dse<Tg,Th,Tt>& market, arma::mat* mu_out, arma::vec* mu_x0_out, arma::vec* mu_0y_out, arma::vec* u_out, arma::vec* v_out, 
                double* val_out, arma::mat* residuals_out, const bool x_first = true);

// wrappers
template<typename Tg, typename Th, typename Tt>
bool oap_lp(const dse<Tg,Th,Tt>& market, arma::mat& mu_out);

template<typename Tg, typename Th, typename Tt>
bool oap_lp(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, const bool x_first_inp);

template<typename Tg, typename Th, typename Tt>
bool oap_lp(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, arma::mat& residuals_out);

template<typename Tg, typename Th, typename Tt>
bool oap_lp(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, arma::mat& residuals_out, const bool x_first_inp);

template<typename Tg, typename Th, typename Tt>
bool oap_lp(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, arma::vec& u_out, arma::vec& v_out);

template<typename Tg, typename Th, typename Tt>
bool oap_lp(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::vec& u_out, arma::vec& v_out, 
            double& val_out, arma::mat& residuals_out, const bool x_first_inp);

#include "oap_lp.tpp"

#endif
