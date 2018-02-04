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
 * Jacobi solver
 *
 * Keith O'Hara
 * 08/25/2016
 *
 * This version:
 * 02/04/2018
 */

#ifndef _trame_jacobi_HPP
#define _trame_jacobi_HPP

// internal
template<typename Tg, typename Th, typename Tt>
bool jacobi_int(const dse<Tg,Th,Tt>& market, const arma::mat* w_low_inp, const arma::mat* w_up_inp, arma::mat* mu_out, arma::vec* mu_x0_out, arma::vec* mu_0y_out, arma::mat* U_out, arma::mat* V_out, 
                const double err_tol = 1E-04, const uint_t max_iter = 5000);

// wrappers
template<typename Tg, typename Th, typename Tt>
bool jacobi(const dse<Tg,Th,Tt>& market, arma::mat& mu_out);

template<typename Tg, typename Th, typename Tt>
bool jacobi(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, const double err_tol_inp, const uint_t max_iter_inp);

template<typename Tg, typename Th, typename Tt>
bool jacobi(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, arma::mat& U_out, arma::mat& V_out);

template<typename Tg, typename Th, typename Tt>
bool jacobi(const dse<Tg,Th,Tt>& market, const arma::mat& w_low_inp, const arma::mat& w_up_inp, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::mat& U_out, arma::mat& V_out, 
            const double* err_tol_inp, const uint_t* max_iter_inp);

// zeroin function

template<typename Tg, typename Th, typename Tt>
double jacobi_zeroin_fn(double z, void* opt_data);

template<typename Tg, typename Th, typename Tt>
struct trame_jacobi_zeroin_data {
    uint_t x_ind;
    uint_t y_ind;

    arma::mat U;
    arma::mat V;

    dse<Tg,Th,Tt> market_obj;
};

#include "jacobi.tpp"

#endif
