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
 * Jacobi solver
 *
 * Keith O'Hara
 * 08/25/2016
 *
 * This version:
 * 03/22/2017
 */

#ifndef _trame_jacobi_HPP
#define _trame_jacobi_HPP

template<typename Tg, typename Th, typename Tm>
struct trame_jacobi_zeroin_data {
    int x_ind;
    int y_ind;

    arma::vec n;
    arma::vec m;

    arma::mat U;
    arma::mat V;

    Tg arums_G;
    Th arums_H;

    Tm trans_obj;
};

// internal
template<typename Tg, typename Th, typename Tm>
bool jacobi_int(const dse<Tg,Th,Tm>& market, const arma::mat* w_low_inp, const arma::mat* w_up_inp, arma::mat* mu_out, arma::vec* mu_x0_out, arma::vec* mu_0y_out, arma::mat* U_out, arma::mat* V_out, const double* tol_inp, const int* max_iter_inp);

// wrappers
template<typename Tg, typename Th, typename Tm>
bool jacobi(const dse<Tg,Th,Tm>& market, arma::mat& mu_out);

template<typename Tg, typename Th, typename Tm>
bool jacobi(const dse<Tg,Th,Tm>& market, arma::mat& mu_out, const double& tol_inp);

template<typename Tg, typename Th, typename Tm>
bool jacobi(const dse<Tg,Th,Tm>& market, arma::mat& mu_out, const int& max_iter_inp);

template<typename Tg, typename Th, typename Tm>
bool jacobi(const dse<Tg,Th,Tm>& market, arma::mat& mu_out, const double& tol_inp, const int& max_iter_inp);

template<typename Tg, typename Th, typename Tm>
bool jacobi(const dse<Tg,Th,Tm>& market, arma::mat& mu_out, arma::mat& U_out, arma::mat& V_out);

template<typename Tg, typename Th, typename Tm>
bool jacobi(const dse<Tg,Th,Tm>& market, const arma::mat& w_low_inp, const arma::mat& w_up_inp, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::mat& U_out, arma::mat& V_out, const double* tol_inp, const int* max_iter_inp);

// zeroin function

template<typename Tg, typename Th, typename Tm>
double jacobi_zeroin_fn(double z, void* opt_data);

#include "jacobi.tpp"

#endif
