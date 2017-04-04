/*################################################################################
  ##
  ##   Copyright (C) 2015 - 2017 the TraME Team:
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
 * auxiliary functions for equilibrium solvers
 *
 * Keith O'Hara
 * 08/23/2016
 *
 * This version:
 * 03/22/2017
 */

#ifndef _trame_aux_solvers_HPP
#define _trame_aux_solvers_HPP

template<typename Tg, typename Th, typename Tm>
struct trame_market_opt_data {
    dse<Tg,Th,Tm> market;
};

template<typename Tm>
struct trame_mfe_opt_data {
    mfe<Tm> market;
};

//

int build_disaggregate_epsilon(arma::vec n, const arums::empirical& arums_emp_inp, arma::mat& epsilon_iy, arma::mat& epsilon0_i, arma::mat& I_ix);

template<typename Tg, typename Th, typename Tm>
arma::mat w_upper_bound(const dse<Tg,Th,Tm>& market);

#include "aux_solvers.tpp"

#endif
